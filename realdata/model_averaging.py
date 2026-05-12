from pathlib import Path

import cvxpy as cp
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from sklearn.linear_model import LassoCV, LinearRegression, lasso_path
from sklearn.metrics import mean_squared_error, r2_score
from sklearn.preprocessing import StandardScaler


SCRIPT_DIR = Path(__file__).resolve().parent
DATA_PATH = SCRIPT_DIR / "cirrhosis.csv"
OUTPUT_PATH = SCRIPT_DIR / "real_res.pdf"


def _scalar(value):
    """Safely convert a length-1 array-like object to a Python float."""
    return float(np.asarray(value).reshape(-1)[0])


def lasso_alpha_bisection(x, y, initial_alpha=1e10):
    alpha0 = initial_alpha
    res_alpha_list = [alpha0]
    _, coef_path, _ = lasso_path(x, y, alphas=[alpha0])
    coef_path = coef_path[0, :, :].T
    nonzero_index = np.where(coef_path[0] != 0)[0]
    current_val_num = len(nonzero_index)
    if current_val_num in (0, 1):
        alpha = alpha0 / 2
    else:
        raise ValueError("Initial alpha does not start from a sparse enough model.")
    while alpha > 1e-10 and len(res_alpha_list) <= x.shape[1]:
        _, coef_path, _ = lasso_path(x, y, alphas=[alpha])
        coef_path = coef_path[0, :, :].T
        nonzero_index = np.where(coef_path[0] != 0)[0]
        if len(nonzero_index) <= current_val_num:
            alpha0 = alpha
            alpha = alpha / 2
        elif len(nonzero_index) == current_val_num + 1:
            res_alpha_list.append(alpha)
            current_val_num += 1
            alpha0 = alpha
            alpha = alpha / 2
        else:
            alpha = alpha + (alpha0 - alpha) / 2
    return res_alpha_list


def fit_model_average(x_sc, y_sc):
    use_alpha = lasso_alpha_bisection(x_sc, y_sc)
    alpha, coef_path, _ = lasso_path(x_sc, y_sc, alphas=use_alpha)
    coef_path = coef_path[0, :, :].T

    current_val_num = 1
    model_index_list = []
    for i in range(len(coef_path)):
        nonzero_index = np.where(coef_path[i] != 0)[0]
        if len(nonzero_index) == current_val_num:
            current_val_num += 1
            model_index_list.append(nonzero_index)

    sample_size = y_sc.shape[0]
    model_coef_list = []
    model_pred = []
    for idx in model_index_list:
        lr_model = LinearRegression()
        lr_model.fit(x_sc[:, idx], y_sc)
        model_coef_list.append(lr_model.coef_)
        model_pred.append((x_sc[:, idx] @ lr_model.coef_.T).reshape(-1))

    model_pred = np.array(model_pred)
    sigma2 = mean_squared_error(y_sc.reshape(-1), model_pred[-1, :]) * sample_size / (
        sample_size - len(model_index_list[-1])
    )
    k = np.array([len(model_index_list[i]) for i in range(len(model_index_list))])
    w = cp.Variable(len(model_index_list), pos=True)
    objective = cp.sum_squares(w @ model_pred - y_sc.reshape(-1)) + 2 * sigma2 * k @ w
    problem = cp.Problem(cp.Minimize(objective), [cp.sum(w) == 1])
    problem.solve()
    weights = w.value

    averaged_pred = np.sum(np.array([weights[i] * model_pred[i, :] for i in range(len(model_index_list))]), axis=0)
    freedom = np.where(weights > 1e-8)[0][-1] + 1
    sigma_ma = mean_squared_error(y_sc.reshape(-1), averaged_pred) * sample_size / (sample_size - freedom)

    return {
        "weights": weights,
        "model_index_list": model_index_list,
        "model_coef_list": model_coef_list,
        "sigma_ma": sigma_ma,
        "averaged_pred": averaged_pred,
    }


data = pd.read_csv(DATA_PATH)
data["Sex"] = data["Sex"].apply(lambda x: 0 if x == "F" else 1)
data["Ascites"] = data["Ascites"].apply(lambda x: 0 if x == "Y" else 1)
data["Hepatomegaly"] = data["Hepatomegaly"].apply(lambda x: 0 if x == "Y" else 1)
data["Spiders"] = data["Spiders"].apply(lambda x: 0 if x == "Y" else 1)
data["N_Days"] = data["N_Days"] / 365
data["Age"] = data["Age"] / 365
data.dropna(inplace=True)

edema1 = []
edema2 = []
for i in range(len(data)):
    if data["Edema"].iloc[i] == "Y":
        edema1.append(0)
        edema2.append(0)
    if data["Edema"].iloc[i] == "N":
        edema1.append(1)
        edema2.append(0)
    if data["Edema"].iloc[i] == "S":
        edema1.append(0)
        edema2.append(1)
data["edema1"] = edema1
data["edema2"] = edema2

add_years = []
average_stat = [25.1 * 0.9, 23.6 * 0.8, 21.1 * 0.7, 13.8 * 0.6]
for i in range(len(data)):
    if data["Status"].iloc[i] == "D":
        add_years.append(0)
    if data["Status"].iloc[i] == "CL":
        add_years.append(1)
    if data["Status"].iloc[i] == "C":
        if data["Stage"].iloc[i] == 1:
            add_years.append(average_stat[0])
        elif data["Stage"].iloc[i] == 2:
            add_years.append(average_stat[1])
        elif data["Stage"].iloc[i] == 3:
            add_years.append(average_stat[2])
        else:
            add_years.append(average_stat[3])
data["N_Days"] = data["N_Days"] + add_years

feature_list = [
    "ID",
    "N_Days",
    "Status",
    "Drug",
    "Age",
    "Sex",
    "Ascites",
    "Spiders",
    "Bilirubin",
    "Cholesterol",
    "Albumin",
    "Copper",
    "Alk_Phos",
    "SGOT",
    "Tryglicerides",
    "Platelets",
    "Prothrombin",
    "Stage",
    "edema1",
    "edema2",
]

data_t = data[data["Drug"] == "D-penicillamine"].dropna()
data_c = data[data["Drug"] == "Placebo"].dropna()

x_t = data_t.loc[:, feature_list]
y_t = data_t.loc[:, ["N_Days"]]
x_t_sc = StandardScaler(with_std=False).fit_transform(x_t.iloc[:, 4:])
y_t_sc = StandardScaler(with_std=False).fit_transform(y_t)

x_c = data_c.loc[:, feature_list]
y_c = data_c.loc[:, ["N_Days"]]
x_c_sc = StandardScaler(with_std=False).fit_transform(x_c.iloc[:, 4:])
y_c_sc = StandardScaler(with_std=False).fit_transform(y_c)

t_fit = fit_model_average(x_t_sc, y_t_sc)
c_fit = fit_model_average(x_c_sc, y_c_sc)

data_all = pd.concat((data_c, data_t))
x_all = data_all.loc[:, feature_list]
x_mu = np.mean(x_all.iloc[:, 4:], axis=0).values.reshape(1, -1)
x_t_mu = np.mean(x_t.iloc[:, 4:], axis=0).values.reshape(1, -1)
x_c_mu = np.mean(x_c.iloc[:, 4:], axis=0).values.reshape(1, -1)
y_t_group_mean = _scalar(y_t.to_numpy().mean())
y_c_group_mean = _scalar(y_c.to_numpy().mean())

model_pred_mu_t = []
model_pred_mu_c = []
for i in range(len(t_fit["model_index_list"])):
    model_pred_mu_t.append(
        _scalar((x_t_mu - x_mu)[:, t_fit["model_index_list"][i]] @ t_fit["model_coef_list"][i].T)
    )
for i in range(len(c_fit["model_index_list"])):
    model_pred_mu_c.append(
        _scalar((x_c_mu - x_mu)[:, c_fit["model_index_list"][i]] @ c_fit["model_coef_list"][i].T)
    )

y_t_mu = y_t_group_mean
y_c_mu = y_c_group_mean
for i in range(len(t_fit["model_index_list"])):
    y_t_mu -= _scalar(t_fit["weights"][i] * model_pred_mu_t[i])
for i in range(len(c_fit["model_index_list"])):
    y_c_mu -= _scalar(c_fit["weights"][i] * model_pred_mu_c[i])
ate_ma = float(y_t_mu - y_c_mu)

sample_size_t = y_t_sc.shape[0]
sample_size_c = y_c_sc.shape[0]
sigma_ma = float((c_fit["sigma_ma"] / sample_size_c + t_fit["sigma_ma"] / sample_size_t) ** 0.5)
sigma_y = float((np.var(data_t["N_Days"]) / sample_size_t + np.var(data_c["N_Days"]) / sample_size_c) ** 0.5)

lasso_t = LassoCV()
lasso_t.fit(x_t_sc, y_t_sc.ravel())
lasso_c = LassoCV()
lasso_c.fit(x_c_sc, y_c_sc.ravel())
sigma_lasso_t = mean_squared_error(y_t_sc.reshape(-1), x_t_sc @ lasso_t.coef_) * sample_size_t / (
    sample_size_t - np.sum(lasso_t.coef_ != 0)
)
sigma_lasso_c = mean_squared_error(y_c_sc.reshape(-1), x_c_sc @ lasso_c.coef_) * sample_size_c / (
    sample_size_c - np.sum(lasso_c.coef_ != 0)
)
y_t_mu_lasso = y_t_group_mean - _scalar((x_t_mu - x_mu) @ lasso_t.coef_)
y_c_mu_lasso = y_c_group_mean - _scalar((x_c_mu - x_mu) @ lasso_c.coef_)
ate_lasso = float(y_t_mu_lasso - y_c_mu_lasso)
sigma_lasso = float((sigma_lasso_c / sample_size_c + sigma_lasso_t / sample_size_t) ** 0.5)

plt.figure(figsize=(4.7, 3.4))
my_xticks = ["Unadjusted", "Lasso", "Model Averaging"]
plt.xticks([0, 1, 2], my_xticks)
plt.xlim(-0.5, 2.5)
plt.vlines(2, ate_ma - sigma_ma * 1.96, ate_ma + sigma_ma * 1.96)
plt.vlines(1, ate_lasso - sigma_lasso * 1.96, ate_lasso + sigma_lasso * 1.96)
ate_unadj = float(np.mean(y_t) - np.mean(y_c))
plt.vlines(0, ate_unadj - sigma_y * 1.96, ate_unadj + sigma_y * 1.96)
plt.hlines(0, -0.5, 2.5, color="black", linestyle="dashed")
plt.plot(0, ate_unadj, marker="o", color="red")
plt.plot(1, ate_lasso, marker="o", color="red")
plt.plot(2, ate_ma, marker="o", color="red")
plt.hlines(ate_unadj - sigma_y * 1.96, -0.05, 0.05)
plt.hlines(ate_unadj + sigma_y * 1.96, -0.05, 0.05)
plt.hlines(ate_lasso - sigma_lasso * 1.96, 1 - 0.05, 1 + 0.05)
plt.hlines(ate_lasso + sigma_lasso * 1.96, 1 - 0.05, 1 + 0.05)
plt.hlines(ate_ma - sigma_ma * 1.96, 2 - 0.05, 2 + 0.05)
plt.hlines(ate_ma + sigma_ma * 1.96, 2 - 0.05, 2 + 0.05)
plt.ylabel("ATE estimates")
plt.tight_layout()
plt.savefig(OUTPUT_PATH, bbox_inches="tight")
plt.close()

print(f"Saved {OUTPUT_PATH}")
print(f"Treatment-group R^2: {r2_score(y_t_sc, t_fit['weights'] @ np.array([(x_t_sc[:, idx] @ coef.T).reshape(-1) for idx, coef in zip(t_fit['model_index_list'], t_fit['model_coef_list'])])):.3f}")
print(f"Control-group R^2: {r2_score(y_c_sc, c_fit['weights'] @ np.array([(x_c_sc[:, idx] @ coef.T).reshape(-1) for idx, coef in zip(c_fit['model_index_list'], c_fit['model_coef_list'])])):.3f}")
