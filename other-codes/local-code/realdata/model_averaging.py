#%%
import numpy as np
import pandas as pd
from sklearn.linear_model import LinearRegression, Lasso, LassoCV, MultiTaskLassoCV,  lasso_path
from sklearn.preprocessing import StandardScaler, MinMaxScaler
from sklearn.metrics import mean_squared_error
from sklearn.metrics import r2_score
import matplotlib.pyplot as plt
import cvxpy as cp
# %%
data = pd.read_csv("/Users/bytedance/pythonfile/RCT/cirrhosis.csv")
# data = data[data["Status"]=="D"]


data["Sex"] = data["Sex"].apply(lambda x: 0 if x=="F" else 1)
data["Ascites"] = data["Ascites"].apply(lambda x: 0 if x=="Y" else 1)
data["Hepatomegaly"] = data["Hepatomegaly"].apply(lambda x: 0 if x=="Y" else 1)
data["Spiders"] = data["Spiders"].apply(lambda x: 0 if x=="Y" else 1)
data["Hepatomegaly"] = data["Hepatomegaly"].apply(lambda x: 0 if x=="Y" else 1)

data["N_Days"] = data["N_Days"]/365
data["Age"] = data["Age"]/365

data.dropna(inplace=True)
#%%
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
average_stat = [25.1*0.9, 23.6*0.8, 21.1*0.7, 13.8*0.6]
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
# %%
col_names = data.columns
feature_list = ['ID', 'N_Days', 'Status', 'Drug', 'Age', 'Sex', 'Ascites',
        'Spiders', 'Bilirubin', 'Cholesterol',
        'Albumin', 'Copper', 'Alk_Phos', 'SGOT', 'Tryglicerides', 'Platelets',
        'Prothrombin', 'Stage', 'edema1', 'edema2']



data_T = data[data["Drug"]=="D-penicillamine"].dropna()
data_C = data[data["Drug"]=="Placebo"].dropna()

#%%
X_t = data_T.loc[:,feature_list]
y_t = data_T.loc[:,['N_Days']]

X_t_scaler = StandardScaler(with_std=False)  
X_t_sc = X_t_scaler.fit_transform(X_t.iloc[:,4:])
y_t_scaler = StandardScaler(with_std=False)
y_t_sc = y_t_scaler.fit_transform(y_t)
#%%
def Lasso_alpha_bisection(X, y, initial_alpha = 1e10):
    alpha0 = initial_alpha
    res_alpha_list = [alpha0]
    _, coef_path, _ = lasso_path(X, y, alphas=[alpha0])
    coef_path = coef_path[0,:,:].T
    nonzero_index = np.where(coef_path[0]!=0)[0]
    current_val_num = len(nonzero_index)
    if current_val_num == 0 or 1:
        print("start bisection")
        alpha = alpha0/2
    else:
        print("initial value error")
        return -1
    while alpha>1e-10 and len(res_alpha_list) <= X.shape[1]:
        _, coef_path, _ = lasso_path(X, y, alphas=[alpha])
        coef_path = coef_path[0,:,:].T
        nonzero_index = np.where(coef_path[0]!=0)[0]
        if len(nonzero_index) <= current_val_num:
            alpha0 = alpha
            alpha = alpha/2
            continue
        elif len(nonzero_index) == current_val_num+1:
            res_alpha_list.append(alpha)
            current_val_num +=1
            alpha0 = alpha
            alpha = alpha/2
            continue
        elif len(nonzero_index) > current_val_num+1:
            alpha = alpha + (alpha0-alpha)/2
            continue
    return res_alpha_list

use_alpha = Lasso_alpha_bisection(X_t_sc,y_t_sc)
alpha, coef_path, _ = lasso_path(X_t_sc,y_t_sc,alphas=use_alpha)
coef_path = coef_path[0,:,:].T
current_val_num = 1
model_index_list_t = []
for i in range(len(coef_path)):
    nonzero_index = np.where(coef_path[i]!=0)[0]
    # print(alpha[i],len(nonzero_index),nonzero_index)
    if len(nonzero_index) == current_val_num :
        current_val_num+=1
        model_index_list_t.append(nonzero_index)

#%% model averaging for t
sample_size_t = y_t_sc.shape[0]
model_list_t = []
model_coef_list_t = []

for idx in model_index_list_t:
    lr_model = LinearRegression()
    lr_model.fit(X_t_sc[:,idx],y_t_sc)
    model_list_t.append(lr_model)
    model_coef_list_t.append(lr_model.coef_)

n = len(model_list_t)
model_pred_t = []
for i in range(n):
    model_pred_t.append((X_t_sc[:,model_index_list_t[i]] @ model_coef_list_t[i].T).reshape(-1) )
model_pred_t = np.array(model_pred_t)

sigma2 = mean_squared_error(y_t_sc.reshape(-1),model_pred_t[-1,:]) * sample_size_t / (sample_size_t-len(model_index_list_t[-1]))  
k = np.array([len(model_index_list_t[i]) for i in range(n)])
W = cp.Variable(n,pos=True)
def objective(W):
    return cp.sum_squares(W @ model_pred_t - y_t_sc.reshape(-1)) + 2*sigma2 * k @ W
constraint = [cp.sum(W)==1]

problem = cp.Problem(cp.Minimize(objective(W)),constraint)
problem.solve()
w_t = W.value

#%%
ma_pred_t_set = []
for i in range(len(model_index_list_t)):
    ma_pred_t_set.append(w_t[i] * model_pred_t[i,:])
ma_pred_t = np.sum(np.array(ma_pred_t_set),axis=0)
t_freedom = (np.where(w_t>1e-8)[0][-1] + 1)
sigma_ma_t = mean_squared_error(y_t_sc.reshape(-1),ma_pred_t) * sample_size_t /(sample_size_t - t_freedom )

#%% model averaging for c --------
X_c = data_C.loc[:,feature_list]
y_c = data_C.loc[:,['N_Days']]

X_c_scaler = StandardScaler(with_std=False)  
X_c_sc = X_c_scaler.fit_transform(X_c.iloc[:,4:])
y_c_scaler = StandardScaler(with_std=False)
y_c_sc = y_c_scaler.fit_transform(y_c)
model = LassoCV()
model.fit(X_c_sc,y_c_sc)

use_alpha = Lasso_alpha_bisection(X_c_sc,y_c_sc)
alpha, coef_path, _ = lasso_path(X_c_sc,y_c_sc,alphas=use_alpha)
coef_path = coef_path[0,:,:].T
current_val_num = 1
model_index_list_c = []
for i in range(len(coef_path)):
    nonzero_index = np.where(coef_path[i]!=0)[0]
    # print(alpha[i],len(nonzero_index),nonzero_index)
    if len(nonzero_index) == current_val_num :
        current_val_num+=1
        model_index_list_c.append(nonzero_index)

sample_size_c = y_c_sc.shape[0]
model_list_c = []
model_coef_list_c = []

for idx in model_index_list_c:
    lr_model = LinearRegression()
    lr_model.fit(X_c_sc[:,idx],y_c_sc)
    model_list_c.append(lr_model)
    model_coef_list_c.append(lr_model.coef_)

n = len(model_list_c)
model_pred_c = []
for i in range(n):
    model_pred_c.append((X_c_sc[:,model_index_list_c[i]] @ model_coef_list_c[i].T).reshape(-1) )
model_pred_c = np.array(model_pred_c)
sigma2 = mean_squared_error(y_c_sc.reshape(-1),model_pred_c[-1,:]) * sample_size_c / (sample_size_c-len(model_index_list_c[-1]))  
k = np.array([len(model_index_list_c[i]) for i in range(n)])
W = cp.Variable(n,pos=True)
def objective(W):
    return cp.sum_squares(W @ model_pred_c - y_c_sc.reshape(-1)) + 2*sigma2 * k @ W
constraint = [cp.sum(W)==1]

problem = cp.Problem(cp.Minimize(objective(W)),constraint)
problem.solve()
w_c = W.value

# %% calculate tau

data_all = pd.concat((data_C,data_T))
X_all = data_all.loc[:,feature_list]
x_mu = np.mean(X_all.iloc[:,4:],axis=0).values.reshape(1,-1)
x_t_mu = np.mean(X_t.iloc[:,4:],axis=0).values.reshape(1,-1)
x_c_mu = np.mean(X_c.iloc[:,4:],axis=0).values.reshape(1,-1)
y_t_mu = np.mean(y_t)
y_c_mu = np.mean(y_c)


ate = y_t_mu -y_c_mu
model_pred_mu_t = []
model_pred_mu_c = []
for i in range(n):
    model_pred_mu_t.append(( (x_t_mu - x_mu)[:,model_index_list_t[i]] @ model_coef_list_t[i].T  ).reshape(-1) )
    model_pred_mu_c.append(( (x_c_mu - x_mu)[:,model_index_list_c[i]] @ model_coef_list_c[i].T  ).reshape(-1) )
for i in range(len(model_index_list_t)):
    y_t_mu -= w_t[i] * model_pred_mu_t[i]
for i in range(len(model_index_list_c)):
    y_c_mu -= w_c[i] * model_pred_mu_c[i]
mu_diff = y_t_mu-y_c_mu
#%%
# for i in range(len(model_list_c)):
#     print(r2_score(y_c_sc,model_pred_c[i,:].reshape(-1,1)))
ma_pred_c_set = []
for i in range(len(model_index_list_c)):
    ma_pred_c_set.append(w_c[i] * model_pred_c[i,:])
ma_pred_c = np.sum(np.array(ma_pred_c_set),axis=0)
c_freedom = (np.where(w_c>1e-8)[0][-1] + 1)
sigma_ma_c = mean_squared_error(y_c_sc.reshape(-1),ma_pred_c) * sample_size_c /(sample_size_c - c_freedom )

# %%
sigma_ma = (sigma_ma_c    / sample_size_c + sigma_ma_t    / sample_size_t) ** 0.5
sigma_y = (np.var(data_T["N_Days"])/ sample_size_t  + np.var(data_C["N_Days"])/ sample_size_c ) ** 0.5
#%% Add lasso method

lasso_t = LassoCV()
lasso_t.fit(X_t_sc,y_t_sc)
lasso_c = LassoCV()
lasso_c.fit(X_c_sc,y_c_sc)
sigma_lasso_t = mean_squared_error(y_t_sc.reshape(-1) , X_t_sc @ lasso_t.coef_) * sample_size_t / (sample_size_t-np.sum(lasso_t.coef_ != 0))
sigma_lasso_c = mean_squared_error(y_c_sc.reshape(-1) , X_c_sc @ lasso_c.coef_) * sample_size_c / (sample_size_c-np.sum(lasso_c.coef_ != 0))

ate_lasso = (y_t_mu - (x_t_mu - x_mu) @ lasso_t.coef_ ) - (y_c_mu - (x_c_mu - x_mu) @ lasso_c.coef_)
sigma_lasso = (sigma_lasso_c  / sample_size_c + sigma_lasso_t  / sample_size_t) ** 0.5


# %%
plt.figure(figsize=(5,5))

my_xticks = ['Unadjusted',"Lasso",'Model Averaging']
plt.xticks([0,1,2], my_xticks)
plt.xlim(-0.5,2.5)
plt.vlines(2,ate-sigma_ma*1.96,ate+sigma_ma*1.96)
plt.vlines(1,ate_lasso-sigma_lasso*1.96,ate_lasso+sigma_lasso*1.96)
plt.vlines(0,mu_diff-sigma_y*1.96,mu_diff+sigma_y*1.96)
plt.hlines(0,-0.5,2.5,color = "black",linestyle ="dashed")
plt.plot(0,mu_diff,marker="o",color="red")
plt.plot(1,ate_lasso,marker="o",color="red")
plt.plot(2,ate,marker="o",color="red")
plt.hlines(mu_diff-sigma_y*1.96,-0.05,0.05)
plt.hlines(mu_diff+sigma_y*1.96,-0.05,0.05)
plt.hlines(ate_lasso-sigma_lasso*1.96,1-0.05,1+0.05)
plt.hlines(ate_lasso+sigma_lasso*1.96,1-0.05,1+0.05)
plt.hlines(ate-sigma_ma*1.96,2-0.05,2+0.05)
plt.hlines(ate+sigma_ma*1.96,2-0.05,2+0.05)
plt.ylabel("ATE estimates")
# plt.yscale('log')
# plt.xlabel("method")
plt.show()
# %%

print(r2_score(y_t_sc,w_c@ model_pred_t)) 
print(r2_score(y_c_sc,w_c@ model_pred_c))
# %%
