# Verification Notes

Last checked: 2026-05-12

## Environment

- `R` and `Rscript` were installed locally on this Mac via `micromamba` and exposed on `PATH` through `~/.local/bin`.
- The verification environment includes the R packages used by this repo (`glmnet`, `snowfall`, `MASS`, `quadprog`, `purrr`, `tidyr`) and the Python packages needed by `realdata/model_averaging.py` (`numpy`, `pandas`, `matplotlib`, `scikit-learn`, `cvxpy`).

## Summary

- `design2/version-240806/*.csv` matches the manuscript asset folder `../MA_RCT_R1/img/version-240806/*.csv` exactly, file by file.
- `realdata/model_averaging.py` required fixes before it could run reproducibly; those fixes are now in this repo.
- `design1/design1_simulation.R` required a small reproducibility fix for parallel execution; that fix is now in this repo.
- The Design 1 paper table is still the least clean part of the project: the current script lineage, `design1_results_reference.md`, and `find_seed_2021.05.27.txt` do not resolve to a single clean one-command pipeline for every published row.

## Detailed Checks

### Demo figures

- The current working-tree version of `demo/article_demo_figure.R` is not the original article script: it has uncommitted local edits, including a different seed, a different treated-sample size, and plotting changes.
- As a result, the currently generated `demo_1.png` and `demo_2.png` in this working tree should not be treated as canonical article reproductions.

### Design 1 table

- `design1/design1_simulation.R` now accepts command-line overrides such as `seed=...`, `n=...`, `p=...`, `s=...`, `rho=...`, `rep_num=...`, and `cpus=...`.
- The script also now exports `script_dir` to the snowfall workers; without that, the parallel run fails because the workers cannot source `algorithms_0.1.r`.
- Historical note: the `ATE_true` value `-0.4460537` recorded in `design1_results_reference.md` for the `n=250, p=50, s=10, rho=0` case corresponds to `seed=1001`.
- Historical note: `design1/find_seed_2021.05.27.txt` separately points to `seed=43` as matching the published variance values for the same low-dimensional case.
- Because those two historical records conflict, the exact row-by-row paper pipeline is not uniquely recoverable from the cleaned code alone. This repo should therefore treat `design1_results_reference.md` and `find_seed_2021.05.27.txt` as provenance notes, not as proof of a single frozen simulation state.

Example verified command:

```bash
cd design1
Rscript design1_simulation.R seed=1001 n=250 p=50 s=10 rho=0 rep_num=1000 cpus=4
```

This command now completes successfully and writes a summary CSV, but its 1000-rep summary is not numerically identical to the published table row. That is consistent with the historical provenance gap described above.

### Design 2 figures

- The numeric inputs for the published coverage and interval-length panels are reproducible: the code-repo CSVs are exact copies of the manuscript asset CSVs.
- The cleaned plotting script `design2/plot_design2_from_csv.R` regenerates the panels from those CSVs, but the resulting PDFs are not pixel-identical to the PDFs under `../MA_RCT_R1/img/version-240806/`.
- This indicates that the numeric results are preserved, but the exact final plotting/export pipeline used for the manuscript PDFs is not fully preserved in the cleaned code repo.

### Real-data figure

- `realdata/model_averaging.py` originally failed because it mixed array-valued predictions with scalar subtraction and also reused the model-averaging-adjusted means inside the Lasso ATE calculation.
- The script now runs successfully and produces `real_res.pdf`.
- The figure size was adjusted to `figsize=(4.7, 3.4)`, which brings the generated PDF much closer to the manuscript asset size.
- The generated PDF is still not byte-identical to `../MA_RCT_R1/img/real_res.pdf`, but the gap is now largely a rendering/layout issue rather than a runtime failure.

Verified output from the current script:

```text
ate_unadj   = -0.1601232589
ate_lasso   =  0.0745959167
ate_ma      = -0.8623474131
sigma_y     =  1.1750230045
sigma_lasso =  0.9782631328
sigma_ma    =  0.7981232671
```
