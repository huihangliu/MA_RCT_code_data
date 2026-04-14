# MA_RCT_code_data

This folder is a cleaned reproducibility subset for the RCT model-averaging paper.

I built it by tracing the current manuscript in `../MA_RCT_R1` back to the legacy code under `/Users/huihang/Documents/Repository/RCT-MA/rct-ma-code`, then copying only the files that appear to support article figures or simulation tables.

## What Is Included

- `demo/article_demo_figure.R`
  Generates `demo_1.png` and `demo_2.png`, which match the conceptual two-panel illustration used in the manuscript.
- `design1/algorithms_0.1.r`
  Shared R helpers for the Design 1 simulations.
- `design1/design1_simulation.R`
  Legacy single-case Design 1 simulation runner for the numerical table. I removed the hard-coded `setwd()` and added a summary CSV export.
- `design1/design1_results_reference.md`
  Historical notes/results that match the Bias/Var/MSE values appearing in the current paper table.
- `design1/find_seed_2021.05.27.txt`
  Historical seed search notes that explain at least part of the Design 1 table provenance.
- `design2/simulation_design2.R`
  Main Design 2 simulation script. It produces the risk curves (`ate_*.pdf`) plus the CI summary CSV files used for the coverage/interval-length panels.
- `design2/plot_design2_from_csv.R`
  Cleaned plotting script for turning the article-used `CI_coverage_*.csv` and `CI_len_*.csv` files into the corresponding PDFs.
- `design2/version-240806/*.csv`
  The exact CSV summaries copied from the manuscript asset folder for the published Design 2 coverage and interval-length figures.
- `realdata/model_averaging.py`
  Self-contained real-data analysis script for the cirrhosis example. I changed the absolute input path to a local relative path and made it save `real_res.pdf`.
- `realdata/cirrhosis.csv`
  Input data for the real-data figure.
- `OUTPUT_MAP.md`
  A figure/table-to-code map with notes on confidence and missing provenance.

## What I Found

- The manuscript's Design 2 figure panels are supported by the `simulation_v2.r` lineage, not by the older `simulation_v1.r`.
- The real-data figure comes from `realdata/model_averaging.py`, but the original file had an absolute local path that would have broken outside the old machine.
- The Design 1 table values are partly recoverable from historical notes (`design1_results_reference.md`) and partly from the legacy simulation code, but the exact final row-by-row paper pipeline is not cleanly preserved in a single script.
- `dag_demo.pdf` is used by the manuscript, but I did not find a script in the legacy code directory that generates it. It looks like a manually prepared manuscript asset rather than a code-generated figure.
- The manuscript's simulation-method and real-data-description tables are typed directly in LaTeX, so there is no separate code file to include for those tables.

## Suggested Usage

From this folder:

```bash
cd design2
Rscript simulation_design2.R
Rscript plot_design2_from_csv.R
```

```bash
cd demo
Rscript article_demo_figure.R
```

```bash
cd realdata
python3 model_averaging.py
```

For Design 1, edit the settings near the top of `design1/design1_simulation.R` before running:

```bash
cd design1
Rscript design1_simulation.R
```

## Dependencies

R side:

- `snowfall`
- `glmnet`
- `MASS`
- `purrr`
- `tidyr`
- `quadprog`
- `RhpcBLASctl`

Python side:

- `numpy`
- `pandas`
- `matplotlib`
- `scikit-learn`
- `cvxpy`

## Limits

- I could not directly `git clone` `https://github.com/huihangliu/MA_ATE` from this shell because the repo is private and the terminal session had no GitHub credentials. To keep moving, I used the local manuscript source at `/Users/huihang/Documents/Repository/RCT-MA/rct-ma-latex` and mirrored it into `../MA_RCT_R1`.
- The code provenance for Design 1 is the least clean part of the project. I included the best-matching script plus the historical result notes instead of pretending the lineage is cleaner than it is.
