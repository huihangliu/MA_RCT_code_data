# MA_RCT_code_data

This folder is a cleaned reproducibility subset for the RCT model-averaging paper.

I built it by tracing the current manuscript in `../MA_RCT_R1` back to the legacy code under `/Users/huihang/Documents/Repository/RCT-MA/rct-ma-code`, then copying only the files that appear to support article figures or simulation tables.

## What Is Included

- `demo/article_demo_figure.R`
  Generates `demo_1.png` and `demo_2.png`, which match the conceptual two-panel illustration used in the manuscript.
- `design1/algorithms_0.1.r`
  Shared R helpers for the Design 1 simulations.
- `design1/design1_simulation.R`
  Legacy single-case Design 1 simulation runner for the numerical table. I removed the hard-coded `setwd()`, added a summary CSV export, and made the main settings runnable from command-line arguments.
- `design1/design1_results_reference.md`
  Historical notes/results that match the Bias/Var/MSE values appearing in the current paper table.
- `design2/simulation_design2.R`
  Main Design 2 simulation script. It produces the risk curves (`ate_*.pdf`) plus the CI summary CSV files used for the coverage/interval-length panels.
- `design2/plot_design2_from_csv.R`
  Cleaned plotting script for turning the article-used `CI_coverage_*.csv` and `CI_len_*.csv` files into the corresponding PDFs.
- `realdata/model_averaging.py`
  Self-contained real-data analysis script for the cirrhosis example. I changed the absolute input path to a local relative path and made it save `real_res.pdf`.
- `realdata/cirrhosis.csv`
  Input data for the real-data figure.

## Usage

For the demo figures, run:
```bash
cd demo
Rscript article_demo_figure.R
```

For Design 1, the current runner is `simulation/sim1.R`. 
The paper setting is `n=250`, balanced treatment/control groups (`125/125`), iid Gaussian covariates (`rho=0`), and `1e5` treatment re-randomization replications. Run a quick smoke
test first:

```bash
scripts/run_design1_local_debug.sh
```

For a local single-case debug run with paper-like dimensions:

```bash
N=250 P=50 S=10 RUN_ONCE=TRUE scripts/run_design1_local_debug.sh
```

For a paper-scale rerun on Precision1, launch the detached remote job:

```bash
scripts/run_design1_remote_paper.sh
```

The remote runner computes the parallel worker count as the available core count
minus 10, and `simulation/sim1.R` restricts BLAS to one thread per worker via
`RhpcBLASctl::blas_set_num_threads(1)`.

Each `sim1.R` run writes to a timestamped output directory:

```text
output/YYYYmmdd_HHMMSS/
  code/
  results/
  metadata/
```

For a multi-command paper-scale run, pass the same `run_dir=output/YYYYmmdd_HHMMSS`
to each command so all table rows land in one experiment directory. After the
remote run finishes, sync the output back to local:

```bash
scripts/sync_output_from_precision1.sh
```

Generate paper-ready Design 1 table files from a run directory:

```bash
Rscript simulation/draw_design1_results.R run_dir=output/YYYYmmdd_HHMMSS
```

This creates `draw/design1_table.tex`, `draw/design1_table_body.tex`,
`draw/design1_table.md`, and `draw/design1_table_values.csv` inside the run
directory.

```bash
cd design2
Rscript simulation_design2.R
Rscript plot_design2_from_csv.R
```

```bash
cd realdata
python3 model_averaging.py
```


## Environment

Create the local conda environment directly at this repository path:

```bash
conda env create -p ./.venv -f environment.yml
conda activate ./.venv
```

Do not copy an existing conda environment into `.venv`; conda environments contain absolute prefix paths, and copied environments can keep pointing back to the original location.

## Dependencies

R side:

- `snowfall`
- `glmnet`
- `MASS`
- `purrr`
- `tidyr`
- `quadprog`
- `RhpcBLASctl`
- `ncvreg`
- `knitr`
- `rmarkdown`
- `bookdown`

Python side:

- `numpy`
- `pandas`
- `matplotlib`
- `scikit-learn`
- `cvxpy`
