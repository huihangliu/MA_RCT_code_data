# MA_RCT_code_data

This folder is the code repo for RCT model-averaging.

## What Is Included

- `environment.yml`
  Conda environment specification for the R and Python dependencies used by the simulations, plotting scripts, and real-data analysis.
- `simulation/algorithms_0.1.r`
  Shared R helper functions used by the Design 1 simulation runner.
- `simulation/article_demo_figure.R`
  Generates the two simple conceptual demo figures, `simulation/demo_1.png` and `simulation/demo_2.png`.
- `simulation/demo_plot.R`
  Generates polished PDF versions of the same demo illustration, `simulation/demo_1.pdf` and `simulation/demo_2.pdf`.
- `simulation/sim1.R`
  Design 1 simulation runner for the paper table. It accepts `key=value` command-line arguments, writes timestamped output directories, and snapshots the code used for each run.
- `simulation/draw_design1_results.R`
  Converts Design 1 summary CSV files into paper-ready table outputs.
- `simulation/sim2.R`
  Design 2 simulation runner for risk, CI length, CI coverage, and ATE histogram outputs. It supports a quick `debug` mode and a full `parallel` mode.
- `simulation/draw_design2_results.py`
  Plots Design 2 `risk_*.csv`, `CI_len_*.csv`, and `CI_coverage_*.csv` summaries into by-seed and seed-average PDF figures.
- `realdata/model_averaging.py`
  Self-contained real-data analysis script for the cirrhosis example. It reads the local CSV and saves `realdata/real_res.pdf`.
- `realdata/cirrhosis.csv`
  Input data for the real-data figure.

## Usage

Create the local environment from the repository root:

```bash
conda env create -p ./.venv -f environment.yml
conda activate ./.venv
```

For the demo figures, run:

```bash
Rscript simulation/article_demo_figure.R
Rscript simulation/demo_plot.R
```

For a quick Design 1 smoke test:

```bash
Rscript simulation/sim1.R seed=43 n=50 p=10 s=2 rho=0 rep_num=1 cpus=1 run_once=TRUE output_root=output
```

For a local Design 1 run with paper-like dimensions:

```bash
Rscript simulation/sim1.R seed=43 n=250 p=50 s=10 rho=0 rep_num=100 cpus=4 run_once=FALSE output_root=output
```

Each `sim1.R` run writes a timestamped directory under `output/` unless `run_dir=...` is supplied:

```text
output/YYYYmmdd_HHMMSS/
  code/
  results/
  metadata/
```

Generate Design 1 table files from a run directory:

```bash
Rscript simulation/draw_design1_results.R run_dir=output/YYYYmmdd_HHMMSS
```

This creates `draw/design1_table.tex`, `draw/design1_table_body.tex`, `draw/design1_table.md`, and `draw/design1_table_values.csv` inside the run directory.

For a quick Design 2 debug run:

```bash
Rscript simulation/sim2.R debug 200 1 2001 0.5 2024
```

For the full Design 2 grid:

```bash
Rscript simulation/sim2.R parallel 2000 2024 8
```

The last argument is the number of worker cores. If omitted, `sim2.R` uses the detected core count minus 10.

Design 2 outputs are written under `output/` as `risk_*.csv`, `risk_*.pdf`, `CI_len_*.csv`, `CI_coverage_*.csv`, `ate_hat_mma_*.csv`, and `ate_hist_*.pdf`.

Plot Design 2 summaries from a directory containing the CSV files:

```bash
python simulation/draw_design2_results.py output output/design2_plots
```

When comparing multiple seeds, place the corresponding summary CSV files in the same input directory before running the plotting script. It will create both `by_seed/` and `seed_average/` PDF figures.

Run the real-data analysis:

```bash
python realdata/model_averaging.py
```

Generated outputs such as PDFs, CSVs, logs, and `output/` directories are ignored by git.
