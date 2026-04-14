# Output Map

This file maps the current manuscript outputs in `../MA_RCT_R1/RCT-MA_manu.tex` to the code or data files collected in this repo.

## Figures

- `fig:demo1` (`img/dag_demo.pdf`)
  Status: manuscript asset only.
  Code found: none in the legacy code folder.
  Confidence: low.

- `fig:demo` (`img/demo_1.png`, `img/demo_2.png`)
  Source: `demo/article_demo_figure.R`
  Legacy origin: `RCT-Demo__Code_1.r` / `RCT-Demo.Rmd`
  Confidence: high.

- `fig:sim2_res`
  Source: `design2/simulation_design2.R`
  Article-used files: mixed `ate_*.pdf` assets under the manuscript's `img/version-240806/`
  Confidence: high on script lineage, medium on exact seed-by-panel selection.

- `fig:sim2_res_CI_coverage`
  Source: `design2/plot_design2_from_csv.R`
  Input CSVs: `design2/version-240806/CI_coverage_*.csv`
  Confidence: high.

- `fig:sim2_res_CI_len`
  Source: `design2/plot_design2_from_csv.R`
  Input CSVs: `design2/version-240806/CI_len_*.csv`
  Confidence: high.

- `fig:real_res`
  Source: `realdata/model_averaging.py`
  Input data: `realdata/cirrhosis.csv`
  Confidence: high.

## Tables

- `tab:sim1_method`
  Source: manual LaTeX table in the manuscript.
  Code found: none needed.
  Confidence: high.

- `tab:sim1_res`
  Source: `design1/design1_simulation.R` plus `design1/algorithms_0.1.r`
  Historical paper-matching values: `design1/design1_results_reference.md`
  Notes: this is the messiest provenance path; the final paper table appears to have been assembled from simulation runs plus historical notes rather than from one preserved end-to-end script.
  Confidence: medium.

- `tab:real_variable`
  Source: manual LaTeX table in the manuscript.
  Code found: none needed.
  Confidence: high.

- `tab:real_lifetime`
  Source: manual LaTeX table in the manuscript.
  Code found: none needed.
  Confidence: high.
