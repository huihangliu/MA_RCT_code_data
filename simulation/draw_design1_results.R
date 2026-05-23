# Generate paper-ready Design 1 tables from archived simulation results.

parse_cli_args <- function(args) {
  if (length(args) == 0) {
    return(list())
  }
  parsed <- list()
  for (arg in args) {
    if (!grepl("=", arg, fixed = TRUE)) {
      next
    }
    key <- sub("=.*$", "", arg)
    value <- sub("^[^=]*=", "", arg)
    parsed[[key]] <- type.convert(value, as.is = TRUE)
  }
  parsed
}

get_script_dir <- function() {
  args <- commandArgs(trailingOnly = FALSE)
  script_arg <- grep("^--file=", args, value = TRUE)
  if (length(script_arg) == 0) {
    return(normalizePath("."))
  }
  dirname(normalizePath(sub("^--file=", "", script_arg[1])))
}

is_absolute_path <- function(path) {
  grepl("^(/|[A-Za-z]:[\\/])", path)
}

resolve_path <- function(path, base_dir) {
  if (is_absolute_path(path)) {
    normalizePath(path, mustWork = FALSE)
  } else {
    normalizePath(file.path(base_dir, path), mustWork = FALSE)
  }
}

latest_output_dir <- function(output_root) {
  dirs <- list.dirs(output_root, recursive = FALSE, full.names = TRUE)
  dirs <- dirs[dir.exists(file.path(dirs, "results"))]
  if (length(dirs) == 0) {
    stop("No output run directories with a results/ subdirectory found.")
  }
  dirs[order(basename(dirs), decreasing = TRUE)][1]
}

parse_result_filename <- function(path) {
  name <- basename(path)
  pattern <- "^design1_summary_n([0-9]+)_p([0-9]+)_s([0-9]+)_rho([^_]+)_seed([0-9]+)_rep([0-9]+)\\.csv$"
  match <- regexec(pattern, name)
  parts <- regmatches(name, match)[[1]]
  if (length(parts) != 7) {
    stop("Unexpected result filename: ", name)
  }
  list(
    n = as.integer(parts[2]),
    p = as.integer(parts[3]),
    s = as.integer(parts[4]),
    rho = parts[5],
    seed = as.integer(parts[6]),
    rep = as.integer(parts[7])
  )
}

format_number <- function(value, digits = 4) {
  formatC(as.numeric(value), format = "f", digits = digits)
}

format_cell <- function(value, metric, ols_infeasible = FALSE) {
  if (ols_infeasible) {
    return("--")
  }
  if (metric == "Coverage_pct") {
    return(paste0(format_number(value, 2), "\\%"))
  }
  format_number(value, 4)
}

latex_line <- function(parts) {
  paste0(paste(parts, collapse = " & "), " \\\\")
}

write_markdown_table <- function(rows, out_file) {
  header <- c("Setting", "Metric", "un-adj", "OLS", "Lasso", "Lasso-OLS", "MA")
  lines <- c(
    paste(header, collapse = " | "),
    paste(rep("---", length(header)), collapse = " | ")
  )
  for (ii in seq_len(nrow(rows))) {
    lines <- c(lines, paste(as.character(rows[ii, header]), collapse = " | "))
  }
  writeLines(lines, out_file)
}

script_dir <- get_script_dir()
repo_dir <- dirname(script_dir)
args <- parse_cli_args(commandArgs(trailingOnly = TRUE))

output_root <- resolve_path(if (!is.null(args$output_root)) args$output_root else "output", repo_dir)
run_dir <- if (!is.null(args$run_dir)) {
  resolve_path(args$run_dir, repo_dir)
} else {
  latest_output_dir(output_root)
}
results_dir <- file.path(run_dir, "results")
draw_dir <- file.path(run_dir, "draw")
code_draw_dir <- file.path(run_dir, "code", "draw")

dir.create(draw_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(code_draw_dir, recursive = TRUE, showWarnings = FALSE)

script_path <- file.path(script_dir, "draw_design1_results.R")
if (file.exists(script_path)) {
  file.copy(script_path, code_draw_dir, overwrite = TRUE)
}

result_files <- list.files(
  results_dir,
  pattern = "^design1_summary_n[0-9]+_p[0-9]+_s[0-9]+_rho[^_]+_seed[0-9]+_rep[0-9]+\\.csv$",
  full.names = TRUE
)
if (length(result_files) == 0) {
  stop("No Design 1 summary CSV files found in ", results_dir)
}

methods <- c("unadj", "ols", "lasso", "lasso_ols", "ma")
method_labels <- c(
  unadj = "un-adj",
  ols = "OLS",
  lasso = "Lasso",
  lasso_ols = "Lasso-OLS",
  ma = "MA"
)
metrics <- c("Bias", "Var", "MSE", "Length", "Coverage_pct")
metric_labels <- c(
  Bias = "Bias",
  Var = "Var",
  MSE = "MSE",
  Length = "{Length}",
  Coverage_pct = "{Coverage}"
)

records <- list()
for (file in result_files) {
  cfg <- parse_result_filename(file)
  table <- read.csv(file, row.names = 1, check.names = FALSE)
  missing_metrics <- setdiff(metrics, rownames(table))
  missing_methods <- setdiff(methods, colnames(table))
  if (length(missing_metrics) > 0 || length(missing_methods) > 0) {
    stop("Missing expected rows or columns in ", file)
  }
  setting <- sprintf("$p=%d,s=%d$", cfg$p, cfg$s)
  plain_setting <- sprintf("p=%d,s=%d", cfg$p, cfg$s)
  ols_infeasible <- cfg$p >= cfg$n / 2
  for (metric in metrics) {
    record <- data.frame(
      n = cfg$n,
      p = cfg$p,
      s = cfg$s,
      rho = cfg$rho,
      seed = cfg$seed,
      rep = cfg$rep,
      Setting = plain_setting,
      Metric = metric,
      check.names = FALSE
    )
    for (method in methods) {
      record[[method]] <- as.numeric(table[metric, method])
    }
    record[["ols_infeasible"]] <- ols_infeasible
    record[["latex_setting"]] <- setting
    records[[length(records) + 1]] <- record
  }
}

combined <- do.call(rbind, records)
combined <- combined[order(combined$p, combined$s, match(combined$Metric, metrics)), ]
write.csv(
  combined[, c("n", "p", "s", "rho", "seed", "rep", "Metric", methods)],
  file.path(draw_dir, "design1_table_values.csv"),
  row.names = FALSE
)

paper_rows <- data.frame(
  Setting = character(0),
  Metric = character(0),
  `un-adj` = character(0),
  OLS = character(0),
  Lasso = character(0),
  `Lasso-OLS` = character(0),
  MA = character(0),
  check.names = FALSE
)

latex_body <- c()
settings <- unique(combined[, c("p", "s", "latex_setting")])
settings <- settings[order(settings$p, settings$s), ]
for (setting_idx in seq_len(nrow(settings))) {
  setting_rows <- combined[combined$p == settings$p[setting_idx] & combined$s == settings$s[setting_idx], ]
  for (metric_idx in seq_along(metrics)) {
    metric <- metrics[metric_idx]
    row <- setting_rows[setting_rows$Metric == metric, ]
    cells <- vapply(
      methods,
      function(method) format_cell(row[[method]], metric, method == "ols" && isTRUE(row$ols_infeasible)),
      character(1)
    )

    setting_cell <- if (metric_idx == 1) {
      sprintf("\\multirow{5}{*}{%s}", row$latex_setting)
    } else {
      ""
    }
    latex_body <- c(latex_body, latex_line(c(setting_cell, metric_labels[[metric]], cells)))

    paper_rows[nrow(paper_rows) + 1, ] <- c(
      if (metric_idx == 1) row$Setting else "",
      gsub("[{}]", "", metric_labels[[metric]]),
      cells[["unadj"]],
      cells[["ols"]],
      cells[["lasso"]],
      cells[["lasso_ols"]],
      cells[["ma"]]
    )
  }
  if (setting_idx < nrow(settings)) {
    latex_body <- c(latex_body, "    \\hline")
  }
}

latex_body <- paste0("    ", latex_body)
writeLines(latex_body, file.path(draw_dir, "design1_table_body.tex"))

latex_table <- c(
  "\\begin{table}[htbp]",
  "  \\setlength{\\tabcolsep}{4.5mm}",
  "  \\caption{Results of numerical simulation in \\emph{Design 1}}",
  "  \\label{tab:sim1_res}",
  "  \\centering",
  "  \\footnotesize",
  "  \\begin{tabular}{ccccccc}",
  "    \\toprule",
  "    & & un-adj & OLS & Lasso & Lasso-OLS & MA \\\\",
  "    \\midrule",
  latex_body,
  "    \\bottomrule",
  "  \\end{tabular}",
  "  \\begin{tablenotes}",
  "    \\footnotesize",
  "    \\item Note: OLS only applies for cases of $n>p$, thus infeasible results are marked by ``--''. Bias represents the mean value of estimation error, Var represents the variance of estimation error, and MSE represents the mean square error of estimates. We keep $4$ digits after the decimal point.",
  "    The tuning parameter in PMA is set as $\\lambda=2$ in the low-dimensional cases and $4\\log(n)$ in the high-dimensional cases.",
  "  \\end{tablenotes}",
  "\\end{table}"
)
writeLines(latex_table, file.path(draw_dir, "design1_table.tex"))

write_markdown_table(paper_rows, file.path(draw_dir, "design1_table.md"))

manifest <- c(
  paste("run_dir:", run_dir),
  paste("results_dir:", results_dir),
  paste("draw_dir:", draw_dir),
  paste("generated_at:", format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z")),
  "outputs:",
  "  design1_table.tex",
  "  design1_table_body.tex",
  "  design1_table.md",
  "  design1_table_values.csv",
  "code:",
  paste(" ", file.path(code_draw_dir, "draw_design1_results.R"))
)
writeLines(manifest, file.path(draw_dir, "MANIFEST.txt"))

cat("Draw output:", draw_dir, "\n")
