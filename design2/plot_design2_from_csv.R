get_script_dir <- function() {
  args <- commandArgs(trailingOnly = FALSE)
  script_arg <- grep("^--file=", args, value = TRUE)
  if (length(script_arg) == 0) {
    return(normalizePath("."))
  }
  dirname(normalizePath(sub("^--file=", "", script_arg[1])))
}

script_dir <- get_script_dir()
input_dir <- file.path(script_dir, "version-240806")

n_list <- c(150, 200, 250, 300)
alpha <- 1.5
rho_list <- c(0.3, 0.6)
num_rep <- 2000
seed <- 43

plot_metric <- function(prefix, y_label, log_scale = FALSE) {
  for (n in n_list) {
    for (rho in rho_list) {
      csv_path <- file.path(input_dir, paste0(prefix, "_", n, "_", alpha, "_", rho, "_", num_rep, "_", seed, ".csv"))
      pdf_path <- file.path(input_dir, paste0(prefix, "_", n, "_", alpha, "_", rho, "_", num_rep, "_", seed, ".pdf"))

      data <- read.csv(csv_path)
      ylim_min <- min(data[[2]], data[[3]], data[[4]], data[[5]], data[[6]], data[[7]], data[[8]], data[[9]])
      ylim_max <- max(data[[2]], data[[3]], data[[4]], data[[5]], data[[6]], data[[7]], data[[8]], data[[9]])

      if (prefix == "CI_coverage") {
        ylim_min <- ylim_min - 0.02
        ylim_max <- ylim_max + 0.02
      } else {
        ylim_min <- max(ylim_min - 0.005, .Machine$double.eps)
        ylim_max <- ylim_max + 0.5
      }

      pdf(pdf_path, width = 6, height = 6)
      plot(
        data[[1]],
        data[[2]],
        type = "o",
        lty = 1,
        pch = 0,
        lwd = 2,
        ylim = c(ylim_min, ylim_max),
        col = "black",
        main = bquote(n == .(n) ~ ", " ~ alpha == .(alpha) ~ " and " ~ rho == .(rho)),
        xlab = expression(R^2),
        ylab = y_label,
        log = if (log_scale) "y" else ""
      )
      lines(data[[1]], data[[3]], type = "b", lty = 2, pch = 2, lwd = 2, col = "blue")
      lines(data[[1]], data[[4]], type = "b", lty = 3, pch = 4, lwd = 2, col = "green")
      lines(data[[1]], data[[5]], type = "b", lty = 4, pch = 6, lwd = 2, col = "orange")
      lines(data[[1]], data[[6]], type = "b", lty = 5, pch = 8, lwd = 2, col = "red")
      lines(data[[1]], data[[7]], type = "b", lty = 6, pch = 10, lwd = 2, col = "darkred")
      lines(data[[1]], data[[8]], type = "b", lty = 7, pch = 12, lwd = 2, col = "purple")
      lines(data[[1]], data[[9]], type = "b", lty = 8, pch = 14, lwd = 2, col = "brown")
      legend(
        "topleft",
        legend = c("MMA", "AIC", "BIC", "S-AIC", "S-BIC", "FULL", "LASSO", "unadj"),
        col = c("black", "blue", "green", "orange", "red", "darkred", "purple", "brown"),
        lty = 1:8,
        lwd = 2,
        pch = c(0, 2, 4, 6, 8, 10, 12, 14)
      )
      if (prefix == "CI_coverage") {
        abline(h = 0.95, lty = 2, col = "gray40")
      }
      dev.off()
    }
  }
}

plot_metric("CI_coverage", "Coverage Rate", log_scale = FALSE)
plot_metric("CI_len", "Interval Length", log_scale = TRUE)
