rm(list = ls())

get_script_dir <- function() {
  args <- commandArgs(trailingOnly = FALSE)
  script_arg <- grep("^--file=", args, value = TRUE)
  if (length(script_arg) == 0) {
    return(normalizePath("."))
  }
  dirname(normalizePath(sub("^--file=", "", script_arg[1])))
}

script_dir <- get_script_dir()

set.seed(2026)
n <- 50
n_T <- 20
n_C <- n - n_T
X <- rnorm(n)
e <- 0.5 * rnorm(n)

idx_treated <- sample(n, n_T, replace = FALSE)
idx_control <- setdiff(seq_len(n), idx_treated)

T <- rep(0, n)
T[idx_treated] <- 1
y <- T + X + e

ATE_unadj <- mean(y[idx_treated]) - mean(y[idx_control])

X_tmp <- X[idx_treated] - mean(X[idx_treated])
y_tmp <- y[idx_treated] - mean(y[idx_treated])
beta_hat_treated <- solve(t(X_tmp) %*% X_tmp) %*% t(X_tmp) %*% y_tmp

X_tmp <- X[idx_control] - mean(X[idx_control])
y_tmp <- y[idx_control] - mean(y[idx_control])
beta_hat_control <- solve(t(X_tmp) %*% X_tmp) %*% t(X_tmp) %*% y_tmp

ATE_adj <- (mean(y[idx_treated]) - (mean(X[idx_treated]) - mean(X)) * beta_hat_treated) -
  (mean(y[idx_control]) - (mean(X[idx_control]) - mean(X)) * beta_hat_control)

# print
print(ATE_unadj)
print(ATE_adj)


# -----------------------------
# nicer plot settings
# -----------------------------
nice_limits <- function(z, pad = 0.08) {
  r <- range(z, finite = TRUE)
  d <- diff(r)
  if (d == 0 || !is.finite(d)) d <- 1
  c(r[1] - pad * d, r[2] + pad * d)
}

col_control <- "#D55E00"   # orange
col_treated <- "#0072B2"   # blue
col_mean    <- "#111111"
col_ref     <- "#6A3D9A"

pch_control <- 16          # filled circle
pch_treated <- 17          # filled triangle

# -----------------------------
# plot 1: unadjusted comparison
# -----------------------------
pdf(
  file.path(script_dir, "demo_1.pdf"),
  width = 7,
  height = 5,
  version = "1.4",
  useDingbats = FALSE
)

par(
  mar = c(4.5, 4.8, 3.8, 1.5),
  las = 1,
  mgp = c(2.6, 0.8, 0),
  cex = 1.05
)

set.seed(2026)
x_control <- jitter(T[idx_control], amount = 0.06)
x_treated <- jitter(T[idx_treated], amount = 0.06)

plot(
  0, 0,
  type = "n",
  xlim = c(-0.45, 1.45),
  ylim = nice_limits(y),
  xaxt = "n",
  xlab = "",
  ylab = "y",
  main = bquote(widehat(ATE)[unadj] == .(round(ATE_unadj, 2))),
  bty = "l"
)

axis(1, at = c(0, 1), labels = c("Control", "Treatment"))
grid(nx = NA, ny = NULL, col = adjustcolor("gray80", alpha.f = 0.6), lty = 1)

abline(v = c(0, 1), col = adjustcolor("gray70", alpha.f = 0.35), lty = 3)

points(
  x_control,
  y[idx_control],
  pch = pch_control,
  col = adjustcolor(col_control, alpha.f = 0.75),
  cex = 1.15
)

points(
  x_treated,
  y[idx_treated],
  pch = pch_treated,
  col = adjustcolor(col_treated, alpha.f = 0.75),
  cex = 1.25
)

mean_control <- mean(y[idx_control])
mean_treated <- mean(y[idx_treated])

segments(
  -0.18, mean_control, 0.18, mean_control,
  col = col_control,
  lwd = 4,
  lend = "round"
)

segments(
  0.82, mean_treated, 1.18, mean_treated,
  col = col_treated,
  lwd = 4,
  lend = "round"
)

points(0, mean_control, pch = 18, cex = 2.1, col = col_mean)
points(1, mean_treated, pch = 18, cex = 2.1, col = col_mean)

legend(
  "topleft",
  legend = c("Control", "Treatment", "Group mean"),
  col = c(col_control, col_treated, col_mean),
  pch = c(pch_control, pch_treated, 18),
  pt.cex = c(1.2, 1.25, 1.8),
  bty = "n"
)

mtext(
  "Horizontal bars mark group means",
  side = 3,
  adj = 1,
  cex = 0.8,
  col = "gray35"
)

dev.off()


# -----------------------------
# plot 2: covariate adjustment
# -----------------------------
fit_band <- function(idx, n_grid = 200) {
  x_center <- X[idx] - mean(X[idx])
  y_center <- y[idx] - mean(y[idx])
  
  df <- data.frame(
    x_center = x_center,
    y_center = y_center
  )
  
  mod <- lm(y_center ~ x_center, data = df)
  
  newx <- seq(min(x_center), max(x_center), length.out = n_grid)
  
  preds <- predict(
    mod,
    newdata = data.frame(x_center = newx),
    interval = "confidence"
  )
  
  data.frame(
    x   = newx + mean(X[idx]),
    fit = preds[, "fit"] + mean(y[idx]),
    lwr = preds[, "lwr"] + mean(y[idx]),
    upr = preds[, "upr"] + mean(y[idx])
  )
}

band_treated <- fit_band(idx_treated)
band_control <- fit_band(idx_control)

pdf(
  file.path(script_dir, "demo_2.pdf"),
  width = 7,
  height = 5,
  version = "1.4",
  useDingbats = FALSE
)

par(
  mar = c(4.8, 4.8, 3.8, 1.5),
  las = 1,
  mgp = c(2.6, 0.8, 0),
  cex = 1.05
)

plot(
  0, 0,
  type = "n",
  xlim = nice_limits(X),
  ylim = nice_limits(y),
  xlab = "x",
  ylab = "y",
  main = bquote("Covariate-adjusted fits, " ~ widehat(ATE)[adj] == .(round(as.numeric(ATE_adj), 2))),
  bty = "l"
)

grid(col = adjustcolor("gray80", alpha.f = 0.6), lty = 1)

# Confidence bands first, so they stay behind the points and lines
polygon(
  c(band_control$x, rev(band_control$x)),
  c(band_control$lwr, rev(band_control$upr)),
  col = adjustcolor(col_control, alpha.f = 0.18),
  border = NA
)

polygon(
  c(band_treated$x, rev(band_treated$x)),
  c(band_treated$lwr, rev(band_treated$upr)),
  col = adjustcolor(col_treated, alpha.f = 0.18),
  border = NA
)

# Regression lines
lines(
  band_control$x,
  band_control$fit,
  col = col_control,
  lwd = 2.8
)

lines(
  band_treated$x,
  band_treated$fit,
  col = col_treated,
  lwd = 2.8
)

# Raw points
points(
  X[idx_control],
  y[idx_control],
  pch = pch_control,
  col = adjustcolor(col_control, alpha.f = 0.75),
  cex = 1.15
)

points(
  X[idx_treated],
  y[idx_treated],
  pch = pch_treated,
  col = adjustcolor(col_treated, alpha.f = 0.75),
  cex = 1.25
)

# Group means
points(
  mean(X[idx_control]),
  mean(y[idx_control]),
  pch = 18,
  cex = 2.3,
  col = col_control
)

points(
  mean(X[idx_treated]),
  mean(y[idx_treated]),
  pch = 18,
  cex = 2.3,
  col = col_treated
)

# Overall mean of X
abline(
  v = mean(X),
  col = col_ref,
  lty = 2,
  lwd = 2
)

legend(
  "topleft",
  legend = c("Control", "Treatment"),
  col = c(col_control, col_treated),
  pch = c(pch_control, pch_treated),
  lty = c(1, 1),
  lwd = c(2.8, 2.8),
  pt.cex = c(1.15, 1.25),
  bty = "n"
)

legend(
  "bottomright",
  legend = c("Group mean", "Overall mean of x"),
  col = c(col_mean, col_ref),
  pch = c(18, NA),
  lty = c(NA, 2),
  lwd = c(NA, 2),
  pt.cex = c(1.8, NA),
  bty = "n"
)

dev.off()