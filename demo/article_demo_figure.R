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

# plot
png(file.path(script_dir, "demo_1.png"), width = 1200, height = 900, res = 150)

plot(
  T[idx_treated],
  y[idx_treated],
  xlim = c(-0.5, 1.5),
  ylim = c(min(y), max(y)),
  ylab = "y",
  xlab = "",  # remove x label
  col = rgb(0, 0, 255, 150, maxColorValue = 255),
  main = paste("ATE hat =", round(ATE_unadj, 2)),
  xaxt = "n"  # hide x axis
)
axis(1, at = c(0, 1), labels = c("Control", "Treatment"))
points(T[idx_control], y[idx_control], col = rgb(0, 0, 0, 150, maxColorValue = 255))
points(1, mean(y[idx_treated]), col = "red", pch = 3, cex = 2)
points(0, mean(y[idx_control]), col = "red", pch = 4, cex = 2)
dev.off()

png(file.path(script_dir, "demo_2.png"), width = 1200, height = 900, res = 150)
plot(
  X[idx_treated],
  y[idx_treated],
  col = rgb(0, 0, 1, 0.8),
  ylab = "y",
  xlab = "x",
  xlim = c(min(X), max(X)),
  ylim = c(min(y), max(y))
)
points(mean(X[idx_treated]), mean(y[idx_treated]), pch = 3, cex = 2, col = "red")
lines(
  c(min(X[idx_treated]), max(X[idx_treated])),
  c(beta_hat_treated) * (c(min(X[idx_treated]), max(X[idx_treated])) - mean(X[idx_treated])) + mean(y[idx_treated]),
  col = rgb(0, 0, 255, 150, maxColorValue = 255)
)
lines(
  c(min(X[idx_control]), max(X[idx_control])),
  c(beta_hat_control) * (c(min(X[idx_control]), max(X[idx_control])) - mean(X[idx_control])) + mean(y[idx_control]),
  col = rgb(0, 0, 0, 150, maxColorValue = 255)
)
points(X[idx_control], y[idx_control], col = rgb(0.5, 0.5, 0.5, 0.8))
points(mean(X[idx_control]), mean(y[idx_control]), pch = 4, cex = 2, col = "red")
abline(v = mean(X), col = rgb(0.5, 0.5, 0, 0.9), lty = "dashed")

df_treated <- data.frame(
  x = X[idx_treated] - mean(X[idx_treated]),
  y = y[idx_treated] - mean(y[idx_treated])
)
mod <- lm(y ~ x, data = df_treated)
newx <- seq(min(X[idx_treated]) - mean(X[idx_treated]), max(X[idx_treated]) - mean(X[idx_treated]), length.out = 100)
preds <- predict(mod, newdata = data.frame(x = newx), interval = "confidence")
polygon(
  c(rev(newx), newx) + mean(X[idx_treated]),
  c(rev(preds[, 3]), preds[, 2]) + mean(y[idx_treated]),
  col = rgb(0, 0, 1, 0.2),
  border = NA
)

df_control <- data.frame(
  x = X[idx_control] - mean(X[idx_control]),
  y = y[idx_control] - mean(y[idx_control])
)
mod <- lm(y ~ x, data = df_control)
newx <- seq(min(X[idx_control]) - mean(X[idx_control]), max(X[idx_control]) - mean(X[idx_control]), length.out = 100)
preds <- predict(mod, newdata = data.frame(x = newx), interval = "confidence")
polygon(
  c(rev(newx), newx) + mean(X[idx_control]),
  c(rev(preds[, 3]), preds[, 2]) + mean(y[idx_control]),
  col = rgb(0.5, 0.5, 0.5, 0.2),
  border = NA
)
legend("topleft", legend = c("Control", "Treatment"), 
       col = c(rgb(0,0,0,150, maxColorValue=255), 
               rgb(0,0,255,150, maxColorValue=255)), 
       pch = 1, bty = "n")
dev.off()
