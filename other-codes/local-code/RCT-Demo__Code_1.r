rm(list = ls())

set.seed(1016)
n = 50
n_T = 15
n_C = n - n_T
X = rnorm(n)
e = 0.5 * rnorm(n)

idx_treated = sample(n, n_T, replace=FALSE)
idx_control = setdiff(1:n,idx_treated)

T = rep(0, n)
T[idx_treated] = 1

y = T + X + e

ATE_unadj = mean(y[idx_treated])-mean(y[idx_control])

X_tmp = X[idx_treated] - mean(X[idx_treated])
y_tmp = y[idx_treated] - mean(y[idx_treated])
beta_hat_treated = solve(t(X_tmp)%*%X_tmp)%*%t(X_tmp)%*%y_tmp
X_tmp = X[idx_control] - mean(X[idx_control])
y_tmp = y[idx_control] - mean(y[idx_control])
beta_hat_control = solve(t(X_tmp)%*%X_tmp)%*%t(X_tmp)%*%y_tmp
ATE_adj = (mean(y[idx_treated])-(mean(X[idx_treated]) - mean(X))*beta_hat_treated) - (mean(y[idx_control])--(mean(X[idx_control]) - mean(X))*beta_hat_control)


# plot
par(mfrow=c(1,2))
plot(T[idx_treated], y[idx_treated], xlim=c(-0.5,1.5), ylim=c(min(y), max(y)), ylab="Outcomes", xlab="T", col=rgb(0, 0, 255, 150, maxColorValue=255), main=paste("ATE hat =", round(ATE_unadj,2)))
points(T[idx_control], y[idx_control], col=rgb(0, 0, 0, 150, maxColorValue=255))
points(1, mean(y[idx_treated]), col="red", pch=3, cex=2)
points(0, mean(y[idx_control]), col="red", pch=4, cex=2)

plot(X[idx_treated], y[idx_treated], col=rgb(0, 0, 1, 0.8), main=paste("ATE hat =",round(ATE_adj,2)), ylab="Outcomes", xlab="X", xlim=c(min(X), max(X)), ylim=c(min(y), max(y)))
points(mean(X[idx_treated]), mean(y[idx_treated]),pch=3, cex=2, col="red")
lines(c(min(X[idx_treated]), max(X[idx_treated])), c(beta_hat_treated)*(c(min(X[idx_treated]), max(X[idx_treated])) - mean(X[idx_treated])) + mean(y[idx_treated]), col=rgb(0, 0, 255, 150, maxColorValue=255))
lines(c(min(X[idx_control]), max(X[idx_control])), c(beta_hat_control)*(c(min(X[idx_control]), max(X[idx_control])) - mean(X[idx_control])) + mean(y[idx_control]), col=rgb(0, 0, 0, 150, maxColorValue=255))
points(X[idx_control], y[idx_control], col=rgb(0.5, 0.5, 0.5, 0.8))
points(mean(X[idx_control]), mean(y[idx_control]),pch=4, cex=2, col="red")
# vertical line: X mean
abline(v=mean(X), col=rgb(0.5, 0.5, 0, 0.9), lty='dashed')

# 95% predict interval
# prediction interval: treated
df_treated <- data.frame(x = X[idx_treated]-mean(X[idx_treated]),
                         y = y[idx_treated]-mean(y[idx_treated]))
mod <- lm(y ~ x, data = df_treated)
newx <- seq(min(X[idx_treated]) - mean(X[idx_treated]), max(X[idx_treated]) - mean(X[idx_treated]), length.out=100)
preds <- predict(mod, newdata = data.frame(x=newx), interval = 'confidence')
polygon(c(rev(newx), newx)+mean(X[idx_treated]), c(rev(preds[ ,3]), preds[ ,2])+mean(y[idx_treated]), col = rgb(0, 0, 1, 0.2), border = NA, )
#lines(newx+mean(X[idx_treated]), preds[ ,3]+mean(y[idx_treated]), lty = 'dashed', col = 'red')
#lines(newx+mean(X[idx_treated]), preds[ ,2]+mean(y[idx_treated]), lty = 'dashed', col = 'red')

# prediction interval: control
df_control <- data.frame(x = X[idx_control]-mean(X[idx_control]),
                         y = y[idx_control]-mean(y[idx_control]))
mod <- lm(y ~ x, data = df_control)
newx <- seq(min(X[idx_control]) - mean(X[idx_control]), max(X[idx_control]) - mean(X[idx_control]), length.out=100)
preds <- predict(mod, newdata = data.frame(x=newx), interval = 'confidence')
polygon(c(rev(newx), newx)+mean(X[idx_control]), c(rev(preds[ ,3]), preds[ ,2])+mean(y[idx_control]), col = rgb(0.5, 0.5, 0.5, 0.2), border = NA, )
#lines(newx+mean(X[idx_control]), preds[ ,3]+mean(y[idx_control]), lty = 'dashed', col = 'red')
#lines(newx+mean(X[idx_control]), preds[ ,2]+mean(y[idx_control]), lty = 'dashed', col = 'red')

# # vertical line: weighted X mean # 我不知道这算的是什么. 所以注释掉了. 
# # abline(v=mean(X), col=rgb(0.5, 0.5, 0, 0.5))
# X_center = (mean(X[idx_treated])%*%solve(t(X[idx_treated])%*%X[idx_treated]) + mean(X[idx_control])%*%solve(t(X[idx_control])%*%X[idx_control])) %*% solve(solve(t(X[idx_treated])%*%X[idx_treated]) + solve(t(X[idx_control])%*%X[idx_control]))
# abline(v=X_center, col=rgb(0.5, 0.5, 0, 0.5))
