###################################################
### chunk number 1: foo
###################################################
options(keep.source = TRUE, width = 60)
foo <- packageDescription("mcmc")


###################################################
### chunk number 2: frequentist
###################################################
foo <- read.table("logit.txt", header = TRUE)
out <- glm(y ~ x1 + x2 + x3 + x4, data = foo,
    family = binomial())
summary(out)


###################################################
### chunk number 3: log.unnormalized.posterior
###################################################
x <- foo
x$y <- NULL
x <- as.matrix(x)
x <- cbind(1, x)
dimnames(x) <- NULL

y <- foo$y

lupost <- function(beta, x, y) {
    eta <- x %*% beta
    p <- 1 / (1 + exp(- eta))   # note: works even when eta is Inf or -Inf
    logl <- sum(log(p[y == 1])) + sum(log(1 - p[y == 0]))
    return(logl + sum(dnorm(beta, 0, 2, log = TRUE)))
}


###################################################
### chunk number 4: metropolis-try-1
###################################################
library(mcmc)
set.seed(42)    # to get reproducible results
beta.init <- as.numeric(coefficients(out))
out <- metrop(lupost, beta.init, 1e3, x = x, y = y)
names(out)
out$accept


###################################################
### chunk number 5: metropolis-try-2
###################################################
out <- metrop(out, scale = 0.1, x = x, y = y)
out$accept
out <- metrop(out, scale = 0.3, x = x, y = y)
out$accept
out <- metrop(out, scale = 0.5, x = x, y = y)
out$accept
out <- metrop(out, scale = 0.4, x = x, y = y)
out$accept


###################################################
### chunk number 6: metropolis-try-3
###################################################
out <- metrop(out, nbatch = 1e4, x = x, y = y)
out$accept
out$time


###################################################
### chunk number 7: fig1too
###################################################
plot(ts(out$batch))


###################################################
### chunk number 8: fig1
###################################################
plot(ts(out$batch))


###################################################
### chunk number 9: fig2too
###################################################
acf(out$batch)


###################################################
### chunk number 10: fig2
###################################################
acf(out$batch)


###################################################
### chunk number 11: metropolis-try-4
###################################################
out <- metrop(out, nbatch = 1e2, blen = 100,
    outfun = function(z, ...) c(z, z^2), x = x, y = y)
out$accept
out$time


###################################################
### chunk number 12: metropolis-batch
###################################################
apply(out$batch, 2, mean)


###################################################
### chunk number 13: metropolis-batch-too
###################################################
foo <- apply(out$batch, 2, mean)
mu <- foo[1:5]
sigmasq <- foo[6:10] - mu^2
mu
sigmasq


###################################################
### chunk number 14: metropolis-mcse-mu
###################################################
mu.mcse <- apply(out$batch[ , 1:5], 2, sd) / sqrt(out$nbatch)
mu.mcse


###################################################
### chunk number 15: metropolis-mcse-sigmasq
###################################################
u <- out$batch[ , 1:5]
v <- out$batch[ , 6:10]
ubar <- apply(u, 2, mean)
vbar <- apply(v, 2, mean)
deltau <- sweep(u, 2, ubar)
deltav <- sweep(v, 2, vbar)
foo <- sweep(deltav, 2, ubar, "*")
sigmasq.mcse <- sqrt(apply((deltav - 2 * foo)^2, 2, mean) / out$nbatch)
sigmasq.mcse


###################################################
### chunk number 16: metropolis-mcse-sigmasq-too
###################################################
sqrt(mean(((v[ , 2] - vbar[2]) - 2 * ubar[2] * (u[ , 2] - ubar[2]))^2) /
    out$nbatch)

stop()


###################################################
### chunk number 17: metropolis-mcse-foo
###################################################
apply(out$batch[ , 6:10], 2, mean) - mu^2
apply(sweep(out$batch[ , 6:10], 2, mu^2), 2, mean)


###################################################
### chunk number 18: metropolis-mcse-sigmasq
###################################################
sigmasq.mcse <- apply(sweep(out$batch[ , 6:10], 2, mu^2), 2, sd) /
    sqrt(out$nbatch)
sigmasq.mcse


###################################################
### chunk number 19: metropolis-mcse-sigma
###################################################
sigma <- sqrt(sigmasq)
sigma.mcse <- sigmasq.mcse / (2 * sigma)
sigma
sigma.mcse


###################################################
### chunk number 20: metropolis-try-5
###################################################
out <- metrop(out, nbatch = 5e2, blen = 400, x = x, y = y)
out$accept
out$time
foo <- apply(out$batch, 2, mean)
mu <- foo[1:5]
sigmasq <- foo[6:10] - mu^2
mu
sigmasq
mu.mcse <- apply(out$batch[ , 1:5], 2, sd) / sqrt(out$nbatch)
mu.mcse
sigmasq.mcse <- apply(sweep(out$batch[ , 6:10], 2, mu^2), 2, sd) /
    sqrt(out$nbatch)
sigmasq.mcse
sigma <- sqrt(sigmasq)
sigma.mcse <- sigmasq.mcse / (2 * sigma)
sigma
sigma.mcse


###################################################
### chunk number 21: tab1
###################################################
foo <- rbind(mu, mu.mcse)
dimnames(foo) <- list(c("estimate", "MCSE"),
    c("constant", paste("\$x_", 1:4, "\$", sep = "")))
library(xtable)
print(xtable(foo, digits = rep(4, 6),
    align = c("l", rep("c", 5))), floating = FALSE)


###################################################
### chunk number 22: tab1
###################################################
foo <- rbind(sigmasq, sigmasq.mcse)
dimnames(foo) <- list(c("estimate", "MCSE"),
    c("constant", paste("\$x_", 1:4, "\$", sep = "")))
library(xtable)
print(xtable(foo, digits = rep(4, 6),
    align = c("l", rep("c", 5))), floating = FALSE)


###################################################
### chunk number 23: tab1
###################################################
foo <- rbind(sigma, sigma.mcse)
dimnames(foo) <- list(c("estimate", "MCSE"),
    c("constant", paste("\$x_", 1:4, "\$", sep = "")))
library(xtable)
print(xtable(foo, digits = rep(4, 6),
    align = c("l", rep("c", 5))), floating = FALSE)


###################################################
### chunk number 24: time
###################################################
cat(out$time[1], "\n")


