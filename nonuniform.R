library(MASS)
library(mvtnorm)
library(mvmesh)
library(cqrReg)
source(paste0(getwd(),"/VCQR2.R"))

p = 1; d = 1; n = 100; t = 50

X = as.matrix(rnorm(n = n, mean = 0, sd = 1))
e = as.matrix(rnorm(n = n, mean = 0, sd = 0.1))
beta = rnorm(n = 1, mean = 2, sd = 1)
y = X * beta + e

nu = as.matrix(rep(1/n, times = n))
mu = as.matrix(rep(1/t, times = t))
Uni = as.matrix(seq(0, 1, length = t))
Z = as.matrix(rnorm(n = t, mean = 6, sd = 1))

resuni = VQRTp(X = X, Y = y, U = Uni, mu = mu, nu = nu)
resZ = VQRTp(X = X, Y = y, U = Z, mu = mu, nu = nu)
rescqruni = cqr.admm(X = X, y = y, tau = Uni)
rescqrZ = cqr.admm(X = X, y = y, tau = Z)
