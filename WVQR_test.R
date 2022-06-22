library(randtoolbox)
library(mvtnorm)
library(quantreg)
library(pls)
library(ggplot2)
library(wordspace)
source(paste0(getwd(),"/WVQR2.R"))

set.seed(123)

### Definition of some parameters
n = 1000 
p = 4 
d = 1 # dimension of Y
t = 50 # Sample size of Uniform distribution
w = 0.4 # concentration mass at w
# m = 5  
X = matrix(rnorm(n * p), n, p) # n x p
beta = matrix(rnorm(p * d, mean = 6, sd = 10), p, d) # p x d
res = matrix(rnorm(n * d), n, d)
Y = X %*% beta + res # n x d
nu = rep(1/n, n)
mu = rep(1/t, t)

# X = matrix(rnorm(n * p, mean = 0, sd = 0.1), n, p) # n x p
# res = matrix(rnorm(n * d, mean = 0, sd = 0.1), n, d)
# Y = X %*% beta + res # n x d

# distance = matrix(0, nrow = rep, ncol = 2)

#for (i in 1:rep){
### 1. Generate a spherical uniform distribution.
  
# zerod = rep(0, d)
# Identity = diag(1, d, d)
# r = runif(t, 0, 1)
# v_1 = rmvnorm(t, mean = zerod, sigma = Identity)
# v = normalize.rows(v_1)
# U =  r * v

U = as.matrix(torus(t, dim = d)) # t * d
# U = as.matrix(seq(1/t, 1 - 1/t, length.out = t))

### 2. Generate weights. 
# mean = rep(w, d)
# Sigma = diag(1e-2, d)
# W = as.matrix(rmvnorm(m, mean = mean, sigma = Sigma)) # m * d

### 3. Solve 
wvqr = WVQR2(X = X, Y = Y, U = U, mu = mu, nu = nu, w = w)
sol_wvqr = wvqr$beta
sol_ls = solve(t(X) %*% X) %*% t(X) %*% Y

if (d == 1) {
  qr = rq(Y ~ X, tau = w)
  sol_qr = as.matrix(qr$coefficients[2:(p + 1)])
}

Norm(sol_wvqr - sol_ls, p = 2)
Norm(sol_wvqr - sol_qr, p = 2)

# Results
# beta
# sol_Wvqr
# sol_ls
# sol_qr # If d = 1

# distance[i, 1] = Norm(sol_wvqr - sol_ls, p = 2)
# distance[i, 2] = Norm(sol_wvqr - sol_qr, p = 2)
# 
# }
# 
# pos = 0
# for (i in 1:rep){
#   if (distance[i, 2] - distance[i, 1] > 0){
#     pos = pos + 1}
# }
# 
# ### 4. Distance
# distance
# pos

