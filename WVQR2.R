library(CVXR)
library(clue)
library(randtoolbox)
library(pracma)
library(pls)

WVQR2 = function(X, Y, U, mu, nu, w){
  n	= dim(Y)[1]
  d	= dim(Y)[2]
  p	= dim(X)[2]
  t	= dim(U)[1]
  N = 80000
  sigma = 5e-3
  
  # W: t * d; P: t *d
  if ((n != dim(X)[1]) | (d != dim(U)[2])){stop("wrong dimensions")}
  
  if (d == 1){
    Rw = pnorm(U, mean = w, sd = 0.05)
    
    psi = Variable(n)
    phi = Variable(t)
    b = Variable(p, d)
    
    obj <- sum(psi * nu) + sum(phi * mu) + t(nu) %*% X %*% b %*% t(Rw) %*% mu
    constraints <- list(psi %*% t(rep(1,t)) + rep(1,n) %*% t(phi) + X %*% b %*% t(Rw) >= Y %*% t(Rw))
    problem <- Problem(Minimize(obj), constraints = constraints)
    solution <- solve(problem)
    beta = solution$getValue(b)
    phi = solution$getValue(phi)
    } else{
      # Estimate vector rank
      iter = 5000
      tau = 0.05
      a = as.matrix(rep(0, t))
      L = as.matrix(rep(0, t))
      # estimate f^c(u)
      for (i in 1:iter){
        dE = as.matrix(rep(-1/t, t))
        z = rmvnorm(1, mean = rep(w, d), sigma = diag(sigma, d))
        for (j in 1:t){
          L[j] = (1/2) * norm(z - U[j, ], type = "2")^2 - a[j]
        }
        dE[which.min(L)] = 1 - 1/t
        a = a + tau * dE
      }
      # estimate f(w)
      rw = rmvnorm(N, mean = rep(w, d), sigma = diag(sigma, d)) # N * d
      W1 = matrix(rep(rowNorms(rw), t), nrow = N) # N * t
      U1 = t(matrix(rep(rowNorms(U), N), nrow = t)) # N* t
      WU = rw %*% t(U)
      fw = (1/2) * W1 + (1/2) * U1 - WU - t(matrix(rep(a, N), nrow = t))
      fhat = apply(fw, 1, FUN = min)
      
      # estimate Rw
      RwMat = U %*% t(rw) - matrix(rep(fhat, t), nrow = t)
      RwIndex = max.col(RwMat) # t * 1
      Rw = rw[RwIndex, ]
      
      psi = Variable(n)
      phi = Variable(t)
      b = Variable(p, d)
      
      obj <- sum(psi * nu) + sum(phi * mu) + t(nu) %*% X %*% b %*% t(Rw) %*% mu
      constraints <- list(psi %*% t(rep(1,t)) + rep(1,n) %*% t(phi) + X %*% b %*% t(Rw) >= Y %*% t(Rw))
      problem <- Problem(Minimize(obj), constraints = constraints)
      solution <- solve(problem)
      beta = solution$getValue(b)
      phi = solution$getValue(phi)
    }
  
  return(list(beta = beta, phi = phi))
}