library(CVXR)
library(Matrix)
library(matrixcalc)
library(gurobi)
library(lpSolve)

# load_VQRData<-function(folderName){
#   VQRDataPath  <- paste0(getwd(),"/Data/",folderName)
#   Xprov 	<- as.matrix(read.csv(paste0(VQRDataPath,"/X.txt"),sep="\t", header=FALSE))
#   Y		<- as.matrix(read.csv(paste0(VQRDataPath,"/Y.txt"),sep="\t", header=FALSE))
#   n		<- dim(Xprov)[1]
#   if (n != dim(Y)[1]) {stop("Numbers of observations of X and Y do not match")}
#   d		<- dim(Y)[2]
#   r		<- dim(Xprov)[2]+1
#   X		<- cbind(matrix(1,n,1),Xprov)
#   return(list(X=X,Y=Y,n=n,d=d,r=r,step=step))
# }

VQRTp <- function(X, Y, U, mu, nu){

  n	<-dim(Y)[1]
  d	<-dim(Y)[2]
  p	<-dim(X)[2]
  t	<-dim(U)[1]
  if ((n != dim(X)[1]) |( d != dim(U)[2] )) {stop("wrong dimensions")}
  xbar	<- t(nu)%*% X # sample means of X

  pi0 <- Variable(rows = t * n)
  c	<- -matrix(U %*% t(Y), nrow=n * t) # -vec(UY^{T})

  # In kro 1t^T and vec(nu)
  A1 <- kronecker(diag(1,n), matrix(1, 1, t))
  f1 <- matrix(nu, nrow=n)
  # 1n^T kro It and vec(mu)
  A2 <- kronecker(matrix(1, 1, n), diag(1, t))
  f2 <- matrix(mu, nrow=t)
  # X^T Kro U^T
  A3 <- kronecker(t(X),t(U))
  f3 <- matrix(t(U) %*% mu %*% xbar,nrow = d * p)

  objective <- Minimize(t(c) %*% pi0)
  constraints <- list(pi0 >= 0, A1 %*% pi0 == f1, A2 %*% pi0 == f2,
                      A3 %*% pi0 == f3)
  prob <- Problem(objective, constraints)
  solution <- solve(prob)

  optimal_value <- solution$value
  beta <- solution$getDualValue(constraints(prob)[[4]])

  return(list(optimal_value=optimal_value, beta=beta))
}

VQRTp1 <-function(X,Y,U,mu,nu){
  # Monge-Kantorovich (transportation version)
  # computes VQR via primal program
  # (c) by Guillaume Carlier, Victor Chernozhukov and Alfred Galichon

  n	<-dim(Y)[1]
  d	<-dim(Y)[2]
  p	<-dim(X)[2]
  t	<-dim(U)[1]
  if ((n != dim(X)[1]) |( d != dim(U)[2] )) {stop("wrong dimensions")}
  xbar	<- t(nu)%*% X # sample means of X

  c	<- -t(matrix(U %*% t(Y),nrow=n * t)) # -vec(UY^{T})
  # c <- t(-kronecker(Y,U) %*% matrix(diag(1,d),nrow=d*d)) ### TO BE REMOVED
  A1 <- kronecker(sparseMatrix(1:n, 1:n), matrix(1, 1, t)) # In kro 1t^T
  A2 <- kronecker(matrix(1, 1, n), sparseMatrix(1:t, 1:t)) # 1n^T kro It 
  A3 <- kronecker(t(X),t(U)) # X^T Kro U^T
  
  f1<- matrix(t(nu), nrow=n)# vec(nu)
  f2 <- matrix(t(mu),nrow=t) # vec(mu)
  f3 <- matrix(t(U) %*% mu %*% xbar,nrow = p * d) # 
  e <- matrix(1,t*n,1)
  #A <- rbind2(A2, A3)
  #f <- rbind2(f2, f3)
  A <- rbind(A1, A2, A3)
  f <- rbind(f1, f2, f3)
  pi_init <- matrix(mu %*% t(nu),nrow=t*n)

  ############### LP SOLVING PHASE ###############
  result <- gurobi(list(A=A,obj=c,modelsense="min",rhs=f,ub=e,sense="=",start=pi_init), params=NULL )
  if(result$status=="OPTIMAL"){pivec <- result$x; Lvec <- t(result$pi) } else {stop("optimization problem with Gurobi")}

  #############################################

  pi <- matrix(pivec,nrow=t)
  L1vec <-Lvec[1 : n]
  L2vec <- Lvec[(n + 1) : (t + n)]
  L3vec <-Lvec[(n + t + 1):(n + t + d * p)]
  L1		<-matrix(L1vec)
  L2		<-matrix(L2vec)
  L3 <- matrix(L3vec)

  psi		<- -t(L1)
  phi <- - t(L2)
  b	<- -t(L3)
  val		<- matrix.trace(t(U) %*% pi %*% Y)

  #############################################


  return(list(psi = psi, phi = phi, b=b,solve = (result$pi)))
}

VQRTp2 <- function(X, Y, U, mu, nu){
  
  n	<-dim(Y)[1]
  d	<-dim(Y)[2]
  p	<-dim(X)[2]
  t	<-dim(U)[1]
  if ((n != dim(X)[1]) |( d != dim(U)[2] )) {stop("wrong dimensions")}
  xbar	<- t(nu)%*% X # sample means of X
  
  pi0 <- Variable(rows = t * n)
  c	<- -matrix(U %*% t(Y), nrow=n * t) # -vec(UY^{T})
  
  # In kro 1t^T and vec(nu)
  A1 <- kronecker(sparseMatrix(1:n, 1:n), matrix(1, 1, t))
  f1 <- matrix(nu, nrow=n)
  # 1n^T kro It and vec(mu)
  A2 <- kronecker(matrix(1, 1, n), sparseMatrix(1:t, 1:t))
  f2 <- matrix(mu, nrow=t)
  # X^T Kro U^T
  A3 <- kronecker(t(X),t(U))
  f3 <- matrix(t(U) %*% mu %*% xbar,nrow = d * p)
  
  A <- rbind2(A1, A2, A3)
  f <- rbind2(f1, f2, f3)
  
  solve <- lp(direction = "min", objective.in = t(c), 
              const.mat = A, const.dir = "=", 
              const.rhs = f)
  
  return(solve$dual)
}
# 
# ###############################################################
# ################### II. 1-dimensional VQR #####################
# ###############################################################
# ComputeBeta1D <- function( mu,b){
#   m <-dim(mu)[1] # number of rows of mu
#   D <-diag(1,m); # m-order identity matrix
#   for (i in 1:(m-1)) {D[i+1,i] <- (-1)} 
#   beta<-diag(c(1/mu))%*% D %*% b
#   return(beta)
# }