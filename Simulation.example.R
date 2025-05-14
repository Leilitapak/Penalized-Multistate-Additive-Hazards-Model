#########################################################################
###### Penalized Variable Selection in Additive Hazards Regression ######
###### for Multi-State Time-to-Event Data                          ######
#########################################################################
# +---------------------------------------------------------------------+
# |                      Checking Required Files                        |
# +---------------------------------------------------------------------+
# Make sure that the file `ahfun.R` and one of the following files 
# (depending on your operating system) are available in R's working directory:
#
#     `ahfun32.dll` – DLL (Dynamic-Link Library) compiled from `ahfun.c` for Windows 32-bit
#     `ahfun64.dll` – DLL (Dynamic-Link Library) compiled from `ahfun.c` for Windows 64-bit
#     `ahfun.so`    – Shared object (SO) compiled from `ahfun.c` for Linux
# 
# Note: On Linux, you must manually compile `ahfun.c` to create `ahfun.so`.
# Please refer to the README file for instructions.
#
# To check the current working directory, run: getwd()

# +---------------------------------------------------------------------+
# |                      Loading External Library                       |
# +---------------------------------------------------------------------+
dyn.load("ahfun64.dll")     # Uncomment this line if you're using 64-bit Windows
#dyn.load("ahfun32.dll")    # Uncomment if you're using 32-bit Windows
#dyn.load("ahfun.so")       # Uncomment if you're on Linux


# +---------------------------------------------------------------------+
# |                    Importing R Wrapper Functions                    |
# +---------------------------------------------------------------------+
source("ahfun.R")


# +---------------------------------------------------------------------+
# |                  Installing and Loading Packages                    |
# +---------------------------------------------------------------------+
# It is recommended to use R version 3.5.0 or higher to ensure that all
# dependencies of the "mstate" package can be installed and that the code
# functions as intended.

packages <- c("mstate", "MASS")

installed <- packages %in% installed.packages()[, "Package"]
if (any(!installed)) install.packages(packages[!installed])

lapply(packages, library, character.only = TRUE)


# +---------------------------------------------------------------------+
# |                  Defining Simulation Parameters                     |
# +---------------------------------------------------------------------+
# General simulation parameters
N <- 500                                      # Sample size
p <- 50                                       # Number of variables
n.rep <- 200                                  # Number of simulation repetitions


# Defining parameters for Multivariate Normal Predictors
times <- 1:p
rho <- 0.5                                    # Correlation value
H <- abs(outer(times, times, "-"))
C <- 1 * rho^H                                # Correlation matrix
sigma <- matrix(C,p,p)
mu <- rep(0,p)



# Defining effect sizes for each transition (Large/Moderate/Low effects)

## Large effects (for informative variables in transitions)
beta12 <- c(.95,.95,.95,.95,.95,rep(0,p - 5)) # Five informative variables in transition 1 → 2
beta13 <- c(.95,.95,.95,rep(0,p - 3))         # Three informative variables in transition 1 → 3
beta23 <- c(.95,.95,.95,.95,.95,rep(0,p - 5)) # Five informative variables in transition 2 → 3

## Moderate effects
beta12 <- c(.5,.5,.5,.5,.5,rep(0,p - 5))
beta13 <- c(.5,.5,.5,rep(0,p - 3))
beta23 <- c(.5,.5,.5,.5,.5,rep(0,p - 5))

## Low effects
beta12 <- c(.1,.1,.1,.1,.1,rep(0,p - 5))
beta13 <- c(.1,.1,.1,rep(0,p - 3))
beta23 <- c(.1,.1,.1,.1,.1,rep(0,p - 5))


# +---------------------------------------------------------------------+
# |          Initializing Matrices for Estimated Coefficients           |
# +---------------------------------------------------------------------+
# Transition 1 → 2
beta.lasso.12 <- matrix(0, ncol = p, nrow = n.rep)
beta.SICA.12 <- matrix(0, ncol = p, nrow = n.rep)
beta.SCAD.12 <- matrix(0, ncol = p, nrow = n.rep)
beta.MCP.12 <- matrix(0, ncol = p, nrow = n.rep)
beta.Enet.12 <- matrix(0, ncol = p, nrow = n.rep)

# Transition 1 → 3
beta.lasso.13 <- matrix(0, ncol = p, nrow = n.rep)
beta.SICA.13 <- matrix(0, ncol = p, nrow = n.rep)
beta.SCAD.13 <- matrix(0, ncol = p, nrow = n.rep)
beta.MCP.13 <- matrix(0, ncol = p, nrow = n.rep)
beta.Enet.13 <- matrix(0, ncol = p, nrow = n.rep)

# Transition 2 → 3
beta.lasso.23 <- matrix(0, ncol = p, nrow = n.rep)
beta.SICA.23 <- matrix(0, ncol = p, nrow = n.rep)
beta.SCAD.23 <- matrix(0, ncol = p, nrow = n.rep)
beta.MCP.23 <- matrix(0, ncol = p, nrow = n.rep)
beta.Enet.23 <- matrix(0, ncol = p, nrow = n.rep)


##############################################################################
# Example Scenario: N <- 500, p <- 50, n.rep <- 200, rho <- 0.5, Low effects #
##############################################################################
# +---------------------------------------------------------------------+
# |                   Running Simulation Repetitions                    |
# +---------------------------------------------------------------------+

for (kk in 1:n.rep){
  
  set.seed(kk)
  # +-------------------------------------------------------------+
  # |          Generating Multivariate Normal Predictors          |
  # +-------------------------------------------------------------+
  x <- mvrnorm(n=N, mu, sigma)
  col_names <- c(sprintf("X%d", seq(1,dim(x)[2])))
  colnames(x) <- col_names
  
  # +-------------------------------------------------------------+
  # |                       Splitting data                        |
  # +-------------------------------------------------------------+
  N12 <- N-2*N%/%5
  N13 <- N-N12
  
  # +-------------------------------------------------------------+
  # |       Defining Multi-State Data Generation Parameters       |
  # +-------------------------------------------------------------+
  hx12 <- 2+(x[1:N12,] %*% beta12)
  ty12 <- rexp(N12,hx12)
  
  hx13 <- 2+(x[(N12+1):N,] %*% beta13)
  ty13 <- rexp(N13,hx13)
  
  hx23 <- 2+(x[1:N12,] %*% beta23)
  ty23 <- rexp(N12,hx23)
  
  t2 <- rep(0,N)
  t2[1:N12] <- ty12
  t2[(N12+1):N] <- ty13
  s2 <- rep(0,N)
  s2[1:N12] <- rep(1,N12)
  
  t3 <- rep(0,N)
  t3[1:N12] <- ty12+ty23
  t3[(N12+1):N] <- ty13
  tcens3 <- rbinom(n=N, prob = 0.4, size = 1) # Censoring rate = %40
  s3 <- 1-tcens3
  
  dt <- data.frame(illt = t2, ills = s2, dt = t3, ds = s3, x)
  indx <- !is.nan(dt$illt) & !is.nan(dt$dt)
  dt2 <- dt[indx,] # Remove rows with missing times
  dim(dt2)
  
  
  # +-------------------------------------------------------------+
  # |         Multi-State Data Preparation (Long Format)          |
  # +-------------------------------------------------------------+
  tmat <- matrix(c(NA,NA,NA,1,NA,NA,2,3,NA), nrow = 3)
  
  longdt <- msprep(time=c(NA,"illt","dt"),
                   status=c(NA,"ills","ds"),
                   keep = col_names, data=dt2,trans=tmat)
  
  msbmt <- longdt # Store long-format data
  
  
  # +-------------------------------------------------------------+
  # |             Extracting Transition-Specific Data             |
  # +-------------------------------------------------------------+
  datatrans12<-msbmt[ which(msbmt$trans==1 ),7:ncol(msbmt)]
  datatrans13<-msbmt[ which(msbmt$trans==2 ),7:ncol(msbmt)]
  datatrans23<-msbmt[ which(msbmt$trans==3 ),7:ncol(msbmt)]
  
  
  # Transition 1 → 2
  z.12 <- as.matrix(datatrans12[,-c(1,2)] )
  cft.12 <- datatrans12$time 		        # Censored failure time
  del.12 <- datatrans12$status 	        # Failure status
  
  
  # Transition 1 → 3
  z.13 <- as.matrix(datatrans13[,-c(1,2)] )
  cft.13 <- datatrans13$time 		        # Censored failure time
  del.13 <- datatrans13$status	        # Failure status
  
  
  # Transition 2 → 3
  z.23 <- as.matrix(datatrans23[,-c(1,2)] )
  cft.23 <- datatrans23$time 		        # Censored failure time
  del.23 <- datatrans23$status	        # Failure status
  
  
  # +-------------------------------------------------------------+
  # |                                                             |
  # |     Fitting Penalized Multi-State Additive Hazards Model    |
  # |                                                             |
  # +-------------------------------------------------------------+
  
  # +-------------------------------------------------------------+
  # |                       Transition 1 → 2                      |
  # +-------------------------------------------------------------+
  ########### Lasso
  ans <- cd(cft.12, del.12, z.12, "Lasso")           # Fitting Lasso using coordinate descent algorithm
  lam <- ans$lam                                     # Lambda values (denoted as γ in the paper)
  sol <- ans$sol                                     # Coefficient solutions over lambda grid
  ans <- cv(cft.12, del.12, z.12, "Lasso", lam =lam) # Performing cross-validation to select optimal lambda
  i <- ans$i                                         # Index of optimal lambda
  beta.lasso.12[kk,] <- sol[,i]                      # Storing selected coefficients
  
  
  
  ########### SICA
  ans <- cd(cft.12, del.12, z.12, "SICA")            # Fitting SICA using coordinate descent algorithm
  lam <- ans$lam                                     # Lambda values (denoted as γ in the paper)
  sol <- ans$sol                                     # Coefficient solutions over lambda grid
  ans <- cv(cft.12, del.12, z.12, "SICA", lam =lam)  # Performing cross-validation to select optimal lambda and a
  i <- ans$i                                         # Index of optimal lambda
  j <- ans$j                                         # Index of optimal a
  beta.SICA.12[kk,] <- sol[, i, j]                   # Storing estimated coefficients
  
  
  
  ########### SCAD
  ans <- cd(cft.12, del.12, z.12, "SCAD")            # Fitting SCAD using coordinate descent algorithm
  lam <- ans$lam                                     # Lambda values (denoted as γ in the paper)
  sol <- ans$sol                                     # Coefficient solutions over lambda grid
  ans <- cv(cft.12, del.12, z.12, "SCAD", lam =lam)  # Performing cross-validation to select optimal lambda
  i <- ans$i                                         # Index of optimal lambda
  beta.SCAD.12[kk,] <- sol[,i]                       # Storing estimated coefficients
  
  
  
  ########### MCP
  ans <- cd(cft.12, del.12, z.12, "MCP")             # Fitting MCP using coordinate descent algorithm
  lam <- ans$lam                                     # Lambda values (denoted as γ in the paper)
  sol <- ans$sol                                     # Coefficient solutions over lambda grid
  ans <- cv(cft.12, del.12, z.12, "MCP", lam =lam)   # Performing cross-validation to select optimal lambda
  i <- ans$i                                         # Index of optimal lambda
  beta.MCP.12[kk,] <- sol[,i]                        # Storing estimated coefficients
  
  
  
  ########### Enet
  ans <- cd(cft.12, del.12, z.12, "Enet")            # Fitting Enet using coordinate descent algorithm
  lam <- ans$lam                                     # Lambda values (denoted as γ in the paper)
  sol <- ans$sol                                     # Coefficient solutions over lambda grid
  ans <- cv(cft.12, del.12, z.12, "Enet", lam =lam)  # Performing cross-validation to select optimal lambda and a
  i <- ans$i                                         # Index of optimal lambda
  j <- ans$j                                         # Index of optimal a
  beta.Enet.12[kk,] <- sol[, i, j]                   # Storing estimated coefficients
  
  
  
  # +-------------------------------------------------------------+
  # |                       Transition 1 → 3                      |
  # +-------------------------------------------------------------+
  ########### Lasso
  ans <- cd(cft.13, del.13, z.13, "Lasso")
  sol <- ans$sol
  lam <- ans$lam
  ans <- cv(cft.13, del.13, z.13, "Lasso", lam =lam)
  i <- ans$i
  beta.lasso.13[kk,] <- sol[,i]
  
  
  
  ########### SICA
  ans <- cd(cft.13, del.13, z.13, "SICA")
  sol <- ans$sol
  lam <- ans$lam
  ans <- cv(cft.13, del.13, z.13, "SICA", lam =lam)
  i <- ans$i
  j <- ans$j
  beta.SICA.13[kk,] <- sol[, i, j]
  
  
  
  ########### SCAD
  ans <- cd(cft.13, del.13, z.13, "SCAD")
  sol <- ans$sol
  lam <- ans$lam
  ans <- cv(cft.13, del.13, z.13, "SCAD", lam =lam)
  i <- ans$i
  beta.SCAD.13[kk,] <- sol[,i]
  
  
  
  ########### MCP
  ans <- cd(cft.13, del.13, z.13, "MCP")
  sol <- ans$sol
  lam <- ans$lam
  ans <- cv(cft.13, del.13, z.13, "MCP", lam =lam)
  i <- ans$i
  beta.MCP.13[kk,] <- sol[,i]
  
  
  
  ########### Enet
  ans <- cd(cft.13, del.13, z.13, "Enet")
  sol <- ans$sol
  lam <- ans$lam
  ans <- cv(cft.13, del.13, z.13, "Enet", lam =lam)
  i <- ans$i
  j <- ans$j
  beta.Enet.13[kk,] <- sol[, i, j]
  
  
  
  
  # +-------------------------------------------------------------+
  # |                       Transition 2 → 3                      |
  # +-------------------------------------------------------------+
  ########### Lasso
  ans <- cd(cft.23, del.23, z.23, "Lasso")
  sol <- ans$sol
  lam <- ans$lam
  ans <- cv(cft.23, del.23, z.23, "Lasso", lam =lam)
  i <- ans$i
  beta.lasso.23[kk,] <- sol[,i]
  
  
  
  ########### SICA
  ans <- cd(cft.23, del.23, z.23, "SICA")
  sol <- ans$sol
  lam <- ans$lam
  ans <- cv(cft.23, del.23, z.23, "SICA", lam =lam)
  i <- ans$i
  j <- ans$j
  beta.SICA.23[kk,] <- sol[, i, j]
  
  
  
  ########### SCAD
  ans <- cd(cft.23, del.23, z.23, "SCAD")
  sol <- ans$sol
  lam <- ans$lam
  ans <- cv(cft.23, del.23, z.23, "SCAD", lam =lam)
  i <- ans$i
  beta.SCAD.23[kk,] <- sol[,i]
  
  
  
  ########### MCP
  ans <- cd(cft.23, del.23, z.23, "MCP")
  sol <- ans$sol
  lam <- ans$lam
  ans <- cv(cft.23, del.23, z.23, "MCP", lam =lam)
  i <- ans$i
  beta.MCP.23[kk,] <- sol[,i]
  
  
  
  ########### Enet
  ans <- cd(cft.23, del.23, z.23, "Enet")
  sol <- ans$sol
  lam <- ans$lam
  ans <- cv(cft.23, del.23, z.23, "Enet", lam =lam)
  i <- ans$i
  j <- ans$j
  beta.Enet.23[kk,] <- sol[, i, j]
  
  
  print(kk)
  
}



# +---------------------------------------------------------------------+
# |    Identifying Selected Variables Based on Non-Zero Coefficients    |
# |    for Each Transition                                              |
# +---------------------------------------------------------------------+
# Transition 1 → 2
beta.lasso.12.freq <- matrix(0, ncol = p, nrow = n.rep)
beta.Enet.12.freq <- matrix(0, ncol = p, nrow = n.rep)
beta.SICA.12.freq <- matrix(0, ncol = p, nrow = n.rep)
beta.SCAD.12.freq <- matrix(0, ncol = p, nrow = n.rep)
beta.MCP.12.freq <- matrix(0, ncol = p, nrow = n.rep)

# Transition 1 → 3
beta.lasso.13.freq <- matrix(0, ncol = p, nrow = n.rep)
beta.Enet.13.freq <- matrix(0, ncol = p, nrow = n.rep)
beta.SICA.13.freq <- matrix(0, ncol = p, nrow = n.rep)
beta.SCAD.13.freq <- matrix(0, ncol = p, nrow = n.rep)
beta.MCP.13.freq <- matrix(0, ncol = p, nrow = n.rep)

# Transition 2 → 3
beta.lasso.23.freq <- matrix(0, ncol = p, nrow = n.rep)
beta.Enet.23.freq <- matrix(0, ncol = p, nrow = n.rep)
beta.SICA.23.freq <- matrix(0, ncol = p, nrow = n.rep)
beta.SCAD.23.freq <- matrix(0, ncol = p, nrow = n.rep)
beta.MCP.23.freq <- matrix(0, ncol = p, nrow = n.rep)



for (jj in 1:n.rep){
  # Transition 1 → 2
  beta.lasso.12.freq[jj,] <- ifelse(beta.lasso.12[jj,]!=0,1,0)
  beta.SICA.12.freq[jj,] <- ifelse(beta.SICA.12[jj,]!=0,1,0)
  beta.SCAD.12.freq[jj,] <- ifelse(beta.SCAD.12[jj,]!=0,1,0)
  beta.MCP.12.freq[jj,] <- ifelse(beta.MCP.12[jj,]!=0,1,0)
  beta.Enet.12.freq[jj,] <- ifelse(beta.Enet.12[jj,]!=0,1,0)
  
  # Transition 1 → 3
  beta.lasso.13.freq[jj,] <- ifelse(beta.lasso.13[jj,]!=0,1,0)
  beta.SICA.13.freq[jj,] <- ifelse(beta.SICA.13[jj,]!=0,1,0)
  beta.SCAD.13.freq[jj,] <- ifelse(beta.SCAD.13[jj,]!=0,1,0)
  beta.MCP.13.freq[jj,] <- ifelse(beta.MCP.13[jj,]!=0,1,0)
  beta.Enet.13.freq[jj,] <- ifelse(beta.Enet.13[jj,]!=0,1,0)
  
  # Transition 2 → 3
  beta.lasso.23.freq[jj,] <- ifelse(beta.lasso.23[jj,]!=0,1,0)
  beta.SICA.23.freq[jj,] <- ifelse(beta.SICA.23[jj,]!=0,1,0)
  beta.SCAD.23.freq[jj,] <- ifelse(beta.SCAD.23[jj,]!=0,1,0)
  beta.MCP.23.freq[jj,] <- ifelse(beta.MCP.23[jj,]!=0,1,0)
  beta.Enet.23.freq[jj,] <- ifelse(beta.Enet.23[jj,]!=0,1,0)
  
}


# +---------------------------------------------------------------------+
# |          Computing Performance Metrics for Each Transition          |
# +---------------------------------------------------------------------+
## Initialize Matrices to Store Performance Metrics

# Transition 1 → 2
TP12 <- matrix(0,ncol = 5, nrow = n.rep)                            # True Positives
colnames(TP12) <- c("LASSO", "SICA" , "SCAD", "MCP", "Enet")

FP12 <- matrix(0,ncol = 5, nrow = n.rep)                            # False Positives
colnames(FP12) <- c("LASSO", "SICA" , "SCAD", "MCP", "Enet")

TN12 <- matrix(0,ncol = 5, nrow = n.rep)                            # True Negatives
colnames(TN12) <- c("LASSO", "SICA" , "SCAD", "MCP", "Enet")

FN12 <- matrix(0,ncol = 5, nrow = n.rep)                            # False Negatives
colnames(FN12) <- c("LASSO", "SICA" , "SCAD", "MCP", "Enet")

L2.Loss12 <-  matrix(0,ncol = 5, nrow = n.rep)                      # L2 Loss
colnames(L2.Loss12) <- c("LASSO", "SICA" , "SCAD", "MCP", "Enet")

L1.Loss12 <-  matrix(0,ncol = 5, nrow = n.rep)                      # L1 Loss
colnames(L1.Loss12) <- c("LASSO", "SICA" , "SCAD", "MCP", "Enet")



# Transition 1 → 3
TP13 <- matrix(0,ncol = 5, nrow = n.rep)
colnames(TP13) <- c("LASSO", "SICA" , "SCAD", "MCP", "Enet")

FP13 <- matrix(0,ncol = 5, nrow = n.rep)
colnames(FP13) <- c("LASSO", "SICA" , "SCAD", "MCP", "Enet")

TN13 <- matrix(0,ncol = 5, nrow = n.rep)
colnames(TN13) <- c("LASSO", "SICA" , "SCAD", "MCP", "Enet")

FN13 <- matrix(0,ncol = 5, nrow = n.rep)
colnames(FN13) <- c("LASSO", "SICA" , "SCAD", "MCP", "Enet")

L2.Loss13 <-  matrix(0,ncol = 5, nrow = n.rep)
colnames(L2.Loss13) <- c("LASSO", "SICA" , "SCAD", "MCP", "Enet")

L1.Loss13 <-  matrix(0,ncol = 5, nrow = n.rep)
colnames(L1.Loss13) <- c("LASSO", "SICA" , "SCAD", "MCP", "Enet")



# Transition 2 → 3
TP23 <- matrix(0,ncol = 5, nrow = n.rep)
colnames(TP23) <- c("LASSO", "SICA" , "SCAD", "MCP", "Enet")

FP23 <- matrix(0,ncol = 5, nrow = n.rep)
colnames(FP23) <- c("LASSO", "SICA" , "SCAD", "MCP", "Enet")

TN23 <- matrix(0,ncol = 5, nrow = n.rep)
colnames(TN23) <- c("LASSO", "SICA" , "SCAD", "MCP", "Enet")

FN23 <- matrix(0,ncol = 5, nrow = n.rep)
colnames(FN23) <- c("LASSO", "SICA" , "SCAD", "MCP", "Enet")

L2.Loss23 <-  matrix(0,ncol = 5, nrow = n.rep)
colnames(L2.Loss23) <- c("LASSO", "SICA" , "SCAD", "MCP", "Enet")

L1.Loss23 <-  matrix(0,ncol = 5, nrow = n.rep)
colnames(L1.Loss23) <- c("LASSO", "SICA" , "SCAD", "MCP", "Enet")





## Number of Informative (Nonzero) and Non-informative (Zero) Variables per Transition

# Transition 1 → 2
infor1   <- sum((beta12!=0))
n.infor1 <- sum((beta12==0))

# Transition 1 → 3
infor2 <- sum((beta13!=0))
n.infor2 <- sum((beta13==0))

# Transition 2 → 3
infor3 <- sum((beta23!=0))
n.infor3 <- sum((beta23==0))



for (jj in 1:n.rep){
  
  # Transition 1 → 2
  TP12[jj,1] <- sum(ifelse((beta.lasso.12[jj,]!=0)&(beta12!=0),1,0)) / infor1
  TP12[jj,2] <- sum(ifelse((beta.SICA.12[jj,]!=0)&(beta12!=0),1,0)) / infor1
  TP12[jj,3] <- sum(ifelse((beta.SCAD.12[jj,]!=0)&(beta12!=0),1,0)) / infor1
  TP12[jj,4] <- sum(ifelse((beta.MCP.12[jj,]!=0)&(beta12!=0),1,0)) / infor1
  TP12[jj,5] <- sum(ifelse((beta.Enet.12[jj,]!=0)&(beta12!=0),1,0)) / infor1
  
  FP12[jj,1] <- sum(ifelse((beta.lasso.12[jj,]!=0)&(beta12==0),1,0)) / n.infor1
  FP12[jj,2] <- sum(ifelse((beta.SICA.12[jj,]!=0)&(beta12==0),1,0)) / n.infor1
  FP12[jj,3] <- sum(ifelse((beta.SCAD.12[jj,]!=0)&(beta12==0),1,0)) / n.infor1
  FP12[jj,4] <- sum(ifelse((beta.MCP.12[jj,]!=0)&(beta12==0),1,0)) / n.infor1
  FP12[jj,5] <- sum(ifelse((beta.Enet.12[jj,]!=0)&(beta12==0),1,0)) / n.infor1
  
  TN12[jj,1] <- sum(ifelse((beta.lasso.12[jj,]==0)&(beta12==0),1,0)) / n.infor1
  TN12[jj,2] <- sum(ifelse((beta.SICA.12[jj,]==0)&(beta12==0),1,0)) / n.infor1
  TN12[jj,3] <- sum(ifelse((beta.SCAD.12[jj,]==0)&(beta12==0),1,0)) / n.infor1
  TN12[jj,4] <- sum(ifelse((beta.MCP.12[jj,]==0)&(beta12==0),1,0)) / n.infor1
  TN12[jj,5] <- sum(ifelse((beta.Enet.12[jj,]==0)&(beta12==0),1,0)) / n.infor1
  
  FN12[jj,1] <- sum(ifelse((beta.lasso.12[jj,]==0)&(beta12!=0),1,0)) / infor1
  FN12[jj,2] <- sum(ifelse((beta.SICA.12[jj,]==0)&(beta12!=0),1,0)) / infor1
  FN12[jj,3] <- sum(ifelse((beta.SCAD.12[jj,]==0)&(beta12!=0),1,0)) / infor1
  FN12[jj,4] <- sum(ifelse((beta.MCP.12[jj,]==0)&(beta12!=0),1,0)) / infor1
  FN12[jj,5] <- sum(ifelse((beta.Enet.12[jj,]==0)&(beta12!=0),1,0)) / infor1
  
  L2.Loss12[jj,1] <- sqrt(sum(beta.lasso.12[jj,]- beta12)^2)
  L2.Loss12[jj,2] <- sqrt(sum(beta.SICA.12[jj,]- beta12)^2)
  L2.Loss12[jj,3] <- sqrt(sum(beta.SCAD.12[jj,]- beta12)^2)
  L2.Loss12[jj,4] <- sqrt(sum(beta.MCP.12[jj,]- beta12)^2)
  L2.Loss12[jj,5] <- sqrt(sum(beta.Enet.12[jj,]- beta12)^2)
  
  L1.Loss12[jj,1] <-sum(abs(beta.lasso.12[jj,]- beta12))
  L1.Loss12[jj,2] <-sum(abs(beta.SICA.12[jj,]- beta12))
  L1.Loss12[jj,3] <-sum(abs(beta.SCAD.12[jj,]- beta12))
  L1.Loss12[jj,4] <-sum(abs(beta.MCP.12[jj,]- beta12))
  L1.Loss12[jj,5] <-sum(abs(beta.Enet.12[jj,]- beta12))
  
  
  
  # Transition 1 → 3
  TP13[jj,1] <- sum(ifelse((beta.lasso.13[jj,]!=0)&(beta13!=0),1,0)) / infor2
  TP13[jj,2] <- sum(ifelse((beta.SICA.13[jj,]!=0)&(beta13!=0),1,0)) / infor2
  TP13[jj,3] <- sum(ifelse((beta.SCAD.13[jj,]!=0)&(beta13!=0),1,0)) / infor2
  TP13[jj,4] <- sum(ifelse((beta.MCP.13[jj,]!=0)&(beta13!=0),1,0)) / infor2
  TP13[jj,5] <- sum(ifelse((beta.Enet.13[jj,]!=0)&(beta13!=0),1,0)) / infor2
  
  FP13[jj,1] <- sum(ifelse((beta.lasso.13[jj,]!=0)&(beta13==0),1,0)) / n.infor2
  FP13[jj,2] <- sum(ifelse((beta.SICA.13[jj,]!=0)&(beta13==0),1,0)) / n.infor2
  FP13[jj,3] <- sum(ifelse((beta.SCAD.13[jj,]!=0)&(beta13==0),1,0)) / n.infor2
  FP13[jj,4] <- sum(ifelse((beta.MCP.13[jj,]!=0)&(beta13==0),1,0)) / n.infor2
  FP13[jj,5] <- sum(ifelse((beta.Enet.13[jj,]!=0)&(beta13==0),1,0)) / n.infor2
  
  TN13[jj,1] <- sum(ifelse((beta.lasso.13[jj,]==0)&(beta13==0),1,0)) / n.infor2
  TN13[jj,2] <- sum(ifelse((beta.SICA.13[jj,]==0)&(beta13==0),1,0)) / n.infor2
  TN13[jj,3] <- sum(ifelse((beta.SCAD.13[jj,]==0)&(beta13==0),1,0)) / n.infor2
  TN13[jj,4] <- sum(ifelse((beta.MCP.13[jj,]==0)&(beta13==0),1,0)) / n.infor2
  TN13[jj,5] <- sum(ifelse((beta.Enet.13[jj,]==0)&(beta13==0),1,0)) / n.infor2
  
  FN13[jj,1] <- sum(ifelse((beta.lasso.13[jj,]==0)&(beta13!=0),1,0)) / infor2
  FN13[jj,2] <- sum(ifelse((beta.SICA.13[jj,]==0)&(beta13!=0),1,0)) / infor2
  FN13[jj,3] <- sum(ifelse((beta.SCAD.13[jj,]==0)&(beta13!=0),1,0)) / infor2
  FN13[jj,4] <- sum(ifelse((beta.MCP.13[jj,]==0)&(beta13!=0),1,0)) / infor2
  FN13[jj,5] <- sum(ifelse((beta.Enet.13[jj,]==0)&(beta13!=0),1,0)) / infor2
  
  L2.Loss13[jj,1] <- sqrt(sum(beta.lasso.13[jj,]- beta13)^2)
  L2.Loss13[jj,2] <- sqrt(sum(beta.SICA.13[jj,]- beta13)^2)
  L2.Loss13[jj,3] <- sqrt(sum(beta.SCAD.13[jj,]- beta13)^2)
  L2.Loss13[jj,4] <- sqrt(sum(beta.MCP.13[jj,]- beta13)^2)
  L2.Loss13[jj,5] <- sqrt(sum(beta.Enet.13[jj,]- beta13)^2)
  
  L1.Loss13[jj,1] <-sum(abs(beta.lasso.13[jj,]- beta13))
  L1.Loss13[jj,2] <-sum(abs(beta.SICA.13[jj,]- beta13))
  L1.Loss13[jj,3] <-sum(abs(beta.SCAD.13[jj,]- beta13))
  L1.Loss13[jj,4] <-sum(abs(beta.MCP.13[jj,]- beta13))
  L1.Loss13[jj,5] <-sum(abs(beta.Enet.13[jj,]- beta13))
  
  
  
  # Transition 2 → 3
  TP23[jj,1] <- sum(ifelse((beta.lasso.23[jj,]!=0)&(beta23!=0),1,0)) / infor3
  TP23[jj,2] <- sum(ifelse((beta.SICA.23[jj,]!=0)&(beta23!=0),1,0)) / infor3
  TP23[jj,3] <- sum(ifelse((beta.SCAD.23[jj,]!=0)&(beta23!=0),1,0)) / infor3
  TP23[jj,4] <- sum(ifelse((beta.MCP.23[jj,]!=0)&(beta23!=0),1,0)) / infor3
  TP23[jj,5] <- sum(ifelse((beta.Enet.23[jj,]!=0)&(beta23!=0),1,0)) / infor3
  
  FP23[jj,1] <- sum(ifelse((beta.lasso.23[jj,]!=0)&(beta23==0),1,0)) / n.infor3
  FP23[jj,2] <- sum(ifelse((beta.SICA.23[jj,]!=0)&(beta23==0),1,0)) / n.infor3
  FP23[jj,3] <- sum(ifelse((beta.SCAD.23[jj,]!=0)&(beta23==0),1,0)) / n.infor3
  FP23[jj,4] <- sum(ifelse((beta.MCP.23[jj,]!=0)&(beta23==0),1,0)) / n.infor3
  FP23[jj,5] <- sum(ifelse((beta.Enet.23[jj,]!=0)&(beta23==0),1,0)) / n.infor3
  
  TN23[jj,1] <- sum(ifelse((beta.lasso.23[jj,]==0)&(beta23==0),1,0)) / n.infor3
  TN23[jj,2] <- sum(ifelse((beta.SICA.23[jj,]==0)&(beta23==0),1,0)) / n.infor3
  TN23[jj,3] <- sum(ifelse((beta.SCAD.23[jj,]==0)&(beta23==0),1,0)) / n.infor3
  TN23[jj,4] <- sum(ifelse((beta.MCP.23[jj,]==0)&(beta23==0),1,0)) / n.infor3
  TN23[jj,5] <- sum(ifelse((beta.Enet.23[jj,]==0)&(beta23==0),1,0)) / n.infor3
  
  FN23[jj,1] <- sum(ifelse((beta.lasso.23[jj,]==0)&(beta23!=0),1,0)) / infor3
  FN23[jj,2] <- sum(ifelse((beta.SICA.23[jj,]==0)&(beta23!=0),1,0)) / infor3
  FN23[jj,3] <- sum(ifelse((beta.SCAD.23[jj,]==0)&(beta23!=0),1,0)) / infor3
  FN23[jj,4] <- sum(ifelse((beta.MCP.23[jj,]==0)&(beta23!=0),1,0)) / infor3
  FN23[jj,5] <- sum(ifelse((beta.Enet.23[jj,]==0)&(beta23!=0),1,0)) / infor3
  
  L2.Loss23[jj,1] <- sqrt(sum(beta.lasso.23[jj,]- beta23)^2)
  L2.Loss23[jj,2] <- sqrt(sum(beta.SICA.23[jj,]- beta23)^2)
  L2.Loss23[jj,3] <- sqrt(sum(beta.SCAD.23[jj,]- beta23)^2)
  L2.Loss23[jj,4] <- sqrt(sum(beta.MCP.23[jj,]- beta23)^2)
  L2.Loss23[jj,5] <- sqrt(sum(beta.Enet.23[jj,]- beta23)^2)
  
  L1.Loss23[jj,1] <-sum(abs(beta.lasso.23[jj,]- beta23))
  L1.Loss23[jj,2] <-sum(abs(beta.SICA.23[jj,]- beta23))
  L1.Loss23[jj,3] <-sum(abs(beta.SCAD.23[jj,]- beta23))
  L1.Loss23[jj,4] <-sum(abs(beta.MCP.23[jj,]- beta23))
  L1.Loss23[jj,5] <-sum(abs(beta.Enet.23[jj,]- beta23))
  
}


# +---------------------------------------------------------------------+
# |    Counting the Number of Selected Variables for Each Transition    |
# +---------------------------------------------------------------------+
# Transition 1 → 2
selected.p12 <- data.frame(rowSums(beta.lasso.12.freq),rowSums(beta.SICA.12.freq),
                           rowSums(beta.SCAD.12.freq), rowSums(beta.MCP.12.freq),rowSums(beta.Enet.12.freq))
colnames(selected.p12) <-  c("LASSO", "SICA" , "SCAD", "MCP", "Enet")



# Transition 1 → 3
selected.p13 <- data.frame(rowSums(beta.lasso.13.freq),rowSums(beta.SICA.13.freq),
                           rowSums(beta.SCAD.13.freq), rowSums(beta.MCP.13.freq),rowSums(beta.Enet.13.freq))
colnames(selected.p13) <-  c("LASSO", "SICA" , "SCAD", "MCP", "Enet")



# Transition 2 → 3
selected.p23 <- data.frame(rowSums(beta.lasso.23.freq),rowSums(beta.SICA.23.freq),
                           rowSums(beta.SCAD.23.freq), rowSums(beta.MCP.23.freq),rowSums(beta.Enet.23.freq))
colnames(selected.p23) <-  c("LASSO", "SICA" , "SCAD", "MCP", "Enet")



# +----------------------------------------------------------------------------------------------+
# |  Computing the Mean and Standard Deviation (SD) of Performance Criteria for Each Transition  |
# +----------------------------------------------------------------------------------------------+
# Transition 1 → 2
T12 <- t(rbind(colMeans(L2.Loss12),colMeans(L1.Loss12),
               colMeans(selected.p12),colMeans(TP12),colMeans(FP12)))
colnames(T12) <- c("L2","L1","#S","TP","FP")                  # Number of selected variables: #S

T12.s <- t(rbind(apply(L2.Loss12,2,sd),apply(L1.Loss12,2,sd),
                 apply(selected.p12,2,sd),apply(TP12,2,sd),apply(FP12,2,sd)))
colnames(T12.s) <- c("L2","L1","#S","TP","FP")

T12;T12.s




# Transition 1 → 3
T13 <- t(rbind(colMeans(L2.Loss13),colMeans(L1.Loss13),
               colMeans(selected.p13),colMeans(TP13),colMeans(FP13)))
colnames(T13) <- c("L2","L1","#S","TP","FP")

T13.s <- t(rbind(apply(L2.Loss13,2,sd),apply(L1.Loss13,2,sd),
                 apply(selected.p13,2,sd),apply(TP13,2,sd),apply(FP13,2,sd)))
colnames(T13.s) <- c("L2","L1","#S","TP","FP")

T13;T13.s




# Transition 2 → 3
T23 <- t(rbind(colMeans(L2.Loss23),colMeans(L1.Loss23),
               colMeans(selected.p23),colMeans(TP23),colMeans(FP23)))
colnames(T23) <- c("L2","L1","#S","TP","FP")

T23.s <- t(rbind(apply(L2.Loss23,2,sd),apply(L1.Loss23,2,sd),
                 apply(selected.p23,2,sd),apply(TP23,2,sd),apply(FP23,2,sd)))
colnames(T23.s) <- c("L2","L1","#S","TP","FP")

T23;T23.s
