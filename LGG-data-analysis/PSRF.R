rm(list=ls())
set.seed(121)
library(treatppmx)
library(parallel)
library(doParallel)
library(doRNG)
library(mcclust)
library(mcclust.ext)
library(ggplot2)
library(reshape2)
library(plotly)
library(dplyr)
library(ggtern)
library(coda)
load("data/LGGdata.rda")
matchRTComp <- matchRTComp[sample(1:nrow(matchRTComp), size = nrow(matchRTComp), replace = F),]
trtsgn <- c(matchRTComp[,10]) + 1
#trtsgn <- ifelse(trtsgn == 1, 2, 1)
npat <- length(trtsgn)
#trtsgn <- sample(1:2, npat, replace = TRUE)

K <- 5#numero di fold

predAPT_all <- matrix(0, nrow = npat, ncol = 9)
nclust_all <- matrix(0, nrow = K, ncol = 6)
gof_all <- matrix(0, nrow = K, ncol = 2)
wk <- c(0, 40, 100)
registerDoParallel(cores = 5)#alloco solo core necessari
Y <- matrix(0, nrow = npat, ncol = max(as.numeric(matchRTComp[,9])))
for(i in 1:nrow(Y)){
  Y[i, as.numeric(matchRTComp[i,9])] <- 1
}
table(matchRTComp[,9:10])
vectf <- c(1, 17, 33, 49, 65, 81, 97, 113, 129, 145, 159)
#load("/home/matt/Dropbox/PHD/study-treatppmx/output/lgg12aprs121.RData")
nout <- 2000
X <- scale(matchRTComp[,16:38])
Z <- scale(matchRTComp[,c(11,13)])

myres0 <- foreach(k = 1:K) %dorng%
  {
    modelpriors <- list()
    modelpriors$hP0_m0 <- rep(0, ncol(Y)); modelpriors$hP0_nu0 <- 10
    modelpriors$hP0_s0 <- ncol(Y) + 2; modelpriors$hP0_Lambda0 <- 10

    #n_aux <- 5 # auxiliary variable for Neal's Algorithm 8
    vec_par <- c(0.0, 1.0, .5, 1.0, 2.0, 2.0, 0.1)
    #double m0=0.0, s20=10.0, v=.5, k0=1.0, nu0=2.0, n0 = 2.0;
    iterations <- 12000
    burnin <- 2000
    thinning <- 5

    nout <- (iterations-burnin)/thinning
    predAPT <- c()

    res0 <- tryCatch(expr = ppmxct(y = data.matrix(Y), X = data.frame(X),
                                   Xpred = data.frame(X[1:2,]), Z = data.frame(Z),
                                   Zpred = data.frame(Z[1:2,]), asstreat = trtsgn, #treatment,
                                   PPMx = 1, cohesion = 2, kappa = 1,
                                   similarity = 2, consim = 2, similparam = vec_par,
                                   calibration = 2, coardegree = 2, modelpriors,
                                   update_hierarchy = F,
                                   hsp = T, iter = iterations, burn = burnin, thin = thinning,
                                   mhtunepar = c(0.05, 0.05), CC = 5, reuse = 1,
                                   nclu_init = 3), error = function(e){FALSE})
    return(res0)
  }

beta11 <- beta12 <- beta13 <- beta21 <- beta22 <- beta23 <- vector("list", 5L)

for(k in 1:K){
  beta11[[k]] <- myres0[[k]]$beta[1,1,]
  beta12[[k]] <- myres0[[k]]$beta[1,2,]
  beta13[[k]] <- myres0[[k]]$beta[1,3,]
  beta21[[k]] <- myres0[[k]]$beta[2,1,]
  beta22[[k]] <- myres0[[k]]$beta[2,2,]
  beta23[[k]] <- myres0[[k]]$beta[2,3,]
}

psrf <- matrix(0, 2, 3)
ml <- as.mcmc.list(lapply(as.data.frame(beta11), mcmc))
mcmc_data <- mcmc.list(ml)
psrf[1,1] <- gelman.diag(mcmc_data)$psrf[1]

ml <- as.mcmc.list(lapply(as.data.frame(beta12), mcmc))
mcmc_data <- mcmc.list(ml)
psrf[1,2] <- gelman.diag(mcmc_data)$psrf[1]

ml <- as.mcmc.list(lapply(as.data.frame(beta13), mcmc))
mcmc_data <- mcmc.list(ml)
psrf[1,3] <- gelman.diag(mcmc_data)$psrf[1]

ml <- as.mcmc.list(lapply(as.data.frame(beta21), mcmc))
mcmc_data <- mcmc.list(ml)
psrf[2,1] <- gelman.diag(mcmc_data)$psrf[1]

ml <- as.mcmc.list(lapply(as.data.frame(beta22), mcmc))
mcmc_data <- mcmc.list(ml)
psrf[2,2] <- gelman.diag(mcmc_data)$psrf[1]

ml <- as.mcmc.list(lapply(as.data.frame(beta23), mcmc))
mcmc_data <- mcmc.list(ml)
psrf[2,3] <- gelman.diag(mcmc_data)$psrf[1]

psrf
