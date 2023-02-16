rm(list=ls())
set.seed(121)

library(rstan)
library(parallel)
library(doParallel)
library(doRNG)

#options(mc.cores = parallel::detectCores())
registerDoParallel(cores = 10)
rstan_options(auto_write = TRUE)

loadRData <- function(fileName){
  #loads an RData file, and returns it
  load(fileName)
  get(ls()[ls() != "fileName"])
}

simdata <- loadRData("~/Dropbox/PHD/study-treatppmx/data/scen3a.RData")

npc2 <- function(output, trtsgn, myoutot){
  K <- dim(output)[3]
  n <- dim(output)[1]
  myctut <- array(0, dim = c(3, 3, K))
  myctutSum <- NULL
  for (i in 1:K) {
    mycurdata <- output[, , i]
    mypre <- NULL
    pretrt1 <- apply(mycurdata[, 1:3], 1, which.max)
    pretrt2 <- apply(mycurdata[, 4:6], 1, which.max)
    mypreTall <- cbind(pretrt1, pretrt2)
    for (j in 1:n) {
      mypre[j] <- mypreTall[j, trtsgn[j]]
    }
    sts <- table(mypre, myoutot)
    mysdls <- as.numeric(rownames(sts))
    str1 <- matrix(0, nrow = 3, ncol = 3)
    str1[mysdls, ] <- sts
    myctut[, , i] <- str1 * diag(3)
    myctutSum[i] <- sum(str1 * diag(3))
  }
  res <- cbind(myctutSum)
  return(res)
}

Yl <- simdata$ymat
Xl <- simdata$pred
Zl <- simdata$prog
trtl <- simdata$trtsgn

R <- length(Yl)
n <- dim(Yl[[1]])[1]
J <- dim(Yl[[1]])[2]
p <- dim(Xl[[1]])[2]
q <- dim(Zl[[1]])[2]
Y <- array(0, dim = c(n, J, R))
X <- array(0, dim = c(n, p, R))
Z <- array(0, dim = c(n, q, R))
Xdis <- array(0, dim = c(n, (p + q + 1), R))#(p + p + q +1), R))
trt <- matrix(0, n, R)


for(r in 1:R){
  Y[,,r] <- Yl[[r]]
  X[,,r] <- Xl[[r]]
  Z[,,r] <- Zl[[r]]
  trt[,r] <- trtl[[r]]-1
  Xdis[,,r] <- cbind(rep(1, n), X[,,r], Z[,,r])#, (trt[,r]), (trt[,r]*X[,,r]))
}

K <- 10#repliche
npat_pred <- 28

#params1 <- params2 <- array(npat_pred, J, K)

predAPT_all <- array(0, dim = c(npat_pred, 9, K))
nclust_all <- matrix(0, nrow = K, ncol = 6)
gof_all <- matrix(0, nrow = K, ncol = 2)

wk <- c(0, 40, 100)

res <- foreach(k = 1:K) %dorng%
  {
    X_train1 <- Xdis[1:124, ,k]
    X_train1 <- X_train1[which(trt[1:124,k] == 0),]
    Y_train1 <- Y[1:124, ,k]
    Y_train1 <- Y_train1[which(trt[1:124,k] == 0),]

    X_train2 <- Xdis[1:124, ,k]
    X_train2 <- X_train2[which(trt[1:124,k] == 1),]
    Y_train2 <- Y[1:124, ,k]
    Y_train2 <- Y_train2[which(trt[1:124,k] == 1),]

    X_test <- Xdis[125:152,,k]
    Y_test <- Y[125:152,,k]

    ss_data1 = list(N = dim(X_train1)[1], M = dim(X_test)[1], J = dim(Y_train1)[2],
                    P = dim(X_train1)[2], X = X_train1, Xp = X_test, Y = Y_train1,
                    sd_prior = 1.0, psi = .25)

    ss_data2 = list(N = dim(X_train2)[1], M = dim(X_test)[1], J = dim(Y_train2)[2],
                    P = dim(X_train2)[2], X = X_train2, Xp = X_test, Y = Y_train2,
                    sd_prior = 1.0, psi = .25)

    fit1 <- rstan::stan(file = "/home/matt/Dropbox/MYCODE/Reproduce-tPPMx/R/dmhs-scripts/model.stan",
                        data = ss_data1, cores = 1, iter = 500,
                        chains = 1, verbose = T, warmup = 100, seed = 121,
                        control = list(max_treedepth = 15, adapt_delta = 0.995))#995))

    fit2 <- rstan::stan(file = "/home/matt/Dropbox/MYCODE/Reproduce-tPPMx/R/dmhs-scripts/model.stan",
                        data = ss_data2, cores = 1, iter = 500,
                        chains = 1, verbose = T, warmup = 100, seed = 121,
                        control = list(max_treedepth = 15, adapt_delta = 0.995))

    #check_hmc_diagnostics(fit1)
    params1 <- rstan::extract(fit1)
    #check_hmc_diagnostics(fit2)
    params2 <- rstan::extract(fit2)

    out <- cbind(pp1 = apply(params1$pipred, c(2,3), median), pp2 = apply(params2$pipred, c(2,3), median))
    return(out)
  }

#summ <- rstan::summary(fit)$summary

for(k in 1:K){
  pp1 <- res[[k]][,1:3]
  A1 <- pp1%*%wk
  pp2 <- res[[k]][,4:6]
  A2 <- pp2%*%wk

  predAPT_all[, 1, k] <- A1
  predAPT_all[, 2, k] <- A2
  myt <- as.numeric(A1 < A2) + 1
  predAPT_all[, 3, k] <- myt
  predAPT_all[, 4:9, k] <- cbind(pp1, pp2)

  myprob <- simdata$prob[[k]]
}
mywk1 <- myprob[[1]][125:152,]%*%wk
mywk2 <- myprob[[2]][125:152,]%*%wk
optrt <- as.numeric(mywk1 < mywk2) + 1
utsum <- sum(abs(mywk2 - mywk1))
utdiff <- abs(as.numeric(mywk2 - mywk1))

#MOT
PPMXCT <- c()
for(k in 1:K){
  PPMXCT[k] <-  sum(abs(predAPT_all[, 3, k] - optrt))
}

MOT <- c(round(mean(PPMXCT), 4), round(sd(PPMXCT), 4))

#MTUg
PPMXpp <- c()
for(k in 1:K){
  PPMXpp[k] <- -(2*sum(abs((predAPT_all[, 3, k] - optrt)) * utdiff) - utsum);
}

MTUg <- c(round(mean(PPMXpp/utsum), 4), round(sd(PPMXpp/utsum), 4))

#NPC
PPMXCUT <- c()
for(k in 1:K){
  trtsgn_test <- simdata$trtsgn[[k]][125:152]
  temp <- array(0, dim = c(28, 6, 1))
  temp[,,1] <- predAPT_all[, 4:9,k]
  myoutot <- simdata$yord[[k]][125:152,]
  PPMXCUT[k] <- npc2(temp, trtsgn_test, myoutot)
}
NPC <- c(round(mean(PPMXCUT), 4), round(sd(PPMXCUT), 4))

#results
resPPMX <- rbind(MOT, MTUg, NPC)#, WAIC, lpml)
colnames(resPPMX) <- c("mean", "sd")
resPPMX
