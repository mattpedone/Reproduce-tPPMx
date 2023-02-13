##### ------ SCENARIO 4S-misspec ------ #####
rm(list = ls())
set.seed(121)

library(treatppmx)
library(parallel)
library(doParallel)
library(mclust)
library(mcclust)
library(mcclust.ext)
library(doRNG)

K <- 50#repliche
npat_pred <- 30

predAPT_all <- array(0, dim = c(npat_pred, 9, K))
nclust_all <- matrix(0, nrow = K, ncol = 6)
gof_all <- matrix(0, nrow = K, ncol = 2)

wk <- c(0, 40, 100)

cor_all <- parallel::detectCores() - 1#cores to be allocated
registerDoParallel(cores = K)

simdata <- list()
for (k in 1:K) {
  simdata[[k]] <- treatppmx::genmech_clu2(npred = 10)
}

myres0 <- foreach(k = 1:K) %dopar%
  {
    pred_mis <- sample(1:10, 3)
    prog_mis <- sample(1:2, 1)
    X_train <- data.frame(cbind(simdata[[k]]$pred[1:170, -pred_mis], simdata[[k]]$prog[1:170, prog_mis]))
    Z_train <- data.frame(cbind(simdata[[k]]$prog[1:170, -prog_mis], simdata[[k]]$pred[1:170, pred_mis]))
    Y_train <- data.frame(simdata[[k]]$Y[1:170, ])

    X_test <- data.frame(simdata[[k]]$pred[171:200, -pred_mis], simdata[[k]]$prog[171:200, prog_mis])
    Z_test <- data.frame(cbind(simdata[[k]]$prog[171:200, -prog_mis], simdata[[k]]$pred[171:200, pred_mis]))
    Y_test <- data.frame(simdata[[k]]$Y[171:200, ])

    trtsgn_train <- simdata[[k]]$treatment[1:170]
    trtsgn_test <- simdata[[k]]$treatment[171:200]

    modelpriors <- list()
    modelpriors$hP0_m0 <- rep(0, ncol(Y_train)); modelpriors$hP0_nu0 <- 10
    modelpriors$hP0_s0 <- ncol(Y_train) + 2; modelpriors$hP0_Lambda0 <- 10

    vec_par <- c(0.0, 1.0, .5, 1.0, 2.0, 2.0, 0.1)
    #double m0=0.0, s20=10.0, v=.5, k0=1.0, nu0=2.0, n0 = 2.0;
    iterations <- 12000
    burnin <- 2000
    thinning <- 5

    nout <- (iterations - burnin) / thinning
    predAPT <- c()

    res0 <-
      tryCatch(expr = ppmxct(y = data.matrix(Y_train), X = data.frame(X_train),
          Xpred = data.frame(X_test), Z = data.frame(Z_train), Zpred = data.frame(Z_test),
          asstreat = trtsgn_train, PPMx = 0, cohesion = 2, kappa = 1, similarity = 2, consim = 2, similparam = vec_par,
          calibration = 2, coardegree = 2, modelpriors, update_hierarchy = T, hsp = T,
          iter = iterations, burn = burnin, thin = thinning, mhtunepar = c(0.05, 0.05),
          CC = 5, reuse = 1, nclu_init = 10), error = function(e) {FALSE})
    return(res0)
  }

for (k in 1:K) {
  res0 <- myres0[[k]]
  #number of a cluster, mean, binder &varinf ----
  mc <- apply(res0$nclu, 1, mean)
  trt <- simdata[[k]]$treatment[1:170]
  num_treat <- table(trt)

  cls1 <- t(as.matrix(res0$label[[1]]))[, c(1:num_treat[1])]
  psm1 <- comp.psm(cls1)
  mc_b1 <- minbinder.ext(psm1)
  mc_vi1 <- minVI(psm1)

  cls2 <- t(as.matrix(res0$label[[2]]))[, c(1:num_treat[2])]
  psm2 <- comp.psm(cls2)
  mc_b2 <- minbinder.ext(psm2)
  mc_vi2 <- minVI(psm2)

  mc_b <- c(max(mc_b1$cl), max(mc_b2$cl))
  mc_vi <- c(max(mc_vi1$cl), max(mc_vi2$cl))

  myres <- apply(res0$pipred, c(1, 2, 3), median, na.rm = TRUE)
  myclu <- rbind(mc, mc_b, mc_vi)
  myfit <- c(res0$WAIC, mean(res0$lpml))
  A1 <- myres[, , 1] %*% wk
  A2 <- myres[, , 2] %*% wk
  predAPT_all[, 1, k] <- A1
  predAPT_all[, 2, k] <- A2
  myt <- as.numeric(A1 < A2) + 1
  predAPT_all[, 3, k] <- myt
  predAPT_all[, 4:9, k] <- cbind(myres[, , 1], myres[, , 2])

  nclust_all[k, ] <- c(t(myclu))
  gof_all[k, ] <- myfit

  myprob <- simdata[[k]]$prob
}

#MOT
PPMXCT <- c()
for (k in 1:K) {
  myprob <- simdata[[k]]$prob
  mywk1 <- myprob[[1]][171:200, ] %*% wk
  mywk2 <- myprob[[2]][171:200, ] %*% wk
  optrt <- as.numeric(mywk1 < mywk2) + 1
  PPMXCT[k] <-  sum(abs(predAPT_all[, 3, k] - optrt))
}

MOT <- c(round(mean(PPMXCT), 4), round(sd(PPMXCT), 4))

#MTUg
PPMXpp <- c()
for (k in 1:K) {
  myprob <- simdata[[k]]$prob
  mywk1 <- myprob[[1]][171:200, ] %*% wk
  mywk2 <- myprob[[2]][171:200, ] %*% wk
  optrt <- as.numeric(mywk1 < mywk2) + 1
  utsum <- sum(abs(mywk2 - mywk1))
  utdiff <- abs(as.numeric(mywk2 - mywk1))
  PPMXpp[k] <-
    -(2 * sum(abs((
      predAPT_all[, 3, k] - optrt
    )) * utdiff) - utsum)

}

MTUg <-
  c(round(mean(PPMXpp / utsum), 4), round(sd(PPMXpp / utsum), 4))

npc2 <- function(output, trtsgn, myoutot) {
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
    str1[mysdls,] <- sts
    myctut[, , i] <- str1 * diag(3)
    myctutSum[i] <- sum(str1 * diag(3))
  }
  res <- cbind(myctutSum)
  return(res)
}

#NPC
PPMXCUT <- c()
for (k in 1:K) {
  trtsgn_test <- simdata[[k]]$treatment[171:200]
  temp <- array(0, dim = c(30, 6, 1))
  temp[, , 1] <- predAPT_all[, 4:9, k]
  myoutot <- simdata[[k]]$Yord[171:200, 1]
  PPMXCUT[k] <- npc2(temp, trtsgn_test, myoutot)
}
NPC <- c(round(mean(PPMXCUT), 4), round(sd(PPMXCUT), 4))

#NCLU
NC <-
  c(apply(nclust_all[, 1:2], 2, mean), apply(nclust_all[, 1:2], 2, sd))
BI <-
  c(apply(nclust_all[, 3:4], 2, mean), apply(nclust_all[, 3:4], 2, sd))
VI <-
  c(apply(nclust_all[, 5:6], 2, mean), apply(nclust_all[, 5:6], 2, sd))

#FIT
WAIC <- c(mean(gof_all[, 1]), sd(gof_all[, 1]))
lpml <- c(mean(gof_all[, 2]), sd(gof_all[, 2]))

#results
resPPMX <- rbind(MOT, MTUg, NPC, WAIC, lpml)
colnames(resPPMX) <- c("mean", "sd")
resPPMX

# fitted clustering
## treatment 1
ari1 <- c()
for (k in 1:K) {
  #tl <- simdata[[k]]$clu[1:170][trt == 1]
  tl <- simdata[[k]]$clu1[1:85]

  cls1 <- t(as.matrix(myres0[[k]]$label[[1]]))[, c(1:num_treat[1])]
  psm1 <- comp.psm(cls1)
  pl <- minVI(psm1)$cl
  #pl <- minbinder.ext(psm1)$cl

  ari1[k] <- adjustedRandIndex(tl, pl)
}
## treatment 2
ari2 <- c()
for (k in 1:K) {
  #tl <- simdata[[k]]$clu[1:170][trt == 2]
  tl <- simdata[[k]]$clu2[1:85]

  cls2 <- t(as.matrix(myres0[[k]]$label[[2]]))[, c(1:num_treat[2])]
  psm2 <- comp.psm(cls2)
  pl <- minVI(psm2)$cl
  #pl <- minbinder.ext(psm2)$cl

  ari2[k] <- adjustedRandIndex(tl, pl)
}
#VI
mean(ari1)
mean(ari2)
sd(ari1)
sd(ari2)

ari1 <- c()
for (k in 1:K) {
  #tl <- simdata[[k]]$clu[1:170][trt == 1]
  tl <- simdata[[k]]$clu1[1:85]

  cls1 <- t(as.matrix(myres0[[k]]$label[[1]]))[, c(1:num_treat[1])]
  psm1 <- comp.psm(cls1)
  #pl <- minVI(psm1)$cl
  pl <- minbinder.ext(psm1)$cl

  ari1[k] <- adjustedRandIndex(tl, pl)
}
## treatment 2
ari2 <- c()
for (k in 1:K) {
  #tl <- simdata[[k]]$clu[1:170][trt == 2]
  tl <- simdata[[k]]$clu2[1:85]

  cls2 <- t(as.matrix(myres0[[k]]$label[[2]]))[, c(1:num_treat[2])]
  psm2 <- comp.psm(cls2)
  #pl <- minVI(psm2)$cl
  pl <- minbinder.ext(psm2)$cl

  ari2[k] <- adjustedRandIndex(tl, pl)
}
#Binder
mean(ari1)
mean(ari2)

##### ------ SCENARIO 5S-misspec ------ #####
rm(list = ls())
set.seed(121)

library(treatppmx)
library(parallel)
library(doParallel)
library(mclust)
library(mcclust)
library(mcclust.ext)
library(doRNG)

K <- 50#repliche
npat_pred <- 30

predAPT_all <- array(0, dim = c(npat_pred, 9, K))
nclust_all <- matrix(0, nrow = K, ncol = 6)
gof_all <- matrix(0, nrow = K, ncol = 2)

wk <- c(0, 40, 100)

cor_all <- parallel::detectCores() - 1#cores to be allocated
registerDoParallel(cores = K)

simdata <- list()
for (k in 1:K) {
  simdata[[k]] <- treatppmx::genmech_clu2(npred = 20)
}

myres0 <- foreach(k = 1:K) %dopar%
  {
    pred_mis <- sample(1:20, 7)
    prog_mis <- sample(1:2, 1)
    X_train <- data.frame(cbind(simdata[[k]]$pred[1:170, -pred_mis], simdata[[k]]$prog[1:170, prog_mis]))
    Z_train <- data.frame(cbind(simdata[[k]]$prog[1:170, -prog_mis], simdata[[k]]$pred[1:170, pred_mis]))
    Y_train <- data.frame(simdata[[k]]$Y[1:170, ])

    X_test <- data.frame(simdata[[k]]$pred[171:200, -pred_mis], simdata[[k]]$prog[171:200, prog_mis])
    Z_test <- data.frame(cbind(simdata[[k]]$prog[171:200, -prog_mis], simdata[[k]]$pred[171:200, pred_mis]))
    Y_test <- data.frame(simdata[[k]]$Y[171:200, ])

    trtsgn_train <- simdata[[k]]$treatment[1:170]
    trtsgn_test <- simdata[[k]]$treatment[171:200]

    modelpriors <- list()
    modelpriors$hP0_m0 <- rep(0, ncol(Y_train)); modelpriors$hP0_nu0 <- 10
    modelpriors$hP0_s0 <- ncol(Y_train) + 2; modelpriors$hP0_Lambda0 <- 10

    vec_par <- c(0.0, 1.0, .5, 1.0, 2.0, 2.0, 0.1)
    #double m0=0.0, s20=10.0, v=.5, k0=1.0, nu0=2.0, n0 = 2.0;
    iterations <- 12000
    burnin <- 2000
    thinning <- 5

    nout <- (iterations - burnin) / thinning
    predAPT <- c()

    res0 <- tryCatch(expr = ppmxct(y = data.matrix(Y_train), X = data.frame(X_train),
          Xpred = data.frame(X_test), Z = data.frame(Z_train), Zpred = data.frame(Z_test),
          asstreat = trtsgn_train, PPMx = 0, cohesion = 2, kappa = 1, similarity = 2, consim = 2, similparam = vec_par,
          calibration = 2, coardegree = 2, modelpriors, update_hierarchy = T,
          hsp = T, iter = iterations, burn = burnin, thin = thinning,
          mhtunepar = c(0.05, 0.05), CC = 5, reuse = 1, nclu_init = 10),
          error = function(e) {FALSE})
    return(res0)
  }

for (k in 1:K) {
  res0 <- myres0[[k]]
  #number of a cluster, mean, binder &varinf ----
  mc <- apply(res0$nclu, 1, mean)
  trt <- simdata[[k]]$treatment[1:170]
  num_treat <- table(trt)

  cls1 <- t(as.matrix(res0$label[[1]]))[, c(1:num_treat[1])]
  psm1 <- comp.psm(cls1)
  mc_b1 <- minbinder.ext(psm1)
  mc_vi1 <- minVI(psm1)

  cls2 <- t(as.matrix(res0$label[[2]]))[, c(1:num_treat[2])]
  psm2 <- comp.psm(cls2)
  mc_b2 <- minbinder.ext(psm2)
  mc_vi2 <- minVI(psm2)

  mc_b <- c(max(mc_b1$cl), max(mc_b2$cl))
  mc_vi <- c(max(mc_vi1$cl), max(mc_vi2$cl))

  myres <- apply(res0$pipred, c(1, 2, 3), median, na.rm = TRUE)
  myclu <- rbind(mc, mc_b, mc_vi)
  myfit <- c(res0$WAIC, mean(res0$lpml))
  A1 <- myres[, , 1] %*% wk
  A2 <- myres[, , 2] %*% wk
  predAPT_all[, 1, k] <- A1
  predAPT_all[, 2, k] <- A2
  myt <- as.numeric(A1 < A2) + 1
  predAPT_all[, 3, k] <- myt
  predAPT_all[, 4:9, k] <- cbind(myres[, , 1], myres[, , 2])

  nclust_all[k, ] <- c(t(myclu))
  gof_all[k, ] <- myfit

  myprob <- simdata[[k]]$prob
}

#MOT
PPMXCT <- c()
for (k in 1:K) {
  myprob <- simdata[[k]]$prob
  mywk1 <- myprob[[1]][171:200, ] %*% wk
  mywk2 <- myprob[[2]][171:200, ] %*% wk
  optrt <- as.numeric(mywk1 < mywk2) + 1
  PPMXCT[k] <-  sum(abs(predAPT_all[, 3, k] - optrt))
}

MOT <- c(round(mean(PPMXCT), 4), round(sd(PPMXCT), 4))

#MTUg
PPMXpp <- c()
for (k in 1:K) {
  myprob <- simdata[[k]]$prob
  mywk1 <- myprob[[1]][171:200, ] %*% wk
  mywk2 <- myprob[[2]][171:200, ] %*% wk
  optrt <- as.numeric(mywk1 < mywk2) + 1
  utsum <- sum(abs(mywk2 - mywk1))
  utdiff <- abs(as.numeric(mywk2 - mywk1))
  PPMXpp[k] <-
    -(2 * sum(abs((
      predAPT_all[, 3, k] - optrt
    )) * utdiff) - utsum)

}

MTUg <-
  c(round(mean(PPMXpp / utsum), 4), round(sd(PPMXpp / utsum), 4))

npc2 <- function(output, trtsgn, myoutot) {
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
    str1[mysdls,] <- sts
    myctut[, , i] <- str1 * diag(3)
    myctutSum[i] <- sum(str1 * diag(3))
  }
  res <- cbind(myctutSum)
  return(res)
}

#NPC
PPMXCUT <- c()
for (k in 1:K) {
  trtsgn_test <- simdata[[k]]$treatment[171:200]
  temp <- array(0, dim = c(30, 6, 1))
  temp[, , 1] <- predAPT_all[, 4:9, k]
  myoutot <- simdata[[k]]$Yord[171:200, 1]
  PPMXCUT[k] <- npc2(temp, trtsgn_test, myoutot)
}
NPC <- c(round(mean(PPMXCUT), 4), round(sd(PPMXCUT), 4))

#NCLU
NC <-
  c(apply(nclust_all[, 1:2], 2, mean), apply(nclust_all[, 1:2], 2, sd))
BI <-
  c(apply(nclust_all[, 3:4], 2, mean), apply(nclust_all[, 3:4], 2, sd))
VI <-
  c(apply(nclust_all[, 5:6], 2, mean), apply(nclust_all[, 5:6], 2, sd))

#FIT
WAIC <- c(mean(gof_all[, 1]), sd(gof_all[, 1]))
lpml <- c(mean(gof_all[, 2]), sd(gof_all[, 2]))

#results
resPPMX <- rbind(MOT, MTUg, NPC, WAIC, lpml)
colnames(resPPMX) <- c("mean", "sd")
resPPMX

# fitted clustering
## treatment 1
ari1 <- c()
for (k in 1:K) {
  #tl <- simdata[[k]]$clu[1:170][trt == 1]
  tl <- simdata[[k]]$clu1[1:85]

  cls1 <- t(as.matrix(myres0[[k]]$label[[1]]))[, c(1:num_treat[1])]
  psm1 <- comp.psm(cls1)
  pl <- minVI(psm1)$cl
  #pl <- minbinder.ext(psm1)$cl

  ari1[k] <- adjustedRandIndex(tl, pl)
}
## treatment 2
ari2 <- c()
for (k in 1:K) {
  #tl <- simdata[[k]]$clu[1:170][trt == 2]
  tl <- simdata[[k]]$clu2[1:85]

  cls2 <- t(as.matrix(myres0[[k]]$label[[2]]))[, c(1:num_treat[2])]
  psm2 <- comp.psm(cls2)
  pl <- minVI(psm2)$cl
  #pl <- minbinder.ext(psm2)$cl

  ari2[k] <- adjustedRandIndex(tl, pl)
}
#VI
mean(ari1)
mean(ari2)
sd(ari1)
sd(ari2)

ari1 <- c()
for (k in 1:K) {
  #tl <- simdata[[k]]$clu[1:170][trt == 1]
  tl <- simdata[[k]]$clu1[1:85]

  cls1 <- t(as.matrix(myres0[[k]]$label[[1]]))[, c(1:num_treat[1])]
  psm1 <- comp.psm(cls1)
  #pl <- minVI(psm1)$cl
  pl <- minbinder.ext(psm1)$cl

  ari1[k] <- adjustedRandIndex(tl, pl)
}
## treatment 2
ari2 <- c()
for (k in 1:K) {
  #tl <- simdata[[k]]$clu[1:170][trt == 2]
  tl <- simdata[[k]]$clu2[1:85]

  cls2 <- t(as.matrix(myres0[[k]]$label[[2]]))[, c(1:num_treat[2])]
  psm2 <- comp.psm(cls2)
  #pl <- minVI(psm2)$cl
  pl <- minbinder.ext(psm2)$cl

  ari2[k] <- adjustedRandIndex(tl, pl)
}
#Binder
mean(ari1)
mean(ari2)

