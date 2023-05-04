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
library(compositions)
library(robCompositions)
library(dirmult)
library(scales)

load("data/LGGdata.rda")
matchRTComp <- matchRTComp[sample(1:nrow(matchRTComp), size = nrow(matchRTComp), replace = F),]
trtsgn <- c(matchRTComp[,10]) + 1
#trtsgn <- ifelse(trtsgn == 1, 2, 1)
npat <- length(trtsgn)
#trtsgn <- sample(1:2, npat, replace = TRUE)

K <- 10#numero di fold

predAPT_all <- matrix(0, nrow = npat, ncol = 11)
nclust_all <- matrix(0, nrow = K, ncol = 6)
gof_all <- matrix(0, nrow = K, ncol = 2)

wk <- c(0, 25, 75, 100)

registerDoParallel(cores = 10)#alloco solo core necessari
my_level <- c("Progressive Disease", "Stable Disease", "Partial Remission/Response", "Complete Remission/Response")
tempy <- factor(matchRTComp[,3], levels = my_level)
tempy <- as.integer(tempy)
Y <- matrix(0, nrow = npat, ncol = max(tempy))
for(i in 1:nrow(Y)){
  Y[i, tempy[i]] <- 1
}

table(matchRTComp[,9:10])
vectf <- c(1, 17, 33, 49, 65, 81, 97, 113, 129, 145, 159)

nout <- 2000
X <- scale(matchRTComp[,16:38])
Z <- scale(matchRTComp[,c(11,13)])

myres0 <- foreach(k = 1:10) %dorng%
  {
    currfold <- (vectf[k]:(vectf[k+1]-1))
    X_train <- data.frame(X[-currfold,])
    Z_train <- data.frame(Z[-currfold,])
    Y_train <- data.frame(Y[-currfold,])

    X_test <- data.frame(X[currfold,])
    Z_test <- data.frame(Z[currfold,])
    Y_test <- data.frame(Y[currfold,])

    trtsgn_train <- trtsgn[-currfold]
    trtsgn_test <- trtsgn[currfold]

    modelpriors <- list()
    modelpriors$hP0_m0 <- rep(0, ncol(Y_train)); modelpriors$hP0_nu0 <- 10
    modelpriors$hP0_s0 <- ncol(Y_train) + 2; modelpriors$hP0_Lambda0 <- 10

    #n_aux <- 5 # auxiliary variable for Neal's Algorithm 8
    vec_par <- c(0.0, 1.0, .5, 1.0, 2.0, 2.0, 0.1)
    #double m0=0.0, s20=10.0, v=.5, k0=1.0, nu0=2.0, n0 = 2.0;
    iterations <- 12000
    burnin <- 2000
    thinning <- 5

    nout <- (iterations-burnin)/thinning
    predAPT <- c()

    res0 <- tryCatch(expr = ppmxct(y = data.matrix(Y_train), X = data.frame(X_train),
                                   Xpred = data.frame(X_test), Z = data.frame(Z_train),
                                   Zpred = data.frame(Z_test), asstreat = trtsgn_train, #treatment,
                                   PPMx = 1, cohesion = 2, kappa = 1,
                                   similarity = 2, consim = 2, similparam = vec_par,
                                   calibration = 2, coardegree = 2, modelpriors = modelpriors,
                                   update_hierarchy = T,
                                   hsp = F, iter = iterations, burn = burnin, thin = thinning,
                                   mhtunepar = c(0.05, 0.05), CC = 5, reuse = 1,
                                   nclu_init = 3), error = function(e){FALSE})
    return(res0)
  }

#save(myres0, file = "output/lgg-da/lgg10fDA_feb9.RData")
#load("output/lgg-da/lgg10fDA_feb9.RData")

for(k in 1:K){
  currfold <- (vectf[k]:(vectf[k+1]-1))
  res0 <- myres0[[k]]
  #number of a cluster, mean, binder &varinf ----
  mc <- apply(res0$nclu, 1, mean)
  trt <- trtsgn[-currfold]#simdata$trtsgn[[k]][1:124]
  num_treat <- table(trt)

  cls1 <- t(as.matrix(res0$label[[1]]))[,c(1:num_treat[1])]
  psm1 <- comp.psm(cls1)
  mc_b1 <- minbinder.ext(psm1)
  mc_vi1 <- minVI(psm1)

  cls2 <- t(as.matrix(res0$label[[2]]))[,c(1:num_treat[2])]
  psm2 <- comp.psm(cls2)
  mc_b2 <- minbinder.ext(psm2)
  mc_vi2 <- minVI(psm2)

  mc_b <- c(max(mc_b1$cl), max(mc_b2$cl))
  mc_vi <- c(max(mc_vi1$cl), max(mc_vi2$cl))

  myres <- apply(res0$pipred, c(1,2,3), median, na.rm=TRUE)
  myclu <- rbind(mc, mc_b, mc_vi)
  myfit <- c(res0$WAIC, mean(res0$lpml))
  A1 <- myres[,, 1]%*%wk
  A2 <- myres[,, 2]%*%wk
  predAPT_all[currfold, 1] <- A1
  predAPT_all[currfold, 2] <- A2
  myt <- as.numeric(A1 < A2) + 1
  predAPT_all[currfold, 3] <- myt
  predAPT_all[currfold, 4:11] <- cbind(myres[,, 1], myres[,, 2])

  nclust_all[k,] <- c(t(myclu))
  gof_all[k,] <- myfit

  #myprob <- simdata$prob[[k]]
}

#NPC
npc_tf <- function(output, trtsgn, myoutot){
  #K <- dim(output)[3]
  n <- dim(output)[1]
  myctut <- matrix(0, nrow = 3, ncol = 3)
  myctutSum <- NULL
  mycurdata <- output
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
  myctut <- str1 * diag(3)
  myctutSum <- sum(str1 * diag(3))
  #res <- cbind(myctutSum)
  return(myctutSum)
}


PPMXCUT <- c()
temp <- matrix(0, nrow = npat, ncol = 6)
for(k in 1:K){
  currfold <- (vectf[k]:(vectf[k+1]-1))
  myoutot <- as.numeric(matchRTComp[currfold,9])#simdata$yord[[k]][131:158,]
  trtsgn_test <- trtsgn[currfold]#simdata$trtsgn[[k]][131:158]
  temp[currfold,] <- predAPT_all[currfold, 4:9]
}
NPC <- npc_tf(temp, trtsgn, as.numeric(matchRTComp[,9]))
#NPC <- c(round(mean(PPMXCUT), 4), round(sd(PPMXCUT), 4))

#ESM
# non ho definito come respondent anche i partial responent
myoutot <- as.numeric(matchRTComp[,9])#simdata$yord[[k]][131:158,]
mytab <- cbind(myass = predAPT_all[,3], rndass = trtsgn, resp = as.numeric(myoutot>2))
pred1 <- subset(mytab, mytab[,1]==1)
table1 <- table(pred1[,3],pred1[,2])
pred2 <- subset(mytab, mytab[,1]==2)
table2 <- table(pred2[,3], pred2[,2])
p1 <- sum(table1)/(sum(table1)+sum(table2))
p2 <- sum(table2)/(sum(table1)+sum(table2))

if(length(table1) == 4){
  crt1 <- table1[2,1]/sum(table1[,1])
}
if(length(table1) < 4){
  crt1 <- as.numeric(row.names(table1))
}

if(length(table2) == 4){
  crt2 <- table2[2,2]/sum(table2[,2])
}
if(length(table2) < 4){
  crt2 <- as.numeric(row.names(table2))
}
ESM <- c(crt1*p1 + crt2*p2 - sum(as.numeric(myoutot>2))/npat)
### summary meaures

#FIT
WAIC <- c(mean(gof_all[,1]), sd(gof_all[,1]))
lpml <- c(mean(gof_all[,2]), sd(gof_all[,2]))
#plot(gof_all[,2], type="l")

#results
resPPMX <- rbind(WAIC, lpml)
colnames(resPPMX) <- c("mean", "sd")

cluPPMX <- nclust_all[,-c(3,4)]
clu <- apply(cluPPMX, 2, mean)
clu <- rbind(clu, apply(cluPPMX, 2, sd))
colnames(clu) <- c("avg # trt 1", "avg # trt 2", "VI trt 1", "VI trt 2")

pred_meas <- c(NPC, ESM)
pred_meas
save(pred_meas, file = "output/lgg-da/pred_meas_lgg_4r.RData")
