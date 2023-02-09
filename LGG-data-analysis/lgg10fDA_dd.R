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

predAPT_all <- matrix(0, nrow = npat, ncol = 9)
nclust_all <- matrix(0, nrow = K, ncol = 6)
gof_all <- matrix(0, nrow = K, ncol = 2)

wk <- c(0, 40, 100)

registerDoParallel(cores = 10)#alloco solo core necessari

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
    modelpriors$hP0_m0 <- rep(0, ncol(Y_train)); modelpriors$hP0_nu0 <- .10
    modelpriors$hP0_s0 <- ncol(Y_train) + 2; modelpriors$hP0_Lambda0 <- .10

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
                                   PPMx = 1, cohesion = 2, kappa = c(1, 20, 10, 1), sigma = c(0.01, .59, 6),
                                   similarity = 2, consim = 2, similparam = vec_par,
                                   calibration = 2, coardegree = 2, modelpriors,
                                   update_hierarchy = T,
                                   hsp = F, iter = iterations, burn = burnin, thin = thinning,
                                   mhtunepar = c(0.05, 0.05), CC = 5, reuse = 1,
                                   nclu_init = 3), error = function(e){FALSE})
    return(res0)
  }

save(myres0, file = "output/lgg-da/lgg10fDA_feb811.RData")
#load("output/lgg-da/lgg10fDA_feb2.RData")

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
  predAPT_all[currfold, 4:9] <- cbind(myres[,, 1], myres[,, 2])

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
#save(pred_meas, file = "output/lgg-da/pred_meas_lgg.RData")
#save(resPPMX, file = "output/lgg-da/resPPMX_lgg.RData")
#save(clu, file = "output/lgg-da/clu_lgg.RData")

lpml <- list()
for(k in 1:10){
  #currfold <- (vectf[k]:(vectf[k+1]-1))
  res0 <- myres0[[k]]
  vec <- res0$lpml
  lpml[[k]] <- coda::mcmc(matrix(vec, nrow = nout, ncol=1))
}

for(k in 1:10){
  df <- ggmcmc::ggs(lpml[[k]])
  pl <- ggplot2::ggplot(df, aes(x = Iteration, y = value)) +
    geom_line() + theme_classic() + ggtitle(paste0("Fold ", k)) +
    xlab("Iterations") + ylab("lpml")# + ylim(c(-250, -50))
  assign(paste("c", k, sep=""), pl)
}

plpml <- ggpubr::ggarrange(c1, c2, c3, c4, c5, c6, c7, c8, c9, c10, ncol = 2, nrow=5,
                           common.legend = TRUE, legend="bottom")
plpml
#ggsave(plpml, device = "pdf", path = "figs", filename = "lpml_plot.pdf")

for(k in 1:10){
  out_ppmx <- myres0[[k]]
  df <- data.frame(t(out_ppmx$nclu))
  colnames(df) <- c("Treatment 1", "Treatment 2")
  df <- cbind(Index = as.numeric(row.names(df)), df)
  df <- reshape2::melt(df, id.vars="Index")
  colnames(df) <- c("Index", "Treatment", "value")
  clu_tp <- ggplot2::ggplot(df, aes(x = Index, y = value, col = Treatment)) +
    geom_line() + theme_classic() + ggtitle(paste0("Fold ", k)) +
    xlab("Iterations") + ylab("# of clusters") + ylim(10, 30)
  assign(paste("clu_tp", k, sep=""), clu_tp)
}

pnc <- ggpubr::ggarrange(clu_tp1, clu_tp2, clu_tp3, clu_tp4, clu_tp5, clu_tp6,
                         clu_tp7, clu_tp8, clu_tp9, clu_tp10, ncol = 2, nrow=5,
                         common.legend = TRUE, legend="bottom")
pnc

#ggsave(pnc, device = "pdf", path = "figs", filename = "nc_plot.pdf")

for(k in 1:10){
  #k=4
  out_ppmx <- myres0[[k]]

  df_sigma <- rbind(data.frame(table(round(out_ppmx$sigmangg[1,],1))),
                    data.frame(table(round(out_ppmx$sigmangg[2,],1))))
  vec <- c(rep("Treatment 1", dim(table(out_ppmx$sigmangg[1,]))), rep("Treatment 2", dim(table(out_ppmx$sigmangg[2,]))))
  df_sigma <- cbind(df_sigma, Treatment = as.factor(vec))
  df <- df_sigma %>%
    group_by(Treatment)
  colnames(df) <- c("sigma", "frequency", "Treatment")

  sigma <- ggplot(df, aes(x=sigma, y=frequency/nout, fill = Treatment)) +
    geom_bar(stat = "identity", position = "dodge", width = 0.1) +
    ylab("proportion") + xlab(expression(sigma)) + ggtitle(paste0("Fold ", k)) +
    ylim(0, 0.75)
  assign(paste("sp", k, sep=""), sigma)

  df_kappa <- rbind(data.frame(table(round(out_ppmx$kappangg[1,],1))),
                    data.frame(table(round(out_ppmx$kappangg[2,],1))))
  vec <- c(rep("Treatment 1", dim(table(out_ppmx$kappangg[1,]))), rep("Treatment 2", dim(table(out_ppmx$kappangg[2,]))))
  df_kappa <- cbind(df_kappa, Treatment = as.factor(vec))
  df <- df_kappa %>%
    group_by(Treatment)
  colnames(df) <- c("kappa", "frequency", "Treatment")

  kappa <- ggplot(df, aes(x=kappa, y=frequency/nout, fill = Treatment)) +
    geom_bar(stat = "identity", position = "dodge", width = 0.1) +
    ylab("proportion") + xlab(expression(kappa)) + ggtitle(paste0("Fold ", k)) +
    ylim(0, .5)
  assign(paste("kp", k, sep=""), kappa)
  #ksp <- ggpubr::ggarrange(sigma, kappa, nrow=1, ncol = 2, legend = "none" )#,
  #                        #common.legend = TRUE, legend="bottom")
  #ksp <- ggpubr::annotate_figure(ksp, top = ggpubr::text_grob(paste("Fold ", k)))
  #assign(paste("ksp", k, sep=""), ksp)
}


pks <- ggpubr::ggarrange(sp1, kp1, sp2, kp2, sp3, kp3, sp4, kp4, sp5, kp5, sp6, kp6,
                         sp7, kp7, sp8, kp8, sp9, kp9, sp10, kp10,
                         ncol = 4, nrow=5,
                         common.legend = TRUE, legend="bottom")
#ggsave(pks, device = "pdf", path = "figs", filename = "ks_plot.pdf")
pks













