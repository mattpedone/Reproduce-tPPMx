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

load("data/LGGdata.rda")
load("~/Dropbox/PHD/study-treatppmx/output/lgg12aprs121.RData")
#load(file = "output/lgg_analysis_24gen_noHS.RData")
matchRTComp <- matchRTComp[sample(1:nrow(matchRTComp), size = nrow(matchRTComp), replace = F),]
trtsgn <- c(matchRTComp[,10]) + 1
#trtsgn <- ifelse(trtsgn == 1, 2, 1)
npat <- length(trtsgn)
#trtsgn <- sample(1:2, npat, replace = TRUE)

K <- 10#numero di fold

wk <- c(0, 40, 100)

Y <- matrix(0, nrow = npat, ncol = max(as.numeric(matchRTComp[,9])))
for(i in 1:nrow(Y)){
  Y[i, as.numeric(matchRTComp[i,9])] <- 1
}

table(matchRTComp[,9:10])
vectf <- c(1, 17, 33, 49, 65, 81, 97, 113, 129, 145, 159)

nout <- 1600
X <- scale(matchRTComp[,16:38])
Z <- scale(matchRTComp[,c(11,13)])

for(k in 1:10){
  out_ppmx <- myres0[[k]]
  df <- data.frame(t(out_ppmx$nclu))
  colnames(df) <- c("Treatment 1", "Treatment 2")
  df <- cbind(Index = as.numeric(row.names(df)), df)
  df <- reshape2::melt(df, id.vars="Index")
  colnames(df) <- c("Index", "Treatment", "value")
  clu_tp <- ggplot2::ggplot(df, aes(x = Index, y = value, col = Treatment)) +
    geom_line() + theme_classic() + ggtitle(paste0("Fold ", k)) +
    xlab("Iterations") + ylab("# of clusters") + ylim(0, 20)
  assign(paste("clu_tp", k, sep=""), clu_tp)
}

pnc <- ggpubr::ggarrange(clu_tp1, clu_tp2, clu_tp3, clu_tp4, clu_tp5, clu_tp6,
                         clu_tp7, clu_tp8, clu_tp9, clu_tp10, ncol = 2, nrow=5,
                         common.legend = TRUE, legend="bottom")
pnc

#ggsave(pnc, device = "pdf", path = "figs", filename = "nc_plot.pdf")

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
    xlab("Iterations") + ylab("lpml") + ylim(c(-200, -50))
  assign(paste("c", k, sep=""), pl)
}

plpml <- ggpubr::ggarrange(c1, c2, c3, c4, c5, c6, c7, c8, c9, c10, ncol = 2, nrow=5,
                           common.legend = TRUE, legend="bottom")
plpml
#ggsave(plpml, device = "pdf", path = "figs", filename = "lpml_plot.pdf")

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
    ylim(0, 0.75)
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

