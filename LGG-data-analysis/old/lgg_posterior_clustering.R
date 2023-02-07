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
#load("output/lgg-da/lgg10fDA_jan31.RData")
load("~/Dropbox/PHD/study-treatppmx/output/lgg12aprs121.RData")
matchRTComp <- matchRTComp[sample(1:nrow(matchRTComp), size = nrow(matchRTComp), replace = F),]

mycn <- c(colnames(matchRTComp)[1:10], "ACVRL1", "Src", "HSP70",
          "HER2", "PAI-1", "Smad1", "FOXO3a",
          "PRAS40", "Cyclin", "SF2", "Bad" ,
          "Lck", "Caspase-7_cleavedD198", "Paxillin", "MYH11",
          "14-3-3_epsilon", "Akt", "Caveolin-1", "Rab25",
          "YAP", "RBM15", "Claudin-7", "ER-alpha",
          "C-Raf", "CD31", "Ku80", "Bcl-2", "GSK3-alpha-beta")
colnames(matchRTComp) <- mycn

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

nout <- 160
X <- scale(matchRTComp[,16:38])
Z <- scale(matchRTComp[,c(11,13)])

##### TRATTAMENTO 1 #####

ccm <- matrix(0, nrow=79, ncol = 79)
for(k in 1:10){

  currfold <- (vectf[k]:(vectf[k+1]-1))
  X_train <- data.frame(X[-currfold,])
  Z_train <- data.frame(Z[-currfold,])
  Y_train <- data.frame(Y[-currfold,])

  X_test <- data.frame(X[currfold,])
  Z_test <- data.frame(Z[currfold,])
  Y_test <- data.frame(Y[currfold,])

  trtsgn_train <- trtsgn[-currfold]
  trtsgn_test <- trtsgn[currfold]

  res0 <- myres0[[k]]
  #number of a cluster, mean, binder &varinf ----
  mc <- apply(res0$nclu, 1, mean)
  trt <- trtsgn[-currfold]#simdata$trtsgn[[k]][1:124]
  num_treat <- table(trt)

  cls2 <- t(as.matrix(res0$label[[1]]))[,c(1:num_treat[1])]
  psm2 <- comp.psm(cls2)

  orp <- as.numeric(rownames(X[which(trtsgn == 1),]))
  id <- as.numeric(rownames(X_train[which(trtsgn_train == 1),]))

  #mc_vi2 <- minVI(psm2)
  #clu <- as.factor(mc_vi2$cl)
  co <- c()
  for(i in 1: num_treat[1]){
    co[i] <- which(orp == id[i])
  }
  #ccm[k,co] <- mc_vi2$cl
  ccm[co,co] <- ccm[co,co] + psm2
}

ccmf <- as.matrix(ccm/9)
colnames(ccmf) <- rownames(ccmf) <- orp

hmccm <- heatmap(ccmf)

labels <- minVI(ccmf)$cl
table(labels)
g1 <- orp[which(labels==1)]
g2 <- orp[which(labels==2)]
g3 <- orp[which(labels==3)]
g4 <- orp[which(labels==4)]
g5 <- orp[which(labels==5)]
g6 <- orp[which(labels==6)]
g7 <- orp[which(labels==7)]
g8 <- orp[which(labels==8)]
g9 <- orp[which(labels==9)]
g10 <- orp[which(labels==10)]

reord <- c()
for(i in 1:max(labels)){
  reord <- c(reord, which(labels == i))
}

# Co-occurence plot ----
data1 <- ccmf[,reord]
data <- data1[reord,]
colnames(data) <- rownames(data) <- NULL
#coincidences<-sapply(1:ncol(data), function(i){ colSums(data[,i]==data) })
mC <- melt(data)
colnames(mC) <- c("Var1", "Var2", "Co-occurrence")
c1 <- ggplot(mC, aes(Var1,Var2, fill=`Co-occurrence`)) + geom_raster() +
  scale_fill_continuous(type = "viridis") +
  theme_classic() +
  theme(legend.position="bottom") +
  xlab("Patients") + ylab("Patients") + ggtitle("Treatment 1")#+ ggtitle("Heatmap of averaged co-occurence matrix")
c1
#ggsave(c1, device = "pdf", path = "figs", filename = "avg_coocc.pdf")

#PREDITTIVE
pred_cov <- matchRTComp[, c(16:38)]
pred_cov_g1 <- pred_cov[as.character(g1),]
pred_cov_g2 <- pred_cov[as.character(g2),]
pred_cov_g3 <- pred_cov[as.character(g3),]
pred_cov_g4 <- pred_cov[as.character(g4),]
pred_cov_g5 <- pred_cov[as.character(g5),]
pred_cov_g6 <- pred_cov[as.character(g6),]
pred_cov_g7 <- pred_cov[as.character(g7),]
pred_cov_g8 <- pred_cov[as.character(g8),]
pred_cov_g9 <- pred_cov[as.character(g9),]
pred_cov_g10 <- pred_cov[as.character(g10),]

vartab <- rbind(apply(pred_cov_g1, 2, var), apply(pred_cov_g2, 2, var),
                apply(pred_cov_g3, 2, var), apply(pred_cov_g4, 2, var),
                apply(pred_cov_g5, 2, var), apply(pred_cov_g6, 2, var),
                apply(pred_cov_g7, 2, var), apply(pred_cov_g8, 2, var),
                apply(pred_cov_g9, 2, var),
                apply(pred_cov_g10, 2, var))
#vartab <- rbind(apply(pred_cov_g1, 2, mean), apply(pred_cov_g2, 2, mean), apply(pred_cov_g3, 2, mean), apply(pred_cov_g4, 2, mean), apply(pred_cov_g5, 2, mean))
df <- data.frame(cbind(Group = c(1: 10), vartab))
df$Group <- as.factor(df$Group)
df <- melt(df, id.vars = "Group")
p <- ggplot(df, aes(x=variable, y=value, group=Group)) +
  geom_line(aes(color=Group))+
  geom_point(aes(color=Group)) +
  ylim(0, 1.5) +
  xlab("Predictive biomarkers") + ylab("Variance") +
  theme_classic() +
  theme(axis.text.x=element_text(angle=45, hjust=1)) + ggtitle("Treatment 2")
p
#ggsave(p, device = "pdf", path = "figs", filename = "var_arr_group.pdf")
table(labels)
vartab <- rbind(apply(pred_cov_g1, 2, mean), apply(pred_cov_g2, 2, mean),
                apply(pred_cov_g3, 2, mean), apply(pred_cov_g4, 2, mean),
                apply(pred_cov_g5, 2, mean), apply(pred_cov_g6, 2, mean),
                apply(pred_cov_g7, 2, mean), apply(pred_cov_g8, 2, mean),
                apply(pred_cov_g9, 2, mean), apply(pred_cov_g10, 2, mean))
df <- data.frame(cbind(Group = c(1:10), vartab))
df$Group <- as.factor(df$Group)
df <- melt(df, id.vars = "Group")
#df <- df[which(df$Group != 3),]
#df <- df[which(df$Group != 6),]
#df <- df[which(df$Group != 7),]
#df <- df[which(df$Group != 8),]
#df <- df[which(df$Group != 9),]
p <- ggplot(df, aes(x=variable, y=value, group=Group)) +
  geom_line(aes(color=Group))+
  geom_point(aes(color=Group)) +
  ylim(-1.0, 1.5) +
  xlab("Predictive biomarkers") + ylab("Mean") +
  theme_classic() +
  theme(axis.text.x=element_text(angle=45, hjust=1)) + ggtitle("Treatment 2")
p

##### TRATTAMENTO 2 #####

ccm <- matrix(0, nrow=79, ncol = 79)
for(k in 1:10){

  currfold <- (vectf[k]:(vectf[k+1]-1))
  X_train <- data.frame(X[-currfold,])
  Z_train <- data.frame(Z[-currfold,])
  Y_train <- data.frame(Y[-currfold,])

  X_test <- data.frame(X[currfold,])
  Z_test <- data.frame(Z[currfold,])
  Y_test <- data.frame(Y[currfold,])

  trtsgn_train <- trtsgn[-currfold]
  trtsgn_test <- trtsgn[currfold]

  res0 <- myres0[[k]]
  #number of a cluster, mean, binder &varinf ----
  mc <- apply(res0$nclu, 1, mean)
  trt <- trtsgn[-currfold]#simdata$trtsgn[[k]][1:124]
  num_treat <- table(trt)

  cls2 <- t(as.matrix(res0$label[[2]]))[,c(1:num_treat[2])]
  psm2 <- comp.psm(cls2)

  orp <- as.numeric(rownames(X[which(trtsgn == 2),]))
  id <- as.numeric(rownames(X_train[which(trtsgn_train == 2),]))

  #mc_vi2 <- minVI(psm2)
  #clu <- as.factor(mc_vi2$cl)
  co <- c()
  for(i in 1: num_treat[2]){
    co[i] <- which(orp == id[i])
  }
  #ccm[k,co] <- mc_vi2$cl
  ccm[co,co] <- ccm[co,co] + psm2
}

ccmf <- as.matrix(ccm/9)
colnames(ccmf) <- rownames(ccmf) <- orp

hmccm <- heatmap(ccmf)

labels <- minVI(ccmf)$cl

table(labels)

g1 <- orp[which(labels==1)]
g2 <- orp[which(labels==2)]
g3 <- orp[which(labels==3)]
g4 <- orp[which(labels==4)]
g5 <- orp[which(labels==5)]
g6 <- orp[which(labels==6)]
g7 <- orp[which(labels==7)]
g8 <- orp[which(labels==8)]
g9 <- orp[which(labels==9)]
g10 <- orp[which(labels==10)]

reord <- c()
for(i in 1:max(labels)){
  reord <- c(reord, which(labels == i))
}

# Co-occurence plot ----
data1 <- ccmf[,reord]
data <- data1[reord,]
colnames(data) <- rownames(data) <- NULL
#coincidences<-sapply(1:ncol(data), function(i){ colSums(data[,i]==data) })
mC <- melt(data)
colnames(mC) <- c("Var1", "Var2", "Co-occurrence")
c1 <- ggplot(mC, aes(Var1,Var2, fill=`Co-occurrence`)) + geom_raster() +
  scale_fill_continuous(type = "viridis") +
  theme_classic() +
  theme(legend.position="bottom") +
  xlab("Patients") + ylab("Patients") + ggtitle("Treatment 1")#+ ggtitle("Heatmap of averaged co-occurence matrix")
c1
#ggsave(c1, device = "pdf", path = "figs", filename = "avg_coocc.pdf")

#PREDITTIVE
pred_cov <- matchRTComp[, c(16:38)]
pred_cov_g1 <- pred_cov[as.character(g1),]
pred_cov_g2 <- pred_cov[as.character(g2),]
pred_cov_g3 <- pred_cov[as.character(g3),]
pred_cov_g4 <- pred_cov[as.character(g4),]
pred_cov_g5 <- pred_cov[as.character(g5),]
pred_cov_g6 <- pred_cov[as.character(g6),]
pred_cov_g7 <- pred_cov[as.character(g7),]
pred_cov_g8 <- pred_cov[as.character(g8),]
pred_cov_g9 <- pred_cov[as.character(g9),]
pred_cov_g10 <- pred_cov[as.character(g10),]

vartab <- rbind(apply(pred_cov_g1, 2, var), apply(pred_cov_g2, 2, var),
                apply(pred_cov_g3, 2, var), apply(pred_cov_g4, 2, var),
                apply(pred_cov_g5, 2, var), apply(pred_cov_g6, 2, var),
                apply(pred_cov_g7, 2, var), apply(pred_cov_g8, 2, var),
                apply(pred_cov_g9, 2, var), apply(pred_cov_g10, 2, var))
#vartab <- rbind(apply(pred_cov_g1, 2, mean), apply(pred_cov_g2, 2, mean), apply(pred_cov_g3, 2, mean), apply(pred_cov_g4, 2, mean), apply(pred_cov_g5, 2, mean))
df <- data.frame(cbind(Group = c(1:10), vartab))
df$Group <- as.factor(df$Group)
df <- melt(df, id.vars = "Group")
p <- ggplot(df, aes(x=variable, y=value, group=Group)) +
  geom_line(aes(color=Group))+
  geom_point(aes(color=Group)) +
  ylim(0, 1.5) +
  xlab("Predictive biomarkers") + ylab("Variance") +
  theme_classic() +
  theme(axis.text.x=element_text(angle=45, hjust=1)) + ggtitle("Treatment 2")
p
#ggsave(p, device = "pdf", path = "figs", filename = "var_arr_group.pdf")
table(labels)
vartab <- rbind(apply(pred_cov_g1, 2, mean), apply(pred_cov_g2, 2, mean),
                apply(pred_cov_g3, 2, mean), apply(pred_cov_g4, 2, mean),
                apply(pred_cov_g5, 2, mean), apply(pred_cov_g6, 2, mean),
                apply(pred_cov_g7, 2, mean), apply(pred_cov_g8, 2, mean),
                apply(pred_cov_g9, 2, mean), apply(pred_cov_g10, 2, mean))
df <- data.frame(cbind(Group = c(1:10), vartab))
df$Group <- as.factor(df$Group)
df <- melt(df, id.vars = "Group")
#df <- df[which(df$Group != 3),]
#df <- df[which(df$Group != 4),]
p <- ggplot(df, aes(x=variable, y=value, group=Group)) +
  geom_line(aes(color=Group))+
  geom_point(aes(color=Group)) +
  ylim(-2.0, 1.5) +
  xlab("Predictive biomarkers") + ylab("Mean") +
  theme_classic() +
  theme(axis.text.x=element_text(angle=45, hjust=1)) + ggtitle("Treatment 2")
p
