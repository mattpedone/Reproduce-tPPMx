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

load("data/LGGdata.rda")
matchRTComp <- matchRTComp[sample(1:nrow(matchRTComp), size = nrow(matchRTComp), replace = F),]
trtsgn <- c(matchRTComp[,10]) + 1
#trtsgn <- ifelse(trtsgn == 1, 2, 1)
npat <- length(trtsgn)
#trtsgn <- sample(1:2, npat, replace = TRUE)

K <- 1#numero di fold

predAPT_all <- matrix(0, nrow = npat, ncol = 9)
nclust_all <- matrix(0, nrow = K, ncol = 6)
gof_all <- matrix(0, nrow = K, ncol = 2)

wk <- c(0, 40, 100)

registerDoParallel(cores = 50)#alloco solo core necessari

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
                                   PPMx = 1, cohesion = 2, kappa = c(1, 20, 10, 1), sigma = c(0.01, .59, 6),
                                   similarity = 2, consim = 2, similparam = vec_par,
                                   calibration = 2, coardegree = 2, modelpriors,
                                   update_hierarchy = T,
                                   hsp = T, iter = iterations, burn = burnin, thin = thinning,
                                   mhtunepar = c(0.05, 0.05), CC = 5, reuse = 1,
                                   nclu_init = 3), error = function(e){FALSE})
    return(res0)
  }

save(myres0, file = "output/lgg-da/lgg_rep4clunks.RData")
#load("output/lgg-da/lgg_rep4clunks.RData")

ccm <- matrix(0, nrow=79, ncol = 79)
for(k in 1:K){
  res0 <- myres0[[k]]
  #number of a cluster, mean, binder &varinf ----
  cls2 <- t(as.matrix(res0$label[[1]]))#[,c(1:num_treat[1])]
  psm2 <- comp.psm(cls2)
  ccm <- ccm + psm2
}

ccmf <- as.matrix(ccm/K)

labels1 <- minVI(ccmf)$cl
table(labels1)

orp1 <- as.numeric(rownames(X[which(trtsgn == 1),]))
g1 <- orp1[which(labels1==1)]
g2 <- orp1[which(labels1==2)]
g3 <- orp1[which(labels1==3)]
g4 <- orp1[which(labels1==4)]
g5 <- orp1[which(labels1==5)]
g6 <- orp1[which(labels1==6)]
g7 <- orp1[which(labels1==7)]
g8 <- orp1[which(labels1==8)]
g9 <- orp1[which(labels1==9)]
g10 <- orp1[which(labels1==10)]

reord <- c()
for(i in 1:max(labels1)){
  reord <- c(reord, which(labels1 == i))
}

# Co-occurence plot ----
data1 <- ccmf[,reord]
data <- data1[reord,]
colnames(data) <- rownames(data) <- NULL
#coincidences<-sapply(1:ncol(data), function(i){ colSums(data[,i]==data) })
mC <- reshape2::melt(data)
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
df <- reshape2::melt(df, id.vars = "Group")
df <- df[which(df$Group != 3),]
df <- df[which(df$Group != 5),]

p <- ggplot(df, aes(x=variable, y=value, group=Group)) +
  geom_line(aes(color=Group))+
  geom_point(aes(color=Group)) +
  ylim(0, 1.5) +
  xlab("Predictive biomarkers") + ylab("Variance") +
  theme_classic() +
  theme(axis.text.x=element_text(angle=45, hjust=1)) + ggtitle("Treatment 2")
p
#ggsave(p, device = "pdf", path = "figs", filename = "var_arr_group.pdf")
vartab <- rbind(apply(pred_cov_g1, 2, mean), apply(pred_cov_g2, 2, mean),
                apply(pred_cov_g3, 2, mean), apply(pred_cov_g4, 2, mean),
                apply(pred_cov_g5, 2, mean), apply(pred_cov_g6, 2, mean),
                apply(pred_cov_g7, 2, mean), apply(pred_cov_g8, 2, mean),
                apply(pred_cov_g9, 2, mean), apply(pred_cov_g10, 2, mean))
#apply(pred_cov_g4, 2, mean), #apply(pred_cov_g8, 2, mean),
#apply(pred_cov_g9, 2, mean))#, apply(pred_cov_g10, 2, mean))
df <- data.frame(cbind(Group = c(1:10), vartab))
df$Group <- as.factor(df$Group)
df <- reshape2::melt(df, id.vars = "Group")
df <- df[which(df$Group != 2),]
df <- df[which(df$Group != 5),]
df <- df[which(df$Group != 6),]
df <- df[which(df$Group != 7),]
df <- df[which(df$Group != 8),]
df <- df[which(df$Group != 10),]
p <- ggplot(df, aes(x=variable, y=value, group=Group)) +
  geom_line(aes(color=Group))+
  geom_point(aes(color=Group)) +
  ylim(-1.0, 1.0) +
  xlab("Predictive biomarkers") + ylab("Mean") +
  theme_classic() +
  theme(axis.text.x=element_text(angle=45, hjust=1)) + ggtitle("Treatment 2")
p
table(matchRTComp[which(labels1 == 1),9])#higher the better response
table(matchRTComp[which(labels1 == 3),9])#higher the better response
table(matchRTComp[which(labels1 == 4),9])#higher the better response
table(matchRTComp[which(labels1 == 9),9])#higher the better response

ccm <- matrix(0, nrow=79, ncol = 79)
for(k in 1:K){
  res0 <- myres0[[k]]
  #number of a cluster, mean, binder &varinf ----
  cls2 <- t(as.matrix(res0$label[[2]]))#[,c(1:num_treat[1])]
  psm2 <- comp.psm(cls2)
  ccm <- ccm + psm2
}

ccmf <- as.matrix(ccm/K)

labels2 <- minVI(ccmf)$cl
table(labels2)

orp2 <- as.numeric(rownames(X[which(trtsgn == 2),]))
g1 <- orp2[which(labels2==1)]
g2 <- orp2[which(labels2==2)]
g3 <- orp2[which(labels2==3)]
g4 <- orp2[which(labels2==4)]
g5 <- orp2[which(labels2==5)]
g6 <- orp2[which(labels2==6)]
g7 <- orp2[which(labels2==7)]
g8 <- orp2[which(labels2==8)]
g9 <- orp2[which(labels2==9)]
g10 <- orp2[which(labels2==10)]

reord <- c()
for(i in 1:max(labels2)){
  reord <- c(reord, which(labels2 == i))
}

# Co-occurence plot ----
data1 <- ccmf[,reord]
data <- data1[reord,]
colnames(data) <- rownames(data) <- NULL
#coincidences<-sapply(1:ncol(data), function(i){ colSums(data[,i]==data) })
mC <- reshape2::melt(data)
colnames(mC) <- c("Var1", "Var2", "Co-occurrence")
c1 <- ggplot(mC, aes(Var1,Var2, fill=`Co-occurrence`)) + geom_raster() +
  scale_fill_continuous(type = "viridis") +
  theme_classic() +
  theme(legend.position="bottom") +
  xlab("Patients") + ylab("Patients") + ggtitle("Treatment 2")#+ ggtitle("Heatmap of averaged co-occurence matrix")
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
                apply(pred_cov_g7, 2, var), apply(pred_cov_g8,  2, var),
                apply(pred_cov_g9, 2, var), apply(pred_cov_g10, 2, var))
#vartab <- rbind(apply(pred_cov_g1, 2, mean), apply(pred_cov_g2, 2, mean), apply(pred_cov_g3, 2, mean), apply(pred_cov_g4, 2, mean), apply(pred_cov_g5, 2, mean))
df <- data.frame(cbind(Group = c(1:10), vartab))
df$Group <- as.factor(df$Group)
#df <- df[which(df$Group != 9),]
df <- reshape2::melt(df, id.vars = "Group")
p <- ggplot(df, aes(x=variable, y=value, group=Group)) +
  geom_line(aes(color=Group))+
  geom_point(aes(color=Group)) +
  ylim(0, 1.5) +
  xlab("Predictive biomarkers") + ylab("Variance") +
  theme_classic() +
  theme(axis.text.x=element_text(angle=45, hjust=1)) + ggtitle("Treatment 2")
p
#ggsave(p, device = "pdf", path = "figs", filename = "var_arr_group.pdf")
table(labels2)
vartab <- rbind(apply(pred_cov_g1, 2, mean), apply(pred_cov_g2, 2, mean),
                apply(pred_cov_g3, 2, mean), apply(pred_cov_g4, 2, mean),
                apply(pred_cov_g5, 2, mean), apply(pred_cov_g6, 2, mean),
                apply(pred_cov_g7, 2, mean), apply(pred_cov_g8, 2, mean),
                apply(pred_cov_g9, 2, mean), apply(pred_cov_g10, 2, mean))
df <- data.frame(cbind(Group = c(1:10), vartab))
df$Group <- as.factor(df$Group)
df <- reshape2::melt(df, id.vars = "Group")
df <- df[which(df$Group != 3),]
df <- df[which(df$Group != 5),]
df <- df[which(df$Group != 6),]
df <- df[which(df$Group != 7),]
df <- df[which(df$Group != 8),]
df <- df[which(df$Group != 9),]
df <- df[which(df$Group != 10),]

p <- ggplot(df, aes(x=variable, y=value, group=Group)) +
  geom_line(aes(color=Group))+
  geom_point(aes(color=Group)) +
  ylim(-2.0, 2.0) +
  xlab("Predictive biomarkers") + ylab("Mean") +
  theme_classic() +
  theme(axis.text.x=element_text(angle=45, hjust=1)) + ggtitle("Treatment 2")
p

table(matchRTComp[which(labels2 == 1),9])#higher the better response
table(matchRTComp[which(labels2 == 2),9])#higher the better response
table(matchRTComp[which(labels2 == 3),9])#higher the better response
table(matchRTComp[which(labels2 == 4),9])#higher the better response
table(matchRTComp[which(labels2 == 9),9])#higher the better response

fix_clu <- cbind(labels1, labels2)
#table(labels1)
#table(labels2)

iterations <- 12000
burnin <- 2000
thinning <- 5
SP <- (iterations-burnin)/thinning
modelpriors <- list()
modelpriors$hP0_m0 <- rep(0, ncol(Y)); modelpriors$hP0_nu0 <- 10
modelpriors$hP0_s0 <- ncol(Y) + 2; modelpriors$hP0_Lambda0 <- 10

vec_par <- c(0.0, 1.0, .5, 1.0, 2.0, 2.0, 0.1)


fclu <- ppmxct_fixed(y = data.matrix(Y), X = data.frame(X),
                     Z = data.frame(Z), asstreat = trtsgn, #treatment,
                     PPMx = 1, cohesion = 1, kappa = c(1, 20, 10, 1), sigma = c(0.01, .59, 6),
                     similarity = 2, consim = 2, similparam = vec_par, calibration = 2,
                     coardegree = 2, modelpriors, update_hierarchy = F, hsp = T,
                     iter = iterations, burn = burnin, thin = thinning,
                     mhtunepar = c(0.05, 0.05), CC = 5, reuse = 1, fix_clu = fix_clu)

save(fclu, file = "output/lgg-da/lggfixedclunks.RData")
load("output/lgg-da/lggfixedclunks.RData")

pi1 <- fclu$pi_out[which(trtsgn == 1),,]
pi11 <- pi1[which(labels1 == 1),,]

pi11 <- t(apply(pi11, c(2, 3), mean))

colnames(pi11) <- c("PD", "PS", "CR")
p11 <- ggtern(data=data.frame(pi11), aes(x=PD, y=PS, z=CR)) + #geom_point()
  geom_density_tern(colour = "#F8766D", bdl = 0.1) +
  ggtitle("Group 1")

pi13 <- pi1[which(labels1 == 3),,]

pi13 <- t(apply(pi13, c(2, 3), mean))
#pi12 <- abs(pi12 + rnorm((dim(pi12)[1]), 0, .01))
#pi12 <- pi12/apply(pi12, 1, sum)

colnames(pi13) <- c("PD", "PS", "CR")
p13 <- ggtern(data=data.frame(pi13), aes(x=PD, y=PS, z=CR)) + #geom_point() +
  geom_density_tern(colour = "#F8766D", bdl = 0.1) +
  ggtitle("Group 3")
p13

pi14 <- pi1[which(labels1 == 4),,]
#pi14 <- abs(pi14 + rnorm((dim(pi14)[1]), 0, .05))

pi14 <- as.data.frame(t(apply(pi14, c(2, 3), mean)))
#pi14 <- pi14[sample(1:2000, 100),]
colnames(pi14) <- c("PD", "PS", "CR")
p14 <- ggtern(data=pi14, aes(x=PD, y=PS, z=CR)) + #geom_point() +
  geom_density_tern(colour = "#F8766D", bdl = 0.1) +
  ggtitle("Group 4")


pi19 <- pi1[which(labels1 == 9),,]

pi19 <- t(apply(pi19, c(2, 3), mean))

colnames(pi19) <- c("PD", "PS", "CR")
p19 <- ggtern(data=data.frame(pi19), aes(x=PD, y=PS, z=CR)) +
  geom_density_tern(colour = "#F8766D", bdl = 0.1) +
  ggtitle("Group 9")

tern1 <- ggtern::grid.arrange(p11, p13, p14, p19, ncol = 2)

pi2 <- fclu$pi_out[which(trtsgn == 2),,]
pi21 <- pi2[which(labels2 == 1),,]

pi21 <- t(apply(pi21, c(2, 3), mean))+.01
#pi21 <- pi21+abs(rnorm(dim(pi21)[1], 0, 0.01))
pi21 <- pi21/apply(pi21,1, sum)

colnames(pi21) <- c("PD", "PS", "CR")
p21 <- ggtern(data=data.frame(pi21), aes(x=PD, y=PS, z=CR)) + #geom_point()
  geom_density_tern(colour = "#F8766D", bdl = 0.1) +
  ggtitle("Group 1")
p21

pi22 <- pi2[which(labels2 == 2),,]

pi22 <- t(apply(pi22, c(2, 3), mean))

colnames(pi22) <- c("PD", "PS", "CR")
p22 <- ggtern(data=data.frame(pi22), aes(x=PD, y=PS, z=CR)) +
  geom_density_tern(colour = "#F8766D", bdl = 0.1) +
  ggtitle("Group 2")

pi23 <- pi2[which(labels2 == 3),,]

pi23 <- as.data.frame(t(apply(pi23, c(2, 3), mean)))
colnames(pi23) <- c("PD", "PS", "CR")
p23 <- ggtern(data=pi23, aes(x=PD, y=PS, z=CR)) +
  geom_density_tern(colour = "#F8766D", bdl = 0.1) +
  ggtitle("Group 3")

pi24 <- pi2[which(labels2 == 4),,]

pi24 <- as.data.frame(t(apply(pi24, c(2, 3), mean)))
colnames(pi24) <- c("PD", "PS", "CR")
p24 <- ggtern(data=pi24, aes(x=PD, y=PS, z=CR)) +
  geom_density_tern(colour = "#F8766D", bdl = 0.1) +
  ggtitle("Group 4")

pi29 <- pi2[which(labels2 == 9),,]
pi29 <- as.data.frame(t(apply(pi29, c(2, 3), mean)))
colnames(pi29) <- c("PD", "PS", "CR")
p29 <- ggtern(data=pi29, aes(x=PD, y=PS, z=CR)) +
  geom_density_tern(colour = "#F8766D", bdl = 0.1) +
  ggtitle("Group 9")

tern2 <- ggtern::grid.arrange(p21, p22, p23, p24, p29, ncol = 3)

