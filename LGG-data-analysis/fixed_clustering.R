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
library(grid)

gg_color_hue <- function(n) {
  hues = seq(15, 375, length = n + 1)
  hcl(h = hues, l = 65, c = 100)[1:n]
}
#n = 7
#cols = gg_color_hue(n)
#dev.new(width = n, height = n)
#plot(1:n, pch = 16, cex = 2, col = cols)
cols1 <- gg_color_hue(3)
cols2 <- gg_color_hue(4)

load("data/LGGdata.rda")
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

K <- 1#numero di fold

wk <- c(0, 40, 100)

registerDoParallel(cores = 50)#alloco solo core necessari

Y <- matrix(0, nrow = npat, ncol = max(as.numeric(matchRTComp[,9])))
for(i in 1:nrow(Y)){
  Y[i, as.numeric(matchRTComp[i,9])] <- 1
}

table(matchRTComp[,9:10])
vectf <- c(1, 17, 33, 49, 65, 81, 97, 113, 129, 145, 159)

nout <- 2000
X <- scale(matchRTComp[,16:38])
Z <- scale(matchRTComp[,c(11,13)])

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

#res0 <- ppmxct(y = data.matrix(Y), X = data.frame(X),
#                               Xpred = data.frame(X[1:2,]), Z = data.frame(Z),
#                               Zpred = data.frame(Z[1:2,]), asstreat = trtsgn, #treatment,
#                               PPMx = 1, cohesion = 2, kappa = 1,
#                               similarity = 2, consim = 2, similparam = vec_par,
#                               calibration = 2, coardegree = 2, modelpriors = modelpriors,
#                               update_hierarchy = T,
#                               hsp = T, iter = iterations, burn = burnin, thin = thinning,
#                               mhtunepar = c(0.05, 0.05), CC = 5, reuse = 1,
#                               nclu_init = 3)
#
#save(res0, file = "output/lgg-da/lgg_rep4clunks.RData")
load("output/lgg-da/lgg_rep4clunks.RData")

cls2 <- t(as.matrix(res0$label[[1]]))#[,c(1:num_treat[1])]
psm2 <- comp.psm(cls2)
ccm <- psm2

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
  theme(legend.position="bottom") +#,
  #legend.title = element_text("Co-occurence")) +
  labs(fill = "Co-occurrence") +
  xlab("Patients") + ylab("Patients") + ggtitle("Treatment 1")#+ ggtitle("Heatmap of averaged co-occurence matrix")

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
df <- data.frame(cbind(Group = c(1:10), vartab))
df$Group <- as.factor(df$Group)
df <- reshape2::melt(df, id.vars = "Group")
df <- df[which(df$Group != 3),]
df <- df[which(df$Group != 5),]

var1 <- ggplot(df, aes(x=variable, y=value, group=Group)) +
  geom_line(aes(color=Group))+
  geom_point(aes(color=Group)) +
  ylim(0, 1.5) +
  xlab("Predictive biomarkers") + ylab("Variance") +
  theme_classic() +
  theme(axis.text.x=element_text(angle=45, hjust=1)) + ggtitle("Treatment 1")

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
df <- df[which(df$Group != 3),]
df <- df[which(df$Group != 5),]
df <- df[which(df$Group != 6),]
df <- df[which(df$Group != 7),]
df <- df[which(df$Group != 8),]
df <- df[which(df$Group != 9),]
df <- df[which(df$Group != 10),]
mean1 <- ggplot(df, aes(x=variable, y=value, group=Group)) +
  geom_line(aes(color=Group))+
  geom_point(aes(color=Group)) +
  scale_colour_manual(values=c(cols1)) +
  ylim(-1.0, 1.0) +
  xlab("Predictive biomarkers") + ylab("Mean") +
  theme_classic() +
  geom_hline(yintercept=0) +
  theme(axis.text.x=element_text(angle=45, hjust=1)) + ggtitle("Treatment 1")

cls2 <- t(as.matrix(res0$label[[2]]))#[,c(1:num_treat[1])]
psm2 <- comp.psm(cls2)
ccm <- psm2

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
c2 <- ggplot(mC, aes(Var1,Var2, fill=`Co-occurrence`)) + geom_raster() +
  scale_fill_continuous(type = "viridis") +
  theme_classic() +
  theme(legend.position="bottom") +
  xlab("Patients") + ylab("Patients") + ggtitle("Treatment 2")#+ ggtitle("Heatmap of averaged co-occurence matrix")

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
df <- df[which(df$Group != 9),]
df <- reshape2::melt(df, id.vars = "Group")
var2 <- ggplot(df, aes(x=variable, y=value, group=Group)) +
  geom_line(aes(color=Group))+
  geom_point(aes(color=Group)) +
  ylim(0, 1.5) +
  xlab("Predictive biomarkers") + ylab("Variance") +
  theme_classic() +
  theme(axis.text.x=element_text(angle=45, hjust=1)) + ggtitle("Treatment 2")

vartab <- rbind(apply(pred_cov_g1, 2, mean), apply(pred_cov_g2, 2, mean),
                apply(pred_cov_g3, 2, mean), apply(pred_cov_g4, 2, mean),
                apply(pred_cov_g5, 2, mean), apply(pred_cov_g6, 2, mean),
                apply(pred_cov_g7, 2, mean), apply(pred_cov_g8, 2, mean),
                apply(pred_cov_g9, 2, mean), apply(pred_cov_g10, 2, mean))
df <- data.frame(cbind(Group = c(1:10), vartab))
df$Group <- as.factor(df$Group)
df <- reshape2::melt(df, id.vars = "Group")
df <- df[which(df$Group != 4),]
df <- df[which(df$Group != 5),]
df <- df[which(df$Group != 7),]
df <- df[which(df$Group != 8),]
df <- df[which(df$Group != 9),]
df <- df[which(df$Group != 10),]

mean2 <- ggplot(df, aes(x=variable, y=value, group=Group)) +
  geom_line(aes(color=Group))+
  geom_point(aes(color=Group)) +
  scale_colour_manual(values=c(cols2)) +
  ylim(-1.0, 1.0) +
  geom_hline(yintercept=0) +
  xlab("Predictive biomarkers") + ylab("Mean") +
  theme_classic() +
  theme(axis.text.x=element_text(angle=45, hjust=1)) + ggtitle("Treatment 2")


grid_arrange_shared_legend <-
  function(...,
           ncol = length(list(...)),
           nrow = 1,
           position = c("bottom", "right")) {

    plots <- list(...)
    position <- match.arg(position)
    g <-
      ggplotGrob(plots[[1]] + theme(legend.position = position))$grobs
    legend <- g[[which(sapply(g, function(x)
      x$name) == "guide-box")]]
    lheight <- sum(legend$height)
    lwidth <- sum(legend$width)
    gl <- lapply(plots, function(x)
      x + theme(legend.position = "none"))
    gl <- c(gl, ncol = ncol, nrow = nrow)

    combined <- switch(
      position,
      "bottom" = arrangeGrob(
        do.call(arrangeGrob, gl),
        legend,
        ncol = 1,
        heights = unit.c(unit(1, "npc") - lheight, lheight)
      ),
      "right" = arrangeGrob(
        do.call(arrangeGrob, gl),
        legend,
        ncol = 2,
        widths = unit.c(unit(1, "npc") - lwidth, lwidth)
      )
    )

    grid.newpage()
    grid.draw(combined)

    # return gtable invisibly
    invisible(combined)

  }

coocc_plot <- grid_arrange_shared_legend(c1, c2)#, nrow = 1, top = "Title of the page")
#ggsave(coocc_plot, device = "pdf", path = "figs", filename = "coocc_plot.pdf")
var_arr_group <- grid.arrange(var1, var2, nrow = 1)#, top = "Title of the page")
#ggsave(var_arr_group, device = "pdf", path = "figs", filename = "var_arr_group.pdf")
mean_arr_group <- grid.arrange(mean1, mean2, nrow = 1)#, top = "Title of the page")
ggsave(mean_arr_group, device = "pdf", path = "figs", filename = "mean_arr_group.pdf")

#fix_clu <- cbind(labels1, labels2)

iterations <- 12000
burnin <- 2000
thinning <- 5
SP <- (iterations-burnin)/thinning
modelpriors <- list()
modelpriors$hP0_m0 <- rep(0, ncol(Y)); modelpriors$hP0_nu0 <- 10
modelpriors$hP0_s0 <- ncol(Y) + 2; modelpriors$hP0_Lambda0 <- 10

vec_par <- c(0.0, 1.0, .5, 1.0, 2.0, 2.0, 0.1)
#
#fclu <- ppmxct_fixed(y = data.matrix(Y), X = data.frame(X),
#                     Z = data.frame(Z), asstreat = trtsgn, #treatment,
#                     PPMx = 1, cohesion = 1, kappa = c(1, 20, 10, 1), sigma = c(0.01, .59, 6),
#                     similarity = 2, consim = 2, similparam = vec_par, calibration = 2,
#                     coardegree = 2, modelpriors, update_hierarchy = T, hsp = T,
#                     iter = iterations, burn = burnin, thin = thinning,
#                     mhtunepar = c(0.05, 0.05), CC = 5, reuse = 1, fix_clu = fix_clu)
#
#save(fclu, file = "output/lgg-da/lggfixedclunks.RData")
load("output/lgg-da/lggfixedclunks.RData")

pi1 <- fclu$pi_out[which(trtsgn == 1),,]
pi11 <- pi1[which(labels1 == 1),,]

pi11 <- t(apply(pi11, c(2, 3), mean))

colnames(pi11) <- c("PD", "PS", "CR")
p11 <- ggtern(data=data.frame(pi11), aes(x=PD, y=PS, z=CR)) + #geom_point()
  geom_density_tern(colour = cols1[1], bdl = 0.01) +
  ggtitle("Group 1")

pi12 <- pi1[which(labels1 == 2),,]

pi12 <- t(apply(pi12, c(2, 3), mean))
colnames(pi12) <- c("PD", "PS", "CR")
p12 <- ggtern(data=data.frame(pi12), aes(x=PD, y=PS, z=CR)) + #geom_point() +
  geom_density_tern(colour = cols1[2], bdl = 0.1) +
  ggtitle("Group 2")

pi14 <- pi1[which(labels1 == 4),,]
pi14 <- abs(pi14 + rnorm((dim(pi14)[1]), 0, .05))

pi14 <- as.data.frame(t(apply(pi14, c(2, 3), mean)))
#pi14 <- pi14[sample(1:2000, 100),]
colnames(pi14) <- c("PD", "PS", "CR")
p14 <- ggtern(data=pi14, aes(x=PD, y=PS, z=CR)) + #geom_point() +
  geom_density_tern(colour = cols1[3], bdl = 0.1) +
  ggtitle("Group 4")

tern1 <- ggtern::grid.arrange(p11, p12, p14, ncol = 3)

pi2 <- fclu$pi_out[which(trtsgn == 2),,]
pi21 <- pi2[which(labels2 == 1),,]

pi21 <- t(apply(pi21, c(2, 3), mean))+.01
pi21 <- pi21/apply(pi21,1, sum)

colnames(pi21) <- c("PD", "PS", "CR")
p21 <- ggtern(data=data.frame(pi21), aes(x=PD, y=PS, z=CR)) + #geom_point()
  geom_density_tern(colour = cols2[1], bdl = 0.05) +
  ggtitle("Group 1")

pi22 <- pi2[which(labels2 == 2),,]

pi22 <- t(apply(pi22, c(2, 3), mean))

colnames(pi22) <- c("PD", "PS", "CR")
p22 <- ggtern(data=data.frame(pi22), aes(x=PD, y=PS, z=CR)) +
  geom_density_tern(colour = cols2[2], bdl = 0.1) +
  ggtitle("Group 2")

pi23 <- pi2[which(labels2 == 3),,]

pi23 <- as.data.frame(t(apply(pi23, c(2, 3), mean)))
colnames(pi23) <- c("PD", "PS", "CR")
p23 <- ggtern(data=pi23, aes(x=PD, y=PS, z=CR)) +
  geom_density_tern(colour = cols2[3], bdl = 0.01) +
  ggtitle("Group 3")

pi26 <- pi2[which(labels2 == 6),,]

pi26 <- as.data.frame(t(apply(pi26, c(2, 3), mean)))
colnames(pi26) <- c("PD", "PS", "CR")
p26 <- ggtern(data=pi26, aes(x=PD, y=PS, z=CR)) +
  geom_density_tern(colour = cols2[4], bdl = 0.1) +
  ggtitle("Group 6")

tern2 <- ggtern::grid.arrange(p21, p22, p23, p26, ncol = 4)

tp <- ggtern::grid.arrange(tern1, tern2, nrow = 2)

#ggsave(tern1, device = "pdf", path = "figs", filename = "tern1.pdf")
#ggsave(tern2, device = "pdf", path = "figs", filename = "tern2.pdf")

#ggsave(tp, device = "pdf", path = "figs", filename = "tp.pdf")


