library(dplyr)
library(tibble)
library(tidyselect)
library(ggplot2)
library(vweights)
library(treatppmx)
set.seed(121)

nobs <- 50

gendata <- function(n = 1000, pro = c(0.2,0.5,0.3), dim = 3){
  if(length(pro) != 3){
    stop("The length of probabilities must be equal to 3!")
  }
  data <- matrix(nrow = n, ncol = dim)
  for(i in 1:n){
    u <- runif(1)
    if(u < pro[1]){
      #cat(i,"ciao1","\n")
      data[i,] <- mvtnorm::rmvnorm(1, mean = rep(-2.1, dim), sigma = diag(.5, dim))
    }
    else{
      if(u < (pro[1] + pro[2])){
        #cat(i,"ciao2","\n")
        data[i,] <- mvtnorm::rmvnorm(1, mean = rep(0, dim), sigma = diag(.5, dim))
      }
      else{
        #cat(i,"ciao3","\n")
        data[i,] <- mvtnorm::rmvnorm(1, mean = rep(2.3, dim), sigma = diag(.5, dim))
      }
    }
  }
  return (data)
}

X <- gendata(n = nobs, dim = 5)

vec_par <- c(0.0, 2.0, .5, 1.0, 2.0, 2.0, 0.1)
#double m0=0.0, s20=10.0, v=.5, k0=1.0, nu0=2.0, n0 = 2.0;
iterations <- 10000; burnin <- 0; thinning <- 1

#beta <- 48.4185; sigma <- .25; theta <- 19.233
beta <- 1.0; sigma <- .7553; theta <- 19.233

nout <- (iterations-burnin)/thinning

####---- Number of Clusters ----#####

ngg_nocal <- prior_ppmx(X = X, PPMx = 1, cohesion = 2, alpha = beta, #*sigma,
                             sigma = sigma, similarity = 2, consim = 2,
                             similparam = vec_par, calibration = 0,
                             coardegree = 1, iter = iterations, burn = burnin,
                             thin = thinning, nclu_init = 30)
res1 <- ngg_nocal$nclu

dp_coa <- prior_ppmx(X = X, PPMx = 1, cohesion = 1, alpha = theta, #1, #2,
                             sigma = sigma, similarity = 2, consim = 2,
                             similparam = vec_par, calibration = 2,
                             coardegree = 2, iter = iterations, burn = burnin,
                             thin = thinning, nclu_init = 30)
res2 <- dp_coa$nclu

ngg_coa <- prior_ppmx(X = X, PPMx = 1, cohesion = 2, alpha = beta, #*sigma,
                             sigma = sigma, similarity = 2, consim = 2,
                             similparam = vec_par, calibration = 2,
                             coardegree = 2, iter = iterations, burn = burnin,
                             thin = thinning, nclu_init = 30)
res3 <- ngg_coa$nclu


dp <- prior_ppmx(X = X, PPMx = 0, cohesion = 1, alpha = theta,
                             sigma = sigma, similarity = 1, consim = 1,
                             similparam = vec_par, calibration = 0,
                             coardegree = 1, iter = iterations, burn = burnin,
                             thin = thinning, nclu_init = 30)
res4 <- dp$nclu

res <- t(rbind(res1/nout, res2/nout, res3/nout, res4/nout, 
               c(round(vweights::computepnclu(nobs, sigma, beta*sigma), 5))))

res <- (res[which(rowSums(res) != 0),])
res[which(res==0)] <- NA
res <- rownames_to_column(as.data.frame(res), var = "cluster")
newres <- cbind(similarity = c(rep("NGGP-nocal", nrow(res)), 
                               #rep("NGG-dd-cal", nrow(res)), 
                               rep("DP-sim", nrow(res)), 
                               rep("NGGP-sim", nrow(res)), 
                               rep("DP", nrow(res)), rep("NGGP", nrow(res))),
                cluster = rep(res[,1], 5),
                freq = c(unlist(res[, c(2:6)])))

dfres <- newres %>%
  as_tibble() %>% # lo trasformo in un dataframe
  #tibble::rownames_to_column(var = "cluster") %>%
  mutate(cluster = as.numeric(cluster)) %>%
  mutate(freq = as.numeric(freq)) #%>%
  #mutate(sd = as.numeric(sd))

plot1 <- ggplot(dfres, aes(x=cluster, y=freq, group=similarity, color=similarity)) +
  #geom_errorbar(aes(ymin=freq-sd, ymax=freq+sd), width=.1) +
  geom_line() +
  geom_point()+
  scale_color_manual(values=c("#E69F00", "#56B4E9", "#009E73", "#0072B2", "#999999"))+
  theme_classic() +
  #xlab(expression(c)) +
  #ylab(expression(P(C[n] == c)))
  labs(x = expression(c), y = expression(P(C[50] == c)), 
       color = expression(' '))

plot1 + theme(legend.position=c(.75,.75), legend.key.size = unit(1, 'cm'))

ggsave("figs/plot_pnc.pdf")

##tentativi grafico su scala di grigi invece che colorato
#ggplot(dfres, aes(x=cluster, y=freq, group=similarity, color=similarity, linetype=similarity)) +
#  geom_line() +
#  scale_colour_grey() +
#  theme_minimal() +
#  scale_linetype_manual(values=c(1, 3, 2, 4, 5))
#  #geom_errorbar(aes(ymin=freq-sd, ymax=freq+sd), width=.1) +
#  geom_line(aes(linetype=similarity)) +
#  #geom_point()+
#  #scale_color_brewer(palette="Paired")+
#  #scale_linetype() +
#  #scale_color_manual(c("darkgrey", "darkgrey", "darkgrey", "darkgrey", "darkgrey")) +
#  #xlab(expression(c)) +
#  #ylab(expression(P(C[n] == c)))
#  labs(x = expression(c), y = expression(P(C[n] == c)), 
#       color = expression(' '))

####---- Cardinalities ----#####
tab <- apply(ngg_nocal$nj, 2, sort, decreasing = T)
tab[which(tab == 0)] <- NA
tab2 <- tab
tab2[which(tab2 != 0)] <- 1
avg_ngg_nocal <- mean(apply(tab2, 2, sum, na.rm = T))
tab[which(tab != 1)] <- 0
ps_ngg_nocal <- mean(apply(tab, 2, mean, na.rm = T))

tab <- apply(dp_coa$nj, 2, sort, decreasing = T)
tab[which(tab == 0)] <- NA
tab2 <- tab
tab2[which(tab2 != 0)] <- 1
avg_dp_coa <- mean(apply(tab2, 2, sum, na.rm = T))
tab[which(tab != 1)] <- 0
ps_dp_coa <- mean(apply(tab, 2, mean, na.rm = T))

tab <- apply(ngg_coa$nj, 2, sort, decreasing = T)
tab[which(tab == 0)] <- NA
tab2 <- tab
tab2[which(tab2 != 0)] <- 1
avg_ngg_coa <- mean(apply(tab2, 2, sum, na.rm = T))
tab[which(tab != 1)] <- 0
ps_ngg_coa <- mean(apply(tab, 2, mean, na.rm = T))

tab <- apply(dp$nj, 2, sort, decreasing = T)
tab[which(tab == 0)] <- NA
tab2 <- tab
tab2[which(tab2 != 0)] <- 1
avg_dp <- mean(apply(tab2, 2, sum, na.rm = T))
tab[which(tab != 1)] <- 0
ps_dp <- mean(apply(tab, 2, mean, na.rm = T))

ngg <- prior_ppmx(X = X, PPMx = 0, cohesion = 2, alpha = beta, #*sigma,
                  sigma = sigma, similarity = 2, consim = 2,
                  similparam = vec_par, calibration = 2,
                  coardegree = 2, iter = iterations, burn = burnin,
                  thin = thinning, nclu_init = 30)

tab <- apply(ngg$nj, 2, sort, decreasing = T)
tab[which(tab == 0)] <- NA
cs_ngg <- cumsum(vweights::computepnclu(nobs, sigma, beta))
avg_ngg <- 26+(1-(cs_ngg[27]-.5)/(cs_ngg[27]-cs_ngg[26]))
tab[which(tab != 1)] <- 0
ps_ngg <- mean(apply(tab, 2, mean, na.rm = T))

avgs <- c(avg_dp, avg_dp_coa, avg_ngg, avg_ngg_nocal, avg_ngg_coa)
pss <- c(ps_dp, ps_dp_coa, ps_ngg, ps_ngg_nocal, ps_ngg_coa)
tab <- rbind(avgs, pss)
colnames(tab) <- c("DP", "DP-sim", "NGG", "NGG-nocal", "NGG-sim")
rownames(tab) <- c("Average number of clusters", "Proportion of singletons")
tab
