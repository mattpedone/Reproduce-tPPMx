rm(list=ls())
set.seed(121)

library(rstan)
library(parallel)
library(doParallel)
library(doRNG)

registerDoParallel(cores = 10)
rstan_options(auto_write = TRUE)

load("data/LGGdata.rda")
matchRTComp <- matchRTComp[sample(1:nrow(matchRTComp), size = nrow(matchRTComp), replace = F),]
trtsgn <- c(matchRTComp[,10]) + 1
npat <- length(trtsgn)

K <- 10#numero di fold

Y <- matrix(0, nrow = npat, ncol = max(as.numeric(matchRTComp[,9])))
for(i in 1:nrow(Y)){
  Y[i, as.numeric(matchRTComp[i,9])] <- 1
}

table(matchRTComp[,9:10])
vectf <- c(1, 17, 33, 49, 65, 81, 97, 113, 129, 145, 159)
X <- scale(matchRTComp[,16:38])
Z <- scale(matchRTComp[,c(11,13)])

npat_pred <- 28

n <- npat-npat_pred#dim(Yl[[1]])[1]
J <- ncol(Y)
p <- ncol(X)
q <- ncol(Z)

wk <- c(0, 25, 75, 100)

my_level <- c("Progressive Disease", "Stable Disease", "Partial Remission/Response", "Complete Remission/Response")
tempy <- factor(matchRTComp[,3], levels = my_level)
tempy <- as.integer(tempy)
Y <- matrix(0, nrow = npat, ncol = max(tempy))
for(i in 1:nrow(Y)){
  Y[i, tempy[i]] <- 1
}


res <- foreach(k = 1:K) %dorng%
  {
    currfold <- (vectf[k]:(vectf[k+1]-1))

    Y_train <- Y[-currfold,]
    ntrain <- dim(Y_train)[1]
    Y_test <- rbind(Y[currfold,], Y[currfold,])
    Y <- rbind(Y_train, Y_test)
    X_train <- X[-currfold,]
    X_test <- rbind(X[currfold,], X[currfold,])
    X <- rbind(X_train, X_test)
    Z_train <- Z[-currfold,]
    Z_test <- rbind(Z[currfold,], Z[currfold,])
    Z <- rbind(Z_train, Z_test)
    trt_train <- trtsgn[-currfold]-1#trtl[[r]][1:124]-1
    trt_test <- c(rep(0, length(currfold)), rep(1, length(currfold)))
    trt <- c(trt_train, trt_test)
    Xdis <- cbind(X, Z, trt, (trt*X))
    ntot <- dim(Xdis)[1]

    X_train <- Xdis[1:ntrain, ]
    Y_train <- Y[1:ntrain, ]

    X_test <- Xdis[(ntrain+1):ntot,]
    Y_test <- Y[(ntrain+1):ntot,]

    ss_data = list(N = dim(X_train)[1], M = dim(X_test)[1], J = dim(Y_train)[2],
                    P = dim(X_train)[2], X = X_train, Xp = X_test, Y = Y_train,
                    sd_prior = 1.0, psi = .25)

    fit <- rstan::stan(file = "R/dmhs-scripts/model.stan", data = ss_data, cores = 1, iter = 1000,
                        chains = 1, verbose = T, warmup = 200, seed = 121,
                        control = list(max_treedepth = 15, adapt_delta = 0.995))#995))


    params <- rstan::extract(fit)
    pp <- apply(params$pipred, c(2,3), median)

    out <- cbind(pp1 = pp[1:length(currfold),], pp2 = pp[(length(currfold)+1):(2*length(currfold)),])
    return(out)
  }

#NPC
npc_tf <- function(output, trtsgn, myoutot){
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
myt <- list()
for(k in 1:K){
  currfold <- (vectf[k]:(vectf[k+1]-1))
  myoutot <- as.numeric(matchRTComp[currfold,9])#simdata$yord[[k]][131:158,]
  trtsgn_test <- trtsgn[currfold]#simdata$trtsgn[[k]][131:158]

  pp1 <- res[[k]][,1:3]
  A1 <- pp1%*%wk
  pp2 <- res[[k]][,4:6]
  A2 <- pp2%*%wk

  myt[[k]] <- as.numeric(A1 < A2) + 1

  temp[currfold,] <- cbind(pp1, pp2)
}
NPC <- npc_tf(temp, trtsgn, as.numeric(matchRTComp[,9]))

#results
myoutot <- as.numeric(matchRTComp[,9])#simdata$yord[[k]][131:158,]
mytab <- cbind(myass = unlist(myt), rndass = trtsgn, resp = as.numeric(myoutot>2))
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

NPC; ESM
