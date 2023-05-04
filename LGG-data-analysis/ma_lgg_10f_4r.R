### Calculate the Training-Test results
rm(list=ls())
set.seed(121)

library(ConsensusClusterPlus)
library("mvtnorm")
library("glmnetcr")
library(parallel)
library(doParallel)


load("data/LGGdata.rda")
#matchRTComp <- matchRTComp[sample(1:nrow(matchRTComp), size = nrow(matchRTComp), replace = F),]
#mypath <- c(paste0("output/simulation-scenarios/train-test/scen-alt-", sc))

predcov <- scale(matchRTComp[,16:38])  ## Predictive features
progcov <- scale(matchRTComp[,c(11,13)])  ## Prognostic features

trtsgn <- as.numeric(paste(matchRTComp[,10]))
yord <- as.numeric(paste(matchRTComp[,9]))
ymat <- matrix(0, nrow = length(yord), ncol = max(yord)+1)
for(i in 1:(length(yord))){
  ymat[i, (yord[i] + 1)] <- 1
}

data <- list()
data$pred <- predcov
data$prog <- progcov
data$trtsgn <- trtsgn
data$yord <- yord
data$ymat <- ymat

################################ Functions ########################################
mymultt <- function(Xtrain, X.pred){
  myln <- length(Xtrain[,1])
  myls <- Xtrain
  mylmu <- apply(myls, 2, mean)
  mymun <- myln/(myln+kappa0)*mylmu
  myS <- cov(myls)*(myln-1)

  ## for the covariates
  myd <- length(Xtrain[1,])       ## numer of covariates
  nu0 <- length(Xtrain[1,])+1       ## numer of covariates +1;
  lambda0 <- diag(myd)                ## identity matrix

  kappan <- kappa0 + myln
  nun <- nu0 + myln
  lambdn <- lambda0 + myS + kappa0*myln/(kappa0 + myln)*(mylmu - mu0)%*%t(mylmu - mu0)
  return2 <- dmvt(x = X.pred, sigma=(kappan + 1)*lambdn/(kappan*(nun - myd + 1)), df = nun - myd + 1, log = FALSE)
}

con.cluster <- function(cons, ymat = s_train_ymat,  yvec = s_train_yord,
                        prog = s_train_prog, trt = s_train_trt){
  ni <- dim(cons)[1]
  utpred1 <- matrix(1, nrow = ni, ncol=12)  ## ut1,ut2,trt

  for (i in 1:ni){
    totut <- ymat*cons[i,]  ### utility
    totpre <- totut[-i,]
    trtpre <- trt[-i]

    l_train_prog <- prog[-i,]
    l_yvec <- yvec[-i];

    mytempt <- cbind(totpre, trtpre, l_yvec, l_train_prog)

    ## for genes
    trt1g <- subset(mytempt[,1:3], mytempt[,4] == 0)
    trt2g <- subset(mytempt[,1:3], mytempt[,4] == 1)
    alphatrt1 <- apply(trt1g, 2, sum) + prior1
    probpre1_temp <- alphatrt1/sum(alphatrt1)
    alphatrt2 <- apply(trt2g, 2, sum) + prior2
    probpre2_temp <- alphatrt2/sum(alphatrt2)

    mycovX <- prog[i,];                  ## observed covariates

    trtc <- subset(mytempt[,-(1:4)])

    trtc0 <- subset(trtc,trtc[,1] == 0)
    trtc0dst <- mymultt(trtc0[,-1], mycovX)
    trtc1 <- subset(trtc,trtc[,1] == 1)
    trtc1dst <- mymultt(trtc1[,-1], mycovX)
    trtc2 <- subset(trtc, trtc[,1] == 2)
    trtc2dst <- mymultt(trtc2[,-1], mycovX)
    trtpfy <- c(trtc0dst, trtc1dst, trtc2dst)   ##prob with cov

    probpre1 <- probpre1_temp*trtpfy/sum(probpre1_temp*trtpfy)
    probpre2 <- probpre2_temp*trtpfy/sum(probpre2_temp*trtpfy)

    ## calculate the utility
    ut1pre <- probpre1%*%wk
    ut2pre <- probpre2%*%wk
    utpred1[i, 1] <- ut1pre
    utpred1[i, 2] <- ut2pre
    if(ut2pre > ut1pre){
      utpred1[i, 3] = 2
    }
    utpred1[i, 4:6] <- probpre1
    utpred1[i, 7:9] <- probpre2
    utpred1[i, 10:12] <- trtpfy/(sum(trtpfy))
    ut1pre<-ut2pre<-totut<-totpre<-mytempt<-trt1g<-trt2g<-alphatrt1<-alphatrt2<-NULL;
  }
  return <- utpred1;
}

PreUt<-function(mth, trt, out.response, SUB.ID){
  #mth<-mth;
  myresults <- cbind(mth, trt+1, out.response, SUB.ID)
  pred1 <- subset(myresults, myresults[,3] == 1)
  table1 <- table(pred1[,14], pred1[,13])
  pred2 <- subset(myresults, myresults[,3] == 2)
  table2 <- table(pred2[,14], pred2[,13])
  p1 <- sum(table1)/(sum(table1) + sum(table2))
  p2 <- sum(table2)/(sum(table1) + sum(table2))
  ## set prob as 0 for cases that none selected patients had response (all 0s);
  ## set prob as 1 for cases that none selected patients had non-response (all 1s);
  if(length(table1) == 4){
    crt1<-table1[2,1]/sum(table1[,1])
  }
  if(length(table1) < 4){
    crt1 <- 0
  }
  #### if table 1 is empty then set crt1=0;
  if(length(table2) == 4){
    crt2 <- table2[2,2]/sum(table2[,2])
  }
  if(length(table2) < 4){
    crt2 <- 0
  }

  #### summary meaures
  return <- crt1*p1 + crt2*p2 - sum(out.response)/length(out.response)
}

################################ Functions ########################################
mymultt <- function(Xtrain, X.pred){
  myln <- tryCatch(expr = length(Xtrain[,1]), error = function(e){return(1)})
  myls <- Xtrain
  mylmu <- tryCatch(expr = apply(myls, 2, mean), error = function(e){return(myls)})
  mymun <- myln/(myln+kappa0)*mylmu
  myS <- tryCatch(expr = cov(myls), error = function(e){return(0)})*(myln-1)

  ## for the covariates
  myd <- tryCatch(expr = length(Xtrain[1,]), error = function(e){return(length(Xtrain))})       ## numer of covariates
  nu0 <- myd+1       ## numer of covariates +1;
  lambda0 <- diag(myd)                ## identity matrix

  kappan <- kappa0 + myln
  nun <- nu0 + myln
  lambdn <- tryCatch(expr = lambda0 + myS + kappa0*myln/(kappa0 + myln)*(mylmu - mu0)%*%t(mylmu - mu0), error = function(e){return(diag(1, 2, 2))})
  return2 <- dmvt(x = X.pred, sigma=(kappan + 1)*lambdn/(kappan*(nun - myd + 1)), df = nun - myd + 1, log = FALSE)
}

# calculate the NPC
countUT <- function(resultsum, myoutot){
  myctut <- matrix(0, nrow = 3, ncol = 3)
  myctutSum <- NULL
    mycurdata <- resultsum
    mypre <- NULL
    pretrt1 <- apply(mycurdata[,4:6], 1, which.max)
    pretrt2 <- apply(mycurdata[,7:9], 1, which.max)
    mypreTall <- cbind(pretrt1, pretrt2)
    for(j in 1:(dim(resultsum)[1])){
      mypre[j] <- mypreTall[j, trtsgn[j]+1]
    }
    sts <- table(mypre, factor(myoutot, levels = c(0:2)))
    mysdls <- as.numeric(rownames(sts))
    str1 <- matrix(0, nrow = 3, ncol = 3)
    str1[mysdls,] <- sts

    myctut <- str1*diag(3)
    myctutSum <- sum(str1*diag(3))

  return <- cbind(myctutSum)
}

################################ setup Parameters ########################################

wk <- c(0,40,100)
prior1 <- c(1/3,1/3,1/3)
prior2 <- c(1/3,1/3,1/3)
kappa0 <- 1
mu0 <- c(0,0)
K <- 10

registerDoParallel(cores = K)

vectf <- c(1, 17, 33, 49, 65, 81, 97, 113, 129, 145, 159)
HC.sum.all <- vector(mode = "list", length = K)
for(k in 1:K){
    currfold <- (vectf[k]:(vectf[k+1]-1))
    train_pred <- data$pred[-currfold,]
    train_prog <- data$prog[-currfold,]
    train_yord <- data$yord[-currfold]#-1
    train_ymat <- data$ymat[-currfold,]
    train_trt <- data$trtsgn[-currfold]#-1

    n <- dim(train_pred)[1]

    HC.sum<-matrix(0,nrow=n,ncol=9)
    #HC.sum<-matrix(0,nrow=n,ncol=4)

    for (mysub in 1:n){
      s_train_pred <- train_pred[-mysub,]
      s_train_prog <- train_prog[-mysub,]
      s_train_yord <- train_yord[-mysub]
      s_train_ymat <- train_ymat[-mysub,]
      s_train_trt <- train_trt[-mysub]

      ### clustering using CONSENSUS MATRIX method ###################################
      con_clu <- ConsensusClusterPlus(t(s_train_pred),maxK=15,reps=500,pItem=0.90,pFeature=1,
                                      #clusterAlg="hc",distance="pearson",
                                      #clusterAlg="km",distance="euclidean",
                                      clusterAlg="pam",distance="manhattan",
                                      seed=126)

      hc2 <- con.cluster(con_clu[[2]][["consensusMatrix"]], ymat = s_train_ymat,  yvec = s_train_yord,
                         prog = s_train_prog, trt = s_train_trt)
      hc3 <- con.cluster(con_clu[[3]][["consensusMatrix"]], ymat = s_train_ymat,  yvec = s_train_yord,
                         prog = s_train_prog, trt = s_train_trt)
      hc4 <- con.cluster(con_clu[[4]][["consensusMatrix"]], ymat = s_train_ymat,  yvec = s_train_yord,
                         prog = s_train_prog, trt = s_train_trt)
      hc5 <- con.cluster(con_clu[[5]][["consensusMatrix"]], ymat = s_train_ymat,  yvec = s_train_yord,
                         prog = s_train_prog, trt = s_train_trt)
      hc6 <- con.cluster(con_clu[[6]][["consensusMatrix"]], ymat = s_train_ymat,  yvec = s_train_yord,
                         prog = s_train_prog, trt = s_train_trt)
      hc7 <- con.cluster(con_clu[[7]][["consensusMatrix"]], ymat = s_train_ymat,  yvec = s_train_yord,
                         prog = s_train_prog, trt = s_train_trt)
      hc8 <- con.cluster(con_clu[[8]][["consensusMatrix"]], ymat = s_train_ymat,  yvec = s_train_yord,
                         prog = s_train_prog, trt = s_train_trt)
      hc9 <- con.cluster(con_clu[[9]][["consensusMatrix"]], ymat = s_train_ymat,  yvec = s_train_yord,
                         prog = s_train_prog, trt = s_train_trt)
      hc10 <- con.cluster(con_clu[[10]][["consensusMatrix"]], ymat = s_train_ymat,  yvec = s_train_yord,
                          prog = s_train_prog, trt = s_train_trt)

      out.response <- as.numeric(s_train_yord > 1)
      SUB.ID <- c(1:(n-1))

      HC.sum[mysub,] <- c(PreUt(hc2, s_train_trt, out.response, SUB.ID),
                          PreUt(hc3, s_train_trt, out.response, SUB.ID),
                          PreUt(hc4, s_train_trt, out.response, SUB.ID),
                          PreUt(hc5, s_train_trt, out.response, SUB.ID),
                          PreUt(hc6, s_train_trt, out.response, SUB.ID),
                          PreUt(hc7, s_train_trt, out.response, SUB.ID),
                          PreUt(hc8, s_train_trt, out.response, SUB.ID),
                          PreUt(hc9, s_train_trt, out.response, SUB.ID),
                          PreUt(hc10, s_train_trt, out.response, SUB.ID))

    }
    HC.sum.all[[k]] <- HC.sum
  }

###capisci da qui!
max_clus_vec <- c()
#max_clus_vec <- rep(5, K)

for(k in 1:K){
  temp <- HC.sum.all[[k]]
  max_clus_vec[k] <- median(apply(temp,1,which.max)+1)
}
#max_clus_vec <- rep(5, K)
utpred1APT.all <- matrix(0, nrow = 158, ncol = 19)

### clustering using CONSENSUS MATRIX method ###################################

rst.hc<-ConsensusClusterPlus(t(data$pred),maxK=15,reps=500,pItem=0.90,pFeature=1,
                             #clusterAlg="hc",distance="pearson",
                             #clusterAlg="km",distance="euclidean",
                             clusterAlg="pam",distance="manhattan",
                             seed=126);


for(myrep in 1:K){
  currfold <- (vectf[myrep]:(vectf[myrep+1]-1))
  trtAPT <- data$trtsgn[currfold]#-1
  Rapp <- data$prog[currfold,]
  outcomAPT <- data$yord[currfold]#-1
  trtsgn <- data$trtsgn[currfold]
  myoutot <- data$ymat[currfold,]

  ntemp <- length(trtAPT)

  utpred1APT<-matrix(1,nrow= ntemp,ncol=19)  ### ut1,ut2,trt,cluster

  ### pick the median rank with the largest summary measure
  #max.clus <-apply(HC.sum.all[,,myrep],1,which.max)+1
  max.clus <- max_clus_vec[myrep]

  for (mysub in 1:ntemp){
    trt <- trtAPT[-mysub]
    outcom <- outcomAPT[-mysub]
    select.sub.n <- length(outcom)
    myRapp <- Rapp[-mysub,]

    mycovXSub <- Rapp[mysub,]                 ## observed covariates
    trtcSub <- cbind(outcom, myRapp)

    trtc0Sub <- subset(trtcSub, trtcSub[,1] == 0)
    trtc0dstSub <- mymultt(trtc0Sub[,-1], mycovXSub)
    trtc1Sub <- subset(trtcSub, trtcSub[,1] == 1)
    trtc1dstSub <- mymultt(trtc1Sub[,-1], mycovXSub)
    trtc2Sub <- subset(trtcSub, trtcSub[,1] == 2)
    trtc2dstSub <- mymultt(trtc2Sub[,-1], mycovXSub)
    trtpfySub <- c(trtc0dstSub, trtc1dstSub, trtc2dstSub)   ##prob with cov

    myyAPT <- matrix(0, nrow = ntemp, ncol = 3)
    for(m in 1:ntemp){
      myyAPT[m, outcomAPT[m]+1] = 1
    }

    mycons1APT <- rst.hc[[max.clus]][["consensusMatrix"]]

    totutAPT <- myyAPT*mycons1APT[vectf[myrep] + mysub - 1, c(currfold)]               ### utility
    totpreAPT <- totutAPT[-mysub,]
    trtpreAPT <- trtAPT[-mysub]
    mytemptAPT <- cbind(totpreAPT, trtpreAPT)

    ## calculate alpha hat
    trt1gAPT <- subset(mytemptAPT[,1:3], mytemptAPT[,4] == 0)
    trt2gAPT <- subset(mytemptAPT[,1:3], mytemptAPT[,4] == 1)
    alphatrt1APT <- apply(trt1gAPT, 2, sum) + prior1
    probpre1APT <- alphatrt1APT/sum(alphatrt1APT)

    alphatrt2APT <- apply(trt2gAPT, 2, sum) + prior2
    probpre2APT <- alphatrt2APT/sum(alphatrt2APT)

    probpre1 <- probpre1APT*trtpfySub/sum(probpre1APT*trtpfySub)
    probpre2 <- probpre2APT*trtpfySub/sum(probpre2APT*trtpfySub)

    ## calculate the utility
    ut1preAPT <- probpre1%*%wk
    ut2preAPT <- probpre2%*%wk
    utpred1APT[mysub, 1] <- ut1preAPT
    utpred1APT[mysub, 2] <- ut2preAPT
    if(ut2preAPT > ut1preAPT){
      utpred1APT[mysub, 3] = 2
    }
    utpred1APT[mysub, 4:6] <- probpre1
    utpred1APT[mysub, 7:9] <- probpre2
    utpred1APT[mysub, 10:12] <- probpre1APT
    utpred1APT[mysub, 13:15] <- probpre2APT
    utpred1APT[mysub, 16:18] <- trtpfySub/sum(trtpfySub)
    utpred1APT[mysub, 19] <- max.clus#[mysub]
  }

  utpred1APT.all[currfold,] <- utpred1APT

}

#####################Save the results#######################################


case2HCppUT <- utpred1APT.all
MOT <- MTUg <- NPC <- c()
for(k in 1:K){
  wk<-c(0,40,100)
  currfold <- (vectf[k]:(vectf[k+1]-1))
  outcomAPT <- data$yord[currfold]
  HCppCUT <- as.vector(countUT(case2HCppUT[currfold,], outcomAPT))
  NPC[k] <- HCppCUT
}


NPC
sum(NPC)

mth<-utpred1APT.all
trtsgn <- as.numeric(paste(matchRTComp[,10]))
myoutot <- as.numeric(matchRTComp[,9])
myresults<-cbind(mth,as.numeric(trtsgn),as.numeric(myoutot>2))
pred1<-subset( myresults,myresults[,3]==1);
table1<-table(pred1[,21],pred1[,20]);
pred2<-subset( myresults,myresults[,3]==2);
table2<-table(pred2[,21],pred2[,20]);
p1<-sum(table1)/(sum(table1)+sum(table2));
p2<-sum(table2)/(sum(table1)+sum(table2));
########################################################################################
## set prob as 0 for cases that none selected patients had response (all 0s);
## set prob as 1 for cases that none selected patients had non-response (all 1s);
########################################################################################
if (length(table1)==4){ crt1<-table1[2,1]/sum(table1[,1])};
if (length(table1)<4) { crt1<-as.numeric(row.names(table1))}
if(length(crt1) == 0) {crt1<-0}

if (length(table2)==4){  crt2<-table2[2,2]/sum(table2[,2])};
if (length(table2)<4) {  crt2<-as.numeric(row.names(table2))};
### summary meaures
return<-crt1*p1+crt2*p2-sum(as.numeric(myoutot>2))/158
return
