



### Real data analysis 
### Generate single data set sF2/sY/sD.

rm(list=ls()) 

#setwd("/users/sbaek/Google Drive/05 real/00 core data/")
setwd("/Users/hoyen/Desktop/Yanyuan/Code/SeungchulCode/real")
tmpd <- read.table('dataWname.dat')
bm <- read.table('gg05182017.dat') 
# 1-Basal(185) 2-Her2(109) 3-LumB(438) 4-LumA(246) 5-Normal(13)

# remove male cases
mind <- which(tmpd[,1]==0)
tmpdata <- tmpd[-mind,]
# 1-Basal(186) 2-Her2(108) 3-LumB(439) 4-LumA(248) 5-Normal(10)
# censoring rate = 880/980=89.8%
negy <- which(tmpdata[,5]<0) # nagative Y
tmpdata <- tmpdata[-negy,] # negative Y removed
# 1-Basal(185) 2-Her2(107) 3-LumB(428) 4-LumA(248) 5-Normal(10)
# censoring rate = 880/980=89.8%

# combine BRCA1 and BRCA2
for (i in 1:978){
  if (tmpdata[i,2]==1 | tmpdata[i,3]==1) {tmpdata[i,2]=1} else {tmpdata[i,2]=0}
}
# from now on, the 2nd column of tmpdata is COMBINED BRCA1 & BRCA2


save(tmpdata, file="full.RData")



setwd("/work/hoyen/Yanyuan/Data/Rerun2019")
source("/work/hoyen/Yanyuan/Code/createFolds.R")

bm <- read.table('gg05182017.dat') 
load("full.RData")

colMax<-function(x){
	ans<-apply(x, 2, max)
	return(ans)
}


gdata<-tmpdata[,8:ncol(tmpdata)]
gdata[1:5,1:5]
gmax<-colMax(gdata)
r<-which(gmax <= 5)
length(r)
gdata2<-gdata[,-r]
dim(gdata2)
#> dim(gdata2)
#[1]   978 18468
gdata3<-t(gdata2)


indtmp <- match(unlist(bm),colnames(tmpdata[,8:19156])) # 12 genes are NA...not found. 
Ai <- indtmp[-which(is.na(indtmp))] # is same for every column...132 genes in total.
### 3/4/2018
geneind <- c(1:19149)
Ai <- geneind[-Ai]
aa=colnames(tmpdata)[8:19156]

######################
### PENALIZED ADMM ###
######################
# In ADMM() function, rho & and lambda values should be tuned based on data. 
# rho & and lambda values should be tuned based on data. 
# (Let's compare means and std.dev. of F0 and original factors.)
# Ai: index set in the paper 

stage2 <- function(x,q,Ai){
  n <- dim(x)[2]
  p <- dim(x)[1]
  
  svdx <- eigen(t(x)%*%x)
  F0 <- sqrt(n)*t(svdx$vector[,1:q])
  B0 <- n^(-1)*x%*%t(F0)
  C0 <- B0
  Y0 <- matrix(0,p,q) #c(rep(1,p))%*%t(c(rep(1,q)))
  rho <- 30;lam <- 30  # 5000000 and 5000000
  tmpC <- matrix(0,p,q);tmpY <- matrix(0,p,q);tmpF<-matrix(0,q,n)
  tolr <- 10000;tols <- 10000;tol3 <- 10000;
  ABSTOL <- 1e-4;RELTOL <- 1e-2;
  eps_pri <- 0;eps_dual <- 0
  maxiter <- 0
#  errnum <- 0
  
   while (maxiter<200 & tols > eps_pri | tolr > eps_dual) {
    maxiter <- maxiter+1
#    if (inherits( try(solve(t(B0)%*%B0), silent=FALSE),"try-error")){errnum <- 1} 
#    else {tmpF <- solve(t(B0)%*%B0) %*% t(B0) %*% x}
    tmpF <- solve(t(B0)%*%B0) %*% t(B0) %*% x
    tmpB <- (2*n^(-1)*x%*%t(tmpF) - Y0 + rho*C0) %*% solve(2*n^(-1)*tmpF%*%t(tmpF)+rho*diag(q))
    for (i in 1:p){
      for (j in 1:q){
        if (is.na(match(i,Ai))) {tmpC[i,j]=tmpB[i,j]+rho^(-1)*Y0[i,j]} # ADDED because of penalzied.
        else { if (tmpB[i,j] > rho^(-1)*(-Y0[i,j]+lam)) {tmpC[i,j] <- tmpB[i,j]+rho^(-1)*(Y0[i,j]-lam)} 
        else { if (tmpB[i,j] < rho^(-1)*(-Y0[i,j]-lam)) {tmpC[i,j] <- tmpB[i,j]+rho^(-1)*(Y0[i,j]+lam)} 
          else {tmpC[i,j]=0} } }
      }
    }
    tmpY <- Y0 + rho*(tmpB-tmpC)
    
    # stopping criterion
    eps_pri <- ABSTOL+RELTOL*max(norm(tmpB,'1'),norm(tmpC,'1'))
    eps_dual <- ABSTOL+RELTOL*norm(rho*tmpY,'1')
    tols <- norm(rho*(tmpC-C0),'1')
    tolr <- norm(tmpB-tmpC,'1')
    tol3 <- norm(tmpF-F0,'1')
    
    B0 <- tmpB
    C0 <- tmpC 
    Y0 <- tmpY
    F0 <- tmpF
   }
  
  results=list(B0,F0)
  return(results)
}



###############################
### Cross-validation for q ###
##############################
cv.q <- NULL
kk <- 10
for (q in 1:40){
  # create k folds
  fold1 <- createFolds(1:(dim(gdata3)[2]), k=kk)  
  
  # prediction average squared error
  pe <- 0
  for (i in 1:kk){
  	print(paste("start fold:  ", i))
    # index for testing group
    te1.idx <- fold1[[i]] 
    # training
    tr1 <- gdata3[,-te1.idx]
    # testing
    te1 <- gdata3[,te1.idx] 
    
    s1r.tr <- stage2(tr1,q,Ai)
    B1.tr<- s1r.tr[[1]]
    
    # Calculte F1 in testing group from loading matrix B1
    F1.te <-solve(t(B1.tr)%*%B1.tr) %*% t(B1.tr) %*% te1
    te1.hat <- B1.tr%*%F1.te
    mean((te1-te1.hat)^2)
    pe <- pe+ mean((te1-te1.hat)^2)
  }
  pmse<-(pe/kk)/q
  cv.q <- cbind(cv.q, c(q, pe/kk, pmse))
}

save(cv.q, file="CVnumq.RData")


###############################
### plotting                ###
##############################

setwd("/Users/hoyen/Desktop/Yanyuan/Data/Rerun2019")
load("CVnumq.RData")
cv.q[,1:10]


pdf("/Users/hoyen/Desktop/Yanyuan/Graph/CVq.pdf", width=10, height=5)
par(mfrow=c(1,2))
plot(cv.q[1,], cv.q[2,], type="b", xlab="q", ylab="MSE")
abline(v=8)
plot(cv.q[1,], cv.q[3,], type="b", xlab="q", ylab="MSE/q")
abline(v=8)
dev.off()


dd1<-diff(cv.q[2,])
dd2<-diff(dd1)
plot(1:length(dd1), dd1)
plot(1:length(dd2), dd2)







