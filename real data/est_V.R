### May 2, 2019
### Estimation of V matrices for each subtype

rm(list=ls()) 

setwd("/users/baek/Google Drive/05 real/00 core data/")
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

indtmp <- match(unlist(bm),colnames(tmpdata[,8:19156])) # 12 genes are NA...not found. 
Ai <- indtmp[-which(is.na(indtmp))] # is same for every column...132 genes in total.
### 3/4/2018
geneind <- c(1:19149)
Ai <- geneind[-Ai]
aa=colnames(tmpdata)[8:19156]
aa[Ai]
write.table(aa[Ai],"132.dat", sep="\t",col.names = F, row.names = TRUE)
###

c1ind <- which(tmpdata[,7]==1)
c2ind <- which(tmpdata[,7]==2)
c3ind <- which(tmpdata[,7]==3)
c4ind <- which(tmpdata[,7]==4)
c5ind <- which(tmpdata[,7]==5)

s1 <- t(tmpdata[c1ind,8:19156]);w1 <- t(tmpdata[c1ind,c(2,4)]);Y1 <- tmpdata[c1ind,5];D1 <- tmpdata[c1ind,6]
s2 <- t(tmpdata[c2ind,8:19156]);w2 <- t(tmpdata[c2ind,c(2,4)]);Y2 <- tmpdata[c2ind,5];D2 <- tmpdata[c2ind,6]
s3 <- t(tmpdata[c3ind,8:19156]);w3 <- t(tmpdata[c3ind,c(2,4)]);Y3 <- tmpdata[c3ind,5];D3 <- tmpdata[c3ind,6]
s4 <- t(tmpdata[c4ind,8:19156]);w4 <- t(tmpdata[c4ind,c(2,4)]);Y4 <- tmpdata[c4ind,5];D4 <- tmpdata[c4ind,6]
s5 <- t(tmpdata[c5ind,8:19156]);w5 <- t(tmpdata[c5ind,c(2,4)]);Y5 <- tmpdata[c5ind,5];D5 <- tmpdata[c5ind,6]

### 5/2/2019
cm <- read.table('cm2.txt')

rname.s1 <- rownames(s1)
rname.cm <- rownames(cm)
tmpind <- NULL
for (i in 1:length(rname.cm)){
  place <- which(rname.s1==rname.cm[i])
  if (length(place) != 0) {tmpind[i] <- place}
  else {tmpind[i] <- 0} # 1099 out of 7840 don't exist among 19149 genes. 
}

# Modify CM matrix
empty <- which(tmpind==0)
new.cm <- cm[-empty,]

setwd("/Users/baek/Dropbox/99 PhD/Research/07 pro3/11 revision2/")
write.table(new.cm,"newcm2.txt") # newcm with row/column names
save(new.cm, file = "newcm2.RData")

# C'C invertible?
C <- matrix(unlist(new.cm),6741,330)
det(t(C) %*% C) # not invertible....

# Rearrange X(s1-s5) matrices




# save X based on category
#write.table(s1,"cat1.dat", sep="\t",col.names = F, row.names = TRUE)
#write.table(s2,"cat2.dat", sep="\t",col.names = F, row.names = TRUE)
#write.table(s3,"cat3.dat", sep="\t",col.names = F, row.names = TRUE)
#write.table(s4,"cat4.dat", sep="\t",col.names = F, row.names = TRUE)
#write.table(s5,"cat5.dat", sep="\t",col.names = F, row.names = TRUE)

########################
##### Main program #####
########################
library(MASS)

# set up some values
q <- 8

ptm <- proc.time()

s1r <- stage2(s1,q,Ai)
s2r <- stage2(s2,q,Ai)
s3r <- stage2(s3,q,Ai)
s4r <- stage2(s4,q,Ai)
s5r <- stage2(s5,q,Ai)

estF <- cbind(s1r[[2]],s2r[[2]],s3r[[2]],s4r[[2]],s5r[[2]])

dataF <- estF 
dataD <- c(D1,D2,D3,D4,D5)
dataY <- c(Y1,Y2,Y3,Y4,Y5)
dataW <- cbind(w1,w2,w3,w4,w5)
dataFs <- rbind(dataF,dataW)

proc.time() - ptm 

write.table(dataY,"sY.dat", sep="\t",col.names = F, row.names = F)
write.table(dataD,"sD.dat", sep="\t",col.names = F, row.names = F)
write.table(dataFs,"sF.dat", sep="\t",col.names = F, row.names = F)

# save factor loadings...19149 by 8
write.table(s1r[[1]],"B1.dat", sep="\t",col.names = F, row.names = TRUE)
write.table(s2r[[1]],"B2.dat", sep="\t",col.names = F, row.names = TRUE)
write.table(s3r[[1]],"B3.dat", sep="\t",col.names = F, row.names = TRUE)
write.table(s4r[[1]],"B4.dat", sep="\t",col.names = F, row.names = TRUE)
write.table(s5r[[1]],"B5.dat", sep="\t",col.names = F, row.names = TRUE)
#########
###### Remove BRCA1 & BRCA2
#########

br <- dataFs[seq(9,nrow(dataFs),10),] # BRCA1 and BRCA2

newestf <- dataFs[-seq(9,nrow(dataFs),10),]
write.table(newestf,"sF2.dat", sep="\t",col.names = F, row.names = F) # gender.BRCA1/2 removed for test....

##############################
##### Defining Functions #####
##############################

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

###################
### NORMAL ADMM ###
###################
# In ADMM() function, rho & and lambda values should be tuned based on data. 
# rho & and lambda values should be tuned based on data. 
# (Let's compare means and std.dev. of F0 and original factors.)

stage1 <- function(x,q){
  n <- dim(x)[2]
  p <- dim(x)[1]
  
  svdx <- eigen(t(x)%*%x)
  F0 <- sqrt(n)*t(svdx$vector[,1:q]) 
  B0 <- n^(-1)*x%*%t(F0)
  C0 <- B0
  Y0 <- matrix(0,p,q) #c(rep(1,p))%*%t(c(rep(1,q)))
  rho <- 30;lam <- 30 
  tmpC <- matrix(0,p,q);tmpY <- matrix(0,p,q);
  tolr <- 10000;tols <- 10000;tol3 <- 10000;
  ABSTOL <- 1e-4;RELTOL <- 1e-2;
  eps_pri <- 0;eps_dual <- 0
  maxiter <- 0
  
  while (maxiter<200 & tols > eps_pri | tolr > eps_dual) {
    maxiter <- maxiter+1
    tmpF <- solve(t(B0)%*%B0) %*% t(B0) %*% x
    tmpB <- (2*n^(-1)*x%*%t(tmpF) - Y0 + rho*C0) %*% solve(2*n^(-1)*tmpF%*%t(tmpF)+rho*diag(q))
    for (i in 1:p){
      for (j in 1:q){
        if (tmpB[i,j] > rho^(-1)*(-Y0[i,j]+lam)) {tmpC[i,j] <- tmpB[i,j]+rho^(-1)*(Y0[i,j]-lam)} 
        else { if (tmpB[i,j] < rho^(-1)*(-Y0[i,j]-lam)) {tmpC[i,j] <- tmpB[i,j]+rho^(-1)*(Y0[i,j]+lam)} 
          else {tmpC[i,j]=0} }
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


#########
###### Remove BRCA1 & BRCA2
#########

rtd=read.table('Y.dat')
rty=read.table('D.dat')
rtf=read.table('F.dat')

rtf[1:20,1:5]

br <- rtf[seq(9,nrow(rtf),10),] # BRCA1 and BRCA2

newestf <- rtf[-seq(9,nrow(rtf),10),]
write.table(newestf,"F2.dat", sep="\t",col.names = F, row.names = F) # gender.BRCA1/2 removed for test....