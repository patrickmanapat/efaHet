# SETUP ##################################################################

#packages
library(MASS) #for mvrnorm
library(gdata) #for unmatrix
library(reshape2) #for melt

#constants
minrep <- 2001
maxrep <- 3000
mEta <- 0

# TESTING ################################################################

# #test runs
# n <- 200
# nEta <- 3
# pEta <- 12
# rEta <- .7
# l <- .45

# SIMULATION #############################################################

#begin simulation function
EFAsim <- function(n,nEta,pEta,rEta,l){

for(d in minrep:maxrep){
  
# SET FACTOR CORRELATION MATRIX ##########################################
  
S <- matrix(rEta, nEta, nEta)
  diag(S) <- 1
  
# SET FACTOR MEAN VECTOR #################################################
  
M <- rep(mEta, 
         times = nEta)

# GENERATE FACTOR SCORES #################################################

etaOG <- mvrnorm(n,
                 mu = M,
                 Sigma = S,
                 empirical = FALSE)
  
etaLong <- etaOG[rep(1:nrow(etaOG),times = pEta),]
  
etaMelt <- melt(etaLong)

eta <- etaMelt$value
  
# GENERATE OBSERVED SCORES ###############################################

#factor loadings
lambdaOG <- matrix(l, 
                   nrow = pEta, 
                   ncol = nEta)

lambdaMelt <- melt(lambdaOG)

lambda <- rep(lambdaMelt$value, 
              each = n)

#residual variances
dev <- rnorm(n*nEta*pEta, 0, 1)

u <- sqrt(1-lambda^2)*dev

#observed scores
y <- lambda * eta + u

yForm <- data.frame(lambda, eta, u, y)

yDF <- matrix(y,
              nrow = n)

# ########################################################################
# #TEST: DATA GENERATION#
# library(lavaan)
# 
# data <- as.data.frame(yDF)
# 
# mod1 <- '
#         eta1 =~ NA*V1 + V2 + V3 + V4 + V5 + V6 + V7 + V8 + V9 + V10 + V12
#         eta2 =~ NA*V13 + V14 + V15 + V16 + V17 + V18 + V19 + V20 + V21 + V22 + V23 + V24
#         eta3 =~ NA*V25 + V26 + V27 + V28 + V29 + V30 + V31 + V32 + V33 + V34 + V35 + V36
#         eta1 ~~ 1*eta1
#         eta2 ~~ 1*eta2
#         eta3 ~~ 1*eta3
#         '
# fit1 <- cfa(mod1,
#             data = data)
# 
# summary(fit1,
#         fit.measures = TRUE)
# ########################################################################
  
# SAVE EVERYTHING ########################################################

#write out inter-factor correlation matrix for every replication (d)
write.table(S,
            file = paste('./s3K/','nEta',nEta,'pEta',pEta,'rEta',rEta,'l',l,'n',n,'d',d,'.dat',sep=''),
            row.names = FALSE,
            col.names = FALSE,
            sep = '\t')

#write out factor scores for every replication (d)
write.table(etaOG,
            file = paste('./eta3K/','nEta',nEta,'pEta',pEta,'rEta',rEta,'l',l,'n',n,'d',d,'.dat',sep=''),
            row.names = FALSE,
            col.names = FALSE,
            sep = '\t')

#write out factor loadings for every replication (d)
write.table(lambdaOG,
            file = paste('./lambda3K/','nEta',nEta,'pEta',pEta,'rEta',rEta,'l',l,'n',n,'d',d,'.dat',sep=''),
            row.names = FALSE,
            col.names = FALSE,
            sep = '\t')

#write out computation of y for every replication (d)
write.table(yForm,
            file = paste('./yForm3K/','nEta',nEta,'pEta',pEta,'rEta',rEta,'l',l,'n',n,'d',d,'.dat',sep=''),
            row.names = FALSE,
            col.names = FALSE,
            sep = '\t')

#write out observed scores for every replication (d)
write.table(yDF,
            file = paste('./y3K/','nEta',nEta,'pEta',pEta,'rEta',rEta,'l',l,'n',n,'d',d,'.dat',sep=''),
            row.names = FALSE,
            col.names = FALSE,
            sep = '\t')

#CLOSE LOOPS##############################################################

}

#end simulation function
}

#SIMULATION CONDITIONS####################################################

set.seed(30143378)

#simulation runs: EFAsim(n,nEta,pEta,rEta,l)
EFAsim( 60,1, 8,0.3,0.45)
EFAsim(100,1, 8,0.3,0.45)
EFAsim(200,1, 8,0.3,0.45)
EFAsim(400,1, 8,0.3,0.45)

EFAsim( 60,3, 8,0.3,0.45)
EFAsim(100,3, 8,0.3,0.45)
EFAsim(200,3, 8,0.3,0.45)
EFAsim(400,3, 8,0.3,0.45)
EFAsim( 60,5, 8,0.3,0.45)
EFAsim(100,5, 8,0.3,0.45)
EFAsim(200,5, 8,0.3,0.45)
EFAsim(400,5, 8,0.3,0.45)

EFAsim( 60,1,12,0.3,0.45)
EFAsim(100,1,12,0.3,0.45)
EFAsim(200,1,12,0.3,0.45)
EFAsim(400,1,12,0.3,0.45)

EFAsim( 60,3,12,0.3,0.45)
EFAsim(100,3,12,0.3,0.45)
EFAsim(200,3,12,0.3,0.45)
EFAsim(400,3,12,0.3,0.45)
EFAsim( 60,5,12,0.3,0.45)
EFAsim(100,5,12,0.3,0.45)
EFAsim(200,5,12,0.3,0.45)
EFAsim(400,5,12,0.3,0.45)

EFAsim( 60,3, 8,0.5,0.45)
EFAsim(100,3, 8,0.5,0.45)
EFAsim(200,3, 8,0.5,0.45)
EFAsim(400,3, 8,0.5,0.45)
EFAsim( 60,5, 8,0.5,0.45)
EFAsim(100,5, 8,0.5,0.45)
EFAsim(200,5, 8,0.5,0.45)
EFAsim(400,5, 8,0.5,0.45)

EFAsim( 60,3,12,0.5,0.45)
EFAsim(100,3,12,0.5,0.45)
EFAsim(200,3,12,0.5,0.45)
EFAsim(400,3,12,0.5,0.45)
EFAsim( 60,5,12,0.5,0.45)
EFAsim(100,5,12,0.5,0.45)
EFAsim(200,5,12,0.5,0.45)
EFAsim(400,5,12,0.5,0.45)

EFAsim( 60,3, 8,0.7,0.45)
EFAsim(100,3, 8,0.7,0.45)
EFAsim(200,3, 8,0.7,0.45)
EFAsim(400,3, 8,0.7,0.45)
EFAsim( 60,5, 8,0.7,0.45)
EFAsim(100,5, 8,0.7,0.45)
EFAsim(200,5, 8,0.7,0.45)
EFAsim(400,5, 8,0.7,0.45)

EFAsim( 60,3,12,0.7,0.45)
EFAsim(100,3,12,0.7,0.45)
EFAsim(200,3,12,0.7,0.45)
EFAsim(400,3,12,0.7,0.45)
EFAsim( 60,5,12,0.7,0.45)
EFAsim(100,5,12,0.7,0.45)
EFAsim(200,5,12,0.7,0.45)
EFAsim(400,5,12,0.7,0.45)

EFAsim( 60,1, 8,0.3,0.71)
EFAsim(100,1, 8,0.3,0.71)
EFAsim(200,1, 8,0.3,0.71)
EFAsim(400,1, 8,0.3,0.71)

EFAsim( 60,3, 8,0.3,0.71)
EFAsim(100,3, 8,0.3,0.71)
EFAsim(200,3, 8,0.3,0.71)
EFAsim(400,3, 8,0.3,0.71)
EFAsim( 60,5, 8,0.3,0.71)
EFAsim(100,5, 8,0.3,0.71)
EFAsim(200,5, 8,0.3,0.71)
EFAsim(400,5, 8,0.3,0.71)

EFAsim( 60,1,12,0.3,0.71)
EFAsim(100,1,12,0.3,0.71)
EFAsim(200,1,12,0.3,0.71)
EFAsim(400,1,12,0.3,0.71)

EFAsim( 60,3,12,0.3,0.71)
EFAsim(100,3,12,0.3,0.71)
EFAsim(200,3,12,0.3,0.71)
EFAsim(400,3,12,0.3,0.71)
EFAsim( 60,5,12,0.3,0.71)
EFAsim(100,5,12,0.3,0.71)
EFAsim(200,5,12,0.3,0.71)
EFAsim(400,5,12,0.3,0.71)

EFAsim( 60,3, 8,0.5,0.71)
EFAsim(100,3, 8,0.5,0.71)
EFAsim(200,3, 8,0.5,0.71)
EFAsim(400,3, 8,0.5,0.71)
EFAsim( 60,5, 8,0.5,0.71)
EFAsim(100,5, 8,0.5,0.71)
EFAsim(200,5, 8,0.5,0.71)
EFAsim(400,5, 8,0.5,0.71)

EFAsim( 60,3,12,0.5,0.71)
EFAsim(100,3,12,0.5,0.71)
EFAsim(200,3,12,0.5,0.71)
EFAsim(400,3,12,0.5,0.71)
EFAsim( 60,5,12,0.5,0.71)
EFAsim(100,5,12,0.5,0.71)
EFAsim(200,5,12,0.5,0.71)
EFAsim(400,5,12,0.5,0.71)

EFAsim( 60,3, 8,0.7,0.71)
EFAsim(100,3, 8,0.7,0.71)
EFAsim(200,3, 8,0.7,0.71)
EFAsim(400,3, 8,0.7,0.71)
EFAsim( 60,5, 8,0.7,0.71)
EFAsim(100,5, 8,0.7,0.71)
EFAsim(200,5, 8,0.7,0.71)
EFAsim(400,5, 8,0.7,0.71)

EFAsim( 60,3,12,0.7,0.71)
EFAsim(100,3,12,0.7,0.71)
EFAsim(200,3,12,0.7,0.71)
EFAsim(400,3,12,0.7,0.71)
EFAsim( 60,5,12,0.7,0.71)
EFAsim(100,5,12,0.7,0.71)
EFAsim(200,5,12,0.7,0.71)
EFAsim(400,5,12,0.7,0.71)

EFAsim( 60,1, 8,0.3,0.89)
EFAsim(100,1, 8,0.3,0.89)
EFAsim(200,1, 8,0.3,0.89)
EFAsim(400,1, 8,0.3,0.89)

EFAsim( 60,3, 8,0.3,0.89)
EFAsim(100,3, 8,0.3,0.89)
EFAsim(200,3, 8,0.3,0.89)
EFAsim(400,3, 8,0.3,0.89)
EFAsim( 60,5, 8,0.3,0.89)
EFAsim(100,5, 8,0.3,0.89)
EFAsim(200,5, 8,0.3,0.89)
EFAsim(400,5, 8,0.3,0.89)

EFAsim( 60,1,12,0.3,0.89)
EFAsim(100,1,12,0.3,0.89)
EFAsim(200,1,12,0.3,0.89)
EFAsim(400,1,12,0.3,0.89)

EFAsim( 60,3,12,0.3,0.89)
EFAsim(100,3,12,0.3,0.89)
EFAsim(200,3,12,0.3,0.89)
EFAsim(400,3,12,0.3,0.89)
EFAsim( 60,5,12,0.3,0.89)
EFAsim(100,5,12,0.3,0.89)
EFAsim(200,5,12,0.3,0.89)
EFAsim(400,5,12,0.3,0.89)

EFAsim( 60,3, 8,0.5,0.89)
EFAsim(100,3, 8,0.5,0.89)
EFAsim(200,3, 8,0.5,0.89)
EFAsim(400,3, 8,0.5,0.89)
EFAsim( 60,5, 8,0.5,0.89)
EFAsim(100,5, 8,0.5,0.89)
EFAsim(200,5, 8,0.5,0.89)
EFAsim(400,5, 8,0.5,0.89)

EFAsim( 60,3,12,0.5,0.89)
EFAsim(100,3,12,0.5,0.89)
EFAsim(200,3,12,0.5,0.89)
EFAsim(400,3,12,0.5,0.89)
EFAsim( 60,5,12,0.5,0.89)
EFAsim(100,5,12,0.5,0.89)
EFAsim(200,5,12,0.5,0.89)
EFAsim(400,5,12,0.5,0.89)

EFAsim( 60,3, 8,0.7,0.89)
EFAsim(100,3, 8,0.7,0.89)
EFAsim(200,3, 8,0.7,0.89)
EFAsim(400,3, 8,0.7,0.89)
EFAsim( 60,5, 8,0.7,0.89)
EFAsim(100,5, 8,0.7,0.89)
EFAsim(200,5, 8,0.7,0.89)
EFAsim(400,5, 8,0.7,0.89)

EFAsim( 60,3,12,0.7,0.89)
EFAsim(100,3,12,0.7,0.89)
EFAsim(200,3,12,0.7,0.89)
EFAsim(400,3,12,0.7,0.89)
EFAsim( 60,5,12,0.7,0.89)
EFAsim(100,5,12,0.7,0.89)
EFAsim(200,5,12,0.7,0.89)
EFAsim(400,5,12,0.7,0.89)
