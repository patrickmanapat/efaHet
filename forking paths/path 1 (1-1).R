# SETUP ##################################################################

#packages
library(psych) #for fa.parallel
library(gtools) #for permutations

#constants
rep <- 1000

# TESTING ################################################################

# #test runs
# n <- 200
# nEta <- 1
# pEta <- 12
# rEta <- 0.7
# l <- 0.71
# d <- 857
# 
# data <- read.table(paste('C:/Users/Patrick Manapat/Dropbox (ASU)/Research/Dissertation/Paper 3/Simulation/Data Generation/','nEta',nEta,'pEta',pEta,'rEta',rEta,'l',l,'n',n,'d',d,'.dat',
#                          sep=''),
#                    header = FALSE,
#                    sep = '\t')

# SIMULATION #############################################################

#begin sim function
EFAsim <- function(n,nEta,pEta,rEta,l){

#open sim loop
for (d in 1:rep){

#read in generated dataset
data <- read.table(paste('E:/ASU/Simulations/RDF QRP Sim/Data Generation/y/','nEta',nEta,'pEta',pEta,'rEta',rEta,'l',l,'n',n,'d',d,'.dat',sep=''),
                   header = FALSE,
                   sep = '\t')

# ANALYSIS (given correct number of factors) #############################

if(nEta == 1){
  
########## Target Rotation ##########

#create target matrix
targ <- matrix(0,
               nrow = pEta*nEta,
               ncol = nEta)

  for (i in 1:nEta) {
    
    targ[(pEta*(i-1)+1):(pEta*i),i] <- NA
    
  }

#EFA with target rotation
efaTarget <- fa(data,
                nfactors = nEta,
                rotate = "TargetQ",
                Target = list(targ), #target matrix must be a list
                fm = "ols")

########## EFA/PCA ##########

# (1-1) oblique quartimin rotation
efa <- fa(data,
          nfactors = nEta,
          rotate = "quartimin",
          fm = "ols")
  
# Collect Loadings and Check for Reflection #####################

#generating loadings
genLoad <- matrix(l,
                  nrow = pEta*nEta,
                  ncol = nEta)

#target loadings
target <- as.matrix(efaTarget$loadings[,1])

#reflection correction for target loadings
targReflect <- 0

if(all(target < 0)){
  
  target <- target * -1
  targReflect <- targReflect + 1
  
}

#estimated loadings
estimated <- as.matrix(efa$loadings[,1])

#reflection correction for estimated loadings
estReflect <- 0

if (all(estimated < 0)){
  
  estimated <- estimated * -1
  estReflect <- estReflect + 1
  
}

# Compute Congruence Coefficient (K) #####################################

#K for target loadings
targCompK <- matrix(NA,
                    nrow = 3,
                    ncol = nEta)
  
targCompK[1,1] <- sum(genLoad*target)
targCompK[2,1] <- sqrt(sum(genLoad^2)*sum(target^2))
targCompK[3,1] <- targCompK[1,1]/targCompK[2,1]
  
targK <- targCompK[3,1]

#K for estimated loadings
estCompK <- matrix(NA,
                   nrow = 3,
                   ncol = nEta)
  
estCompK[1,1] <- sum(genLoad*estimated)
estCompK[2,1] <- sqrt(sum(genLoad^2)*sum(estimated^2))
estCompK[3,1] <- estCompK[1,1]/estCompK[2,1]

estK <- estCompK[3,1]

# NA Output for 1-Factor #################################################

#scalars
countTargPerm <- NA
countEstPerm <- NA
targResMin <- NA
targAbsMin <- NA
targResAbs <- NA
estResMin <- NA
estAbsMin <- NA
estResAbs <- NA

#matrices
targPhi <- NA 
efaPhi <- NA
targTop1 <- NA
targTop2 <- NA
estTop1 <- NA
estTop2 <- NA
targLoadRes <- NA
targLoadAbs <- NA
estLoadRes <- NA
estLoadAbs <- NA

} else {
  
# ANALYSIS (given correct number of factors) #############################
  
########## Target Rotation ##########

#create target matrix
targ <- matrix(0,
               nrow = pEta*nEta,
               ncol = nEta)

for (i in 1:nEta) {
  
  targ[(pEta*(i-1)+1):(pEta*i),i] <- NA
  
}

#EFA with target rotation
efaTarget <- fa(data,
                nfactors = nEta,
                rotate = "TargetQ",
                Target = list(targ), #target matrix must be a list
                fm = "ols")

targPhi <- efaTarget$Phi

########## EFA/PCA ##########

# (1-1) oblique quartimin rotation
efa <- fa(data,
          nfactors = nEta,
          rotate = "quartimin",
          fm = "ols")

efaPhi <- efa$Phi #only for nEta > 1

# Collect Loadings and Check for Reflection #####################

#generating loadings
genLoad <- matrix(0,
                  nrow = pEta*nEta,
                  ncol = nEta)

for (i in 1:nEta) {
  
  genLoad[(pEta*(i-1)+1):(pEta*i),i] <- l
  
}

#target loadings
targLoad <- matrix(NA,
                   nrow = pEta*nEta,
                   ncol = nEta)

targTop1 <- matrix(NA,
                   nrow = pEta/2,
                   ncol = nEta)

targTop2 <- matrix(NA,
                   nrow = pEta/2,
                   ncol = nEta)

for (i in 1:nEta) {
  
  targLoad[,i] <- efaTarget$loadings[,i]
  targTop1[,i] <- tail(sort(abs(targLoad[,i])),pEta/2)
  targTop2[,i] <- head(sort(targLoad[,i]),pEta/2)
  
}

#reflection correction for target loadings
targTop1 <- apply(targTop1, 2, rev)
targReflect <- 0

for (i in 1:nEta) {
  if (mean(targTop1[,i] + targTop2[,i]) == 0){
    
    targLoad[,i] <- targLoad[,i] * -1
    targReflect <- targReflect + 1
    
  }
}

#estimated loadings
estLoad <- matrix(NA,
                  nrow = pEta*nEta,
                  ncol = nEta)

estTop1 <- matrix(NA,
                  nrow = pEta/2,
                  ncol = nEta)

estTop2 <- matrix(NA,
                  nrow = pEta/2,
                  ncol = nEta)

for (i in 1:nEta) {
  
  estLoad[,i] <- efa$loadings[,i]
  estTop1[,i] <- tail(sort(abs(estLoad[,i])),pEta/2)
  estTop2[,i] <- head(sort(estLoad[,i]),pEta/2)
  
}

#reflection correction for estimated loadings
estTop1 <- apply(estTop1, 2, rev)
estReflect <- 0

for (i in 1:nEta) {
  if (mean(estTop1[,i] + estTop2[,i]) == 0){
    
    estLoad[,i] <- estLoad[,i] * -1
    estReflect <- estReflect + 1
    
  }
}

# Check for Consistent Column Ordering ###################################

perm <- permutations(n = nEta, r = nEta, repeats.allowed = FALSE)

########## Target Loadings ##########

#list that stores every possible column permutation of the loading matrix
targLoadPerm <- list()
for (i in 1:nrow(perm)) {
  
  targLoadPerm[[i]] <- targLoad[,perm[i,]]
  
}

countTargPerm <- length(targLoadPerm)

#lists for SS residual (SSresid) and mean absolute bias (MAB)
targLoadRes <- list()
targLoadAbs <- list()
for (i in 1:nrow(perm)) {
  
  targLoadRes[i] <- sum((genLoad-targLoadPerm[[i]])^2)
  targLoadAbs[i] <- mean(abs(genLoad-targLoadPerm[[i]]))
  
}

#identify matrix with minimum SSresid and MAB
targResMin <- which.min(targLoadRes)
targAbsMin <- which.min(targLoadAbs)
#do SSresid and MAB agree?
if (targResMin == targAbsMin){targResAbs <- 1} else {targResAbs <- 0}

#final loading matrix
target <- targLoadPerm[[targResMin]]

########## Estimated Loadings ##########

#list that stores every possible column permutation of the loading matrix
estLoadPerm <- list()
for (i in 1:nrow(perm)) {
  
  estLoadPerm[[i]] <- estLoad[,perm[i,]]
  
}

countEstPerm <- length(estLoadPerm)

#lists for SS residual (SSresid) and mean absolute bias (MAB)
estLoadRes <- list()
estLoadAbs <- list()
for (i in 1:nrow(perm)) {
  
  estLoadRes[i] <- sum((genLoad-estLoadPerm[[i]])^2)
  estLoadAbs[i] <- mean(abs(genLoad-estLoadPerm[[i]]))
  
}

#identify matrix with minimum SSresid and MAB
estResMin <- which.min(estLoadRes)
estAbsMin <- which.min(estLoadAbs)
#do SSresid and MAB agree?
if (estResMin == estAbsMin){estResAbs <- 1} else {estResAbs <- 0}

#final loading matrix
estimated <- estLoadPerm[[estResMin]]

# Compute Congruence Coefficient (K) #####################################

#K for target loadings
targCompK <- matrix(NA,
                    nrow = 3,
                    ncol = nEta)

for (i in 1:nEta) {
  
  targCompK[1,i] <- sum(genLoad[,i]*target[,i])
  targCompK[2,i] <- sqrt(sum(genLoad[,i]^2)*sum(target[,i]^2))
  targCompK[3,i] <- targCompK[1,i]/targCompK[2,i]
  
}

targK <- mean(targCompK[3,])

#K for estimated loadings
estCompK <- matrix(NA,
                   nrow = 3,
                   ncol = nEta)

for (i in 1:nEta) {
  
  estCompK[1,i] <- sum(genLoad[,i]*estimated[,i])
  estCompK[2,i] <- sqrt(sum(genLoad[,i]^2)*sum(estimated[,i]^2))
  estCompK[3,i] <- estCompK[1,i]/estCompK[2,i]
  
}

estK <- mean(estCompK[3,])
  
}

# Collect Output #########################################################

#checks
outCheck <- cbind(n,nEta,pEta,rEta,l,targReflect,estReflect,countTargPerm,countEstPerm,targResMin,targAbsMin,targResAbs,estResMin,estAbsMin,estResAbs)

if(d == 1){finalCheck <- outCheck}
if(d > 1){finalCheck <- rbind(finalCheck,outCheck)}

#recovery measures
outRecov <- cbind(n,nEta,pEta,rEta,l,targK,estK)

if(d == 1){finalRecov <- outRecov}
if(d > 1){finalRecov <- rbind(finalRecov,outRecov)}

#objects
write.table(targ,
            paste('E:/ASU/Simulations/RDF QRP Sim/Estimation & Recovery/path 1/saved objects/target rotation/','nEta',nEta,'pEta',pEta,'rEta',rEta,'l',l,'n',n,'d',d,'.dat',sep=''),
            row.names = FALSE,
            col.names = FALSE,
            sep = '\t')

write.table(targPhi,
            paste('E:/ASU/Simulations/RDF QRP Sim/Estimation & Recovery/path 1/saved objects/factor correlations (target)/','nEta',nEta,'pEta',pEta,'rEta',rEta,'l',l,'n',n,'d',d,'.dat',sep=''),
            row.names = FALSE,
            col.names = FALSE,
            sep = '\t')

write.table(efaPhi,
            paste('E:/ASU/Simulations/RDF QRP Sim/Estimation & Recovery/path 1/saved objects/factor correlations/','nEta',nEta,'pEta',pEta,'rEta',rEta,'l',l,'n',n,'d',d,'.dat',sep=''),
            row.names = FALSE,
            col.names = FALSE,
            sep = '\t')

write.table(genLoad,
            paste('E:/ASU/Simulations/RDF QRP Sim/Estimation & Recovery/path 1/saved objects/generating loadings/','nEta',nEta,'pEta',pEta,'rEta',rEta,'l',l,'n',n,'d',d,'.dat',sep=''),
            row.names = FALSE,
            col.names = FALSE,
            sep = '\t')

write.table(targTop1,
            paste('E:/ASU/Simulations/RDF QRP Sim/Estimation & Recovery/path 1/saved objects/greatest absolute loadings (target)/','nEta',nEta,'pEta',pEta,'rEta',rEta,'l',l,'n',n,'d',d,'.dat',sep=''),
            row.names = FALSE,
            col.names = FALSE,
            sep = '\t')

write.table(targTop2,
            paste('E:/ASU/Simulations/RDF QRP Sim/Estimation & Recovery/path 1/saved objects/greatest negative loadings (target)/','nEta',nEta,'pEta',pEta,'rEta',rEta,'l',l,'n',n,'d',d,'.dat',sep=''),
            row.names = FALSE,
            col.names = FALSE,
            sep = '\t')

write.table(estTop1,
            paste('E:/ASU/Simulations/RDF QRP Sim/Estimation & Recovery/path 1/saved objects/greatest absolute loadings/','nEta',nEta,'pEta',pEta,'rEta',rEta,'l',l,'n',n,'d',d,'.dat',sep=''),
            row.names = FALSE,
            col.names = FALSE,
            sep = '\t')

write.table(estTop2,
            paste('E:/ASU/Simulations/RDF QRP Sim/Estimation & Recovery/path 1/saved objects/greatest negative loadings/','nEta',nEta,'pEta',pEta,'rEta',rEta,'l',l,'n',n,'d',d,'.dat',sep=''),
            row.names = FALSE,
            col.names = FALSE,
            sep = '\t')

write.table(targLoadRes,
            paste('E:/ASU/Simulations/RDF QRP Sim/Estimation & Recovery/path 1/saved objects/ss residual (target)/','nEta',nEta,'pEta',pEta,'rEta',rEta,'l',l,'n',n,'d',d,'.dat',sep=''),
            row.names = FALSE,
            col.names = FALSE,
            sep = '\t')

write.table(targLoadAbs,
            paste('E:/ASU/Simulations/RDF QRP Sim/Estimation & Recovery/path 1/saved objects/mean absolute bias (target)/','nEta',nEta,'pEta',pEta,'rEta',rEta,'l',l,'n',n,'d',d,'.dat',sep=''),
            row.names = FALSE,
            col.names = FALSE,
            sep = '\t')

write.table(target,
            paste('E:/ASU/Simulations/RDF QRP Sim/Estimation & Recovery/path 1/saved objects/estimated loadings (target)/','nEta',nEta,'pEta',pEta,'rEta',rEta,'l',l,'n',n,'d',d,'.dat',sep=''),
            row.names = FALSE,
            col.names = FALSE,
            sep = '\t')

write.table(estLoadRes,
            paste('E:/ASU/Simulations/RDF QRP Sim/Estimation & Recovery/path 1/saved objects/ss residual/','nEta',nEta,'pEta',pEta,'rEta',rEta,'l',l,'n',n,'d',d,'.dat',sep=''),
            row.names = FALSE,
            col.names = FALSE,
            sep = '\t')

write.table(estLoadAbs,
            paste('E:/ASU/Simulations/RDF QRP Sim/Estimation & Recovery/path 1/saved objects/mean absolute bias/','nEta',nEta,'pEta',pEta,'rEta',rEta,'l',l,'n',n,'d',d,'.dat',sep=''),
            row.names = FALSE,
            col.names = FALSE,
            sep = '\t')

write.table(estimated,
            paste('E:/ASU/Simulations/RDF QRP Sim/Estimation & Recovery/path 1/saved objects/estimated loadings/','nEta',nEta,'pEta',pEta,'rEta',rEta,'l',l,'n',n,'d',d,'.dat',sep=''),
            row.names = FALSE,
            col.names = FALSE,
            sep = '\t')

write.table(targCompK,
            paste('E:/ASU/Simulations/RDF QRP Sim/Estimation & Recovery/path 1/saved objects/congruence coefficient (target)/','nEta',nEta,'pEta',pEta,'rEta',rEta,'l',l,'n',n,'d',d,'.dat',sep=''),
            row.names = FALSE,
            col.names = FALSE,
            sep = '\t')

write.table(estCompK,
            paste('E:/ASU/Simulations/RDF QRP Sim/Estimation & Recovery/path 1/saved objects/congruence coefficient/','nEta',nEta,'pEta',pEta,'rEta',rEta,'l',l,'n',n,'d',d,'.dat',sep=''),
            row.names = FALSE,
            col.names = FALSE,
            sep = '\t')

#close sim loop
}

write.table(finalCheck, 
            paste('E:/ASU/Simulations/RDF QRP Sim/Estimation & Recovery/path 1/checks/','nEta',nEta,'pEta',pEta,'rEta',rEta,'l',l,'n',n,'.dat',sep=''),
            row.names = FALSE,
            col.names = TRUE,
            sep = '\t')

write.table(finalRecov, 
            paste('E:/ASU/Simulations/RDF QRP Sim/Estimation & Recovery/path 1/recovery/','nEta',nEta,'pEta',pEta,'rEta',rEta,'l',l,'n',n,'.dat',sep=''),
            row.names = FALSE,
            col.names = TRUE,
            sep = '\t')

#close sim function
}

# SIMULATION CONDITIONS ##################################################

#EFAsim(n,nEta,pEta,rEta,l)
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

EFAsim( 60,1, 8,0.5,0.45)
EFAsim(100,1, 8,0.5,0.45)
EFAsim(200,1, 8,0.5,0.45)
EFAsim(400,1, 8,0.5,0.45)
EFAsim( 60,3, 8,0.5,0.45)
EFAsim(100,3, 8,0.5,0.45)
EFAsim(200,3, 8,0.5,0.45)
EFAsim(400,3, 8,0.5,0.45)
EFAsim( 60,5, 8,0.5,0.45)
EFAsim(100,5, 8,0.5,0.45)
EFAsim(200,5, 8,0.5,0.45)
EFAsim(400,5, 8,0.5,0.45)

EFAsim( 60,1,12,0.5,0.45)
EFAsim(100,1,12,0.5,0.45)
EFAsim(200,1,12,0.5,0.45)
EFAsim(400,1,12,0.5,0.45)
EFAsim( 60,3,12,0.5,0.45)
EFAsim(100,3,12,0.5,0.45)
EFAsim(200,3,12,0.5,0.45)
EFAsim(400,3,12,0.5,0.45)
EFAsim( 60,5,12,0.5,0.45)
EFAsim(100,5,12,0.5,0.45)
EFAsim(200,5,12,0.5,0.45)
EFAsim(400,5,12,0.5,0.45)

EFAsim( 60,1, 8,0.7,0.45)
EFAsim(100,1, 8,0.7,0.45)
EFAsim(200,1, 8,0.7,0.45)
EFAsim(400,1, 8,0.7,0.45)
EFAsim( 60,3, 8,0.7,0.45)
EFAsim(100,3, 8,0.7,0.45)
EFAsim(200,3, 8,0.7,0.45)
EFAsim(400,3, 8,0.7,0.45)
EFAsim( 60,5, 8,0.7,0.45)
EFAsim(100,5, 8,0.7,0.45)
EFAsim(200,5, 8,0.7,0.45)
EFAsim(400,5, 8,0.7,0.45)

EFAsim( 60,1,12,0.7,0.45)
EFAsim(100,1,12,0.7,0.45)
EFAsim(200,1,12,0.7,0.45)
EFAsim(400,1,12,0.7,0.45)
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

EFAsim( 60,1, 8,0.5,0.71)
EFAsim(100,1, 8,0.5,0.71)
EFAsim(200,1, 8,0.5,0.71)
EFAsim(400,1, 8,0.5,0.71)
EFAsim( 60,3, 8,0.5,0.71)
EFAsim(100,3, 8,0.5,0.71)
EFAsim(200,3, 8,0.5,0.71)
EFAsim(400,3, 8,0.5,0.71)
EFAsim( 60,5, 8,0.5,0.71)
EFAsim(100,5, 8,0.5,0.71)
EFAsim(200,5, 8,0.5,0.71)
EFAsim(400,5, 8,0.5,0.71)

EFAsim( 60,1,12,0.5,0.71)
EFAsim(100,1,12,0.5,0.71)
EFAsim(200,1,12,0.5,0.71)
EFAsim(400,1,12,0.5,0.71)
EFAsim( 60,3,12,0.5,0.71)
EFAsim(100,3,12,0.5,0.71)
EFAsim(200,3,12,0.5,0.71)
EFAsim(400,3,12,0.5,0.71)
EFAsim( 60,5,12,0.5,0.71)
EFAsim(100,5,12,0.5,0.71)
EFAsim(200,5,12,0.5,0.71)
EFAsim(400,5,12,0.5,0.71)

EFAsim( 60,1, 8,0.7,0.71)
EFAsim(100,1, 8,0.7,0.71)
EFAsim(200,1, 8,0.7,0.71)
EFAsim(400,1, 8,0.7,0.71)
EFAsim( 60,3, 8,0.7,0.71)
EFAsim(100,3, 8,0.7,0.71)
EFAsim(200,3, 8,0.7,0.71)
EFAsim(400,3, 8,0.7,0.71)
EFAsim( 60,5, 8,0.7,0.71)
EFAsim(100,5, 8,0.7,0.71)
EFAsim(200,5, 8,0.7,0.71)
EFAsim(400,5, 8,0.7,0.71)

EFAsim( 60,1,12,0.7,0.71)
EFAsim(100,1,12,0.7,0.71)
EFAsim(200,1,12,0.7,0.71)
EFAsim(400,1,12,0.7,0.71)
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

EFAsim( 60,1, 8,0.5,0.89)
EFAsim(100,1, 8,0.5,0.89)
EFAsim(200,1, 8,0.5,0.89)
EFAsim(400,1, 8,0.5,0.89)
EFAsim( 60,3, 8,0.5,0.89)
EFAsim(100,3, 8,0.5,0.89)
EFAsim(200,3, 8,0.5,0.89)
EFAsim(400,3, 8,0.5,0.89)
EFAsim( 60,5, 8,0.5,0.89)
EFAsim(100,5, 8,0.5,0.89)
EFAsim(200,5, 8,0.5,0.89)
EFAsim(400,5, 8,0.5,0.89)

EFAsim( 60,1,12,0.5,0.89)
EFAsim(100,1,12,0.5,0.89)
EFAsim(200,1,12,0.5,0.89)
EFAsim(400,1,12,0.5,0.89)
EFAsim( 60,3,12,0.5,0.89)
EFAsim(100,3,12,0.5,0.89)
EFAsim(200,3,12,0.5,0.89)
EFAsim(400,3,12,0.5,0.89)
EFAsim( 60,5,12,0.5,0.89)
EFAsim(100,5,12,0.5,0.89)
EFAsim(200,5,12,0.5,0.89)
EFAsim(400,5,12,0.5,0.89)

EFAsim( 60,1, 8,0.7,0.89)
EFAsim(100,1, 8,0.7,0.89)
EFAsim(200,1, 8,0.7,0.89)
EFAsim(400,1, 8,0.7,0.89)
EFAsim( 60,3, 8,0.7,0.89)
EFAsim(100,3, 8,0.7,0.89)
EFAsim(200,3, 8,0.7,0.89)
EFAsim(400,3, 8,0.7,0.89)
EFAsim( 60,5, 8,0.7,0.89)
EFAsim(100,5, 8,0.7,0.89)
EFAsim(200,5, 8,0.7,0.89)
EFAsim(400,5, 8,0.7,0.89)

EFAsim( 60,1,12,0.7,0.89)
EFAsim(100,1,12,0.7,0.89)
EFAsim(200,1,12,0.7,0.89)
EFAsim(400,1,12,0.7,0.89)
EFAsim( 60,3,12,0.7,0.89)
EFAsim(100,3,12,0.7,0.89)
EFAsim(200,3,12,0.7,0.89)
EFAsim(400,3,12,0.7,0.89)
EFAsim( 60,5,12,0.7,0.89)
EFAsim(100,5,12,0.7,0.89)
EFAsim(200,5,12,0.7,0.89)
EFAsim(400,5,12,0.7,0.89)