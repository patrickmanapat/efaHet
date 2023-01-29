# SETUP ##################################################################

#packages
library(psych) #for fa.parallel

#constants
minrep <- 4001
maxrep <- 5000
K <- "5K"

# TESTING ################################################################
# 
# #test runs
# n <- 60
# nEta <- 3
# pEta <- 12
# rEta <- 0.3
# l <- 0.45
# d <- 4001
# PA <- 1
# EFA <- 2
# rot <- 1
# 
# data <- read.table(paste('C:/Users/Patrick Manapat/Dropbox (ASU)/Research/Dissertation/Paper 3/Simulation/Data Generation/','nEta',nEta,'pEta',pEta,'rEta',rEta,'l',l,'n',n,'d',d,'.dat',
#                          sep=''),
#                    header = FALSE,
#                    sep = '\t')

# SIMULATION #############################################################

#begin sim function
EFAsim <- function(n,nEta,pEta,rEta,l){

#open sim loop  
for (d in minrep:maxrep){

#read in generated dataset
data <- read.table(paste('../y',K,'/','nEta',nEta,'pEta',pEta,'rEta',rEta,'l',l,'n',n,'d',d,'.dat',sep=''),
                   header = FALSE,
                   sep = '\t')

# NUMBER OF FACTORS ######################################################

########## Parallel Analysis ########## 
pa <- fa.parallel(data, 
                  n.iter = 100, #Garrido et al. (2013)
                  sim = FALSE,
                  plot = FALSE)

PA_nFact <- pa$nfact
PA_nComp <- pa$ncomp

#did PA suggest correct number of factors?
if(PA_nFact == nEta){PA_nFactCorrect <- 1} else {PA_nFactCorrect <- 0}
if(PA_nComp == nEta){PA_nCompCorrect <- 1} else {PA_nCompCorrect <- 0}

########## K1 Criterion ##########
eig <- eigen(cor(data))

K1_nFact <- sum(eig$values > 1)

#did K1 suggest correct number of factors?
if(K1_nFact == nEta) {K1_nFactCorrect <- 1} else {K1_nFactCorrect <- 0}

# Collect Output #########################################################

out <- cbind(n,nEta,pEta,rEta,l,PA_nFact,PA_nComp,K1_nFact,PA_nFactCorrect,PA_nCompCorrect,K1_nFactCorrect)

#compile output from every iteration
if(d == minrep){final <- out}
if(d > minrep){final <- rbind(final,out)}

#close sim loop
}

write.table(final, 
            paste('./results',K,'/','nEta',nEta,'pEta',pEta,'rEta',rEta,'l',l,'n',n,'.dat',sep=''),
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