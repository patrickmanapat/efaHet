library(dplyr)

K <- "1K"

#read in full data set
dfRep <- read.table(paste0('./results/anovaData',K,'.dat'),
                    header = TRUE,
                    sep = "\t")

#select variables and create index for missing (missK)
dfRep2 <- dfRep %>% select(-c(mean25, max25))

#filter not missing
df1 <- filter(dfRep2, is.na(estK) == FALSE)

#filter missing
dfNA <- filter(dfRep2, is.na(estK) == TRUE)

#path 7
df7 <- filter(dfNA, efa == 1 & rot == 1)

for (i in 1:nrow(df7)) {
  
  #set index variables (levels/conditions)
  n <- as.numeric(df7$n[i])
  nEta <- as.numeric(df7$nEta[i])
  pEta <- as.numeric(df7$pEta[i]) 
  rEta <- as.numeric(df7$rEta[i])
  l <- as.numeric(df7$l[i])
  d <- as.numeric(df7$d[i])
  
  #read in generating loadings
  genLamb <- read.table(paste0('./path 7/',K,'/objects/generating loadings (drop)/','nEta',nEta,'pEta',pEta,'rEta',rEta,'l',l,'n',n,'d',d,'.dat'),
                        header = FALSE,
                        sep = '\t')
    dropCol <- which(colSums(genLamb != 0) == 0)
    genLamb <- genLamb[,-dropCol]
  
  #read in estimated loadings    
  estLamb <- read.table(paste0('./path 7/',K,'/objects/estimated loadings (drop)/','nEta',nEta,'pEta',pEta,'rEta',rEta,'l',l,'n',n,'d',d,'.dat'),
                        header = FALSE,
                        sep = '\t')
    estLamb <- estLamb[,-dropCol]
    
  #compute congruence coefficient (K)
  estCompK <- matrix(NA,
                     nrow = 3,
                     ncol = nEta-1)
  
  for (j in 1:(nEta-1)) {
    
    estCompK[1,j] <- sum(genLamb[,j]*estLamb[,j])
    estCompK[2,j] <- sqrt(sum(genLamb[,j]^2)*sum(estLamb[,j]^2))
    estCompK[3,j] <- estCompK[1,j]/estCompK[2,j]
    
  }
  
  df7$estK[i] <- mean(estCompK[3,])
    
}

#path 8
df8 <- filter(dfNA, efa == 1 & rot == 2)

for (i in 1:nrow(df8)) {
  
  #set index variables (levels/conditions)
  n <- as.numeric(df8$n[i])
  nEta <- as.numeric(df8$nEta[i])
  pEta <- as.numeric(df8$pEta[i]) 
  rEta <- as.numeric(df8$rEta[i])
  l <- as.numeric(df8$l[i])
  d <- as.numeric(df8$d[i])
  
  #read in generating loadings
  genLamb <- read.table(paste0('./path 8/',K,'/objects/generating loadings (drop)/','nEta',nEta,'pEta',pEta,'rEta',rEta,'l',l,'n',n,'d',d,'.dat'),
                        header = FALSE,
                        sep = '\t')
  dropCol <- which(colSums(genLamb != 0) == 0)
  genLamb <- genLamb[,-dropCol]
  
  #read in estimated loadings    
  estLamb <- read.table(paste0('./path 8/',K,'/objects/estimated loadings (drop)/','nEta',nEta,'pEta',pEta,'rEta',rEta,'l',l,'n',n,'d',d,'.dat'),
                        header = FALSE,
                        sep = '\t')
  estLamb <- estLamb[,-dropCol]
  
  #compute congruence coefficient (K)
  estCompK <- matrix(NA,
                     nrow = 3,
                     ncol = nEta-1)
  
  for (j in 1:(nEta-1)) {
    
    estCompK[1,j] <- sum(genLamb[,j]*estLamb[,j])
    estCompK[2,j] <- sqrt(sum(genLamb[,j]^2)*sum(estLamb[,j]^2))
    estCompK[3,j] <- estCompK[1,j]/estCompK[2,j]
    
  }
  
  df8$estK[i] <- mean(estCompK[3,])
  
}

#path 9
df9 <- filter(dfNA, efa == 1 & rot == 3)

for (i in 1:nrow(df9)) {
  
  #set index variables (levels/conditions)
  n <- as.numeric(df9$n[i])
  nEta <- as.numeric(df9$nEta[i])
  pEta <- as.numeric(df9$pEta[i]) 
  rEta <- as.numeric(df9$rEta[i])
  l <- as.numeric(df9$l[i])
  d <- as.numeric(df9$d[i])
  
  #read in generating loadings
  genLamb <- read.table(paste0('./path 9/',K,'/objects/generating loadings (drop)/','nEta',nEta,'pEta',pEta,'rEta',rEta,'l',l,'n',n,'d',d,'.dat'),
                        header = FALSE,
                        sep = '\t')
  dropCol <- which(colSums(genLamb != 0) == 0)
  genLamb <- genLamb[,-dropCol]
  
  #read in estimated loadings    
  estLamb <- read.table(paste0('./path 9/',K,'/objects/estimated loadings (drop)/','nEta',nEta,'pEta',pEta,'rEta',rEta,'l',l,'n',n,'d',d,'.dat'),
                        header = FALSE,
                        sep = '\t')
  estLamb <- estLamb[,-dropCol]
  
  #compute congruence coefficient (K)
  estCompK <- matrix(NA,
                     nrow = 3,
                     ncol = nEta-1)
  
  for (j in 1:(nEta-1)) {
    
    estCompK[1,j] <- sum(genLamb[,j]*estLamb[,j])
    estCompK[2,j] <- sqrt(sum(genLamb[,j]^2)*sum(estLamb[,j]^2))
    estCompK[3,j] <- estCompK[1,j]/estCompK[2,j]
    
  }
  
  df9$estK[i] <- mean(estCompK[3,])
  
}

#path 10
df10 <- filter(dfNA, efa == 2 & rot == 1)

for (i in 1:nrow(df10)) {
  
  #set index variables (levels/conditions)
  n <- as.numeric(df10$n[i])
  nEta <- as.numeric(df10$nEta[i])
  pEta <- as.numeric(df10$pEta[i]) 
  rEta <- as.numeric(df10$rEta[i])
  l <- as.numeric(df10$l[i])
  d <- as.numeric(df10$d[i])
  
  #read in generating loadings
  genLamb <- read.table(paste0('./path 10/',K,'/objects/generating loadings (drop)/','nEta',nEta,'pEta',pEta,'rEta',rEta,'l',l,'n',n,'d',d,'.dat'),
                        header = FALSE,
                        sep = '\t')
  dropCol <- which(colSums(genLamb != 0) == 0)
  genLamb <- genLamb[,-dropCol]
  
  #read in estimated loadings    
  estLamb <- read.table(paste0('./path 10/',K,'/objects/estimated loadings (drop)/','nEta',nEta,'pEta',pEta,'rEta',rEta,'l',l,'n',n,'d',d,'.dat'),
                        header = FALSE,
                        sep = '\t')
  estLamb <- estLamb[,-dropCol]
  
  #compute congruence coefficient (K)
  estCompK <- matrix(NA,
                     nrow = 3,
                     ncol = nEta-1)
  
  for (j in 1:(nEta-1)) {
    
    estCompK[1,j] <- sum(genLamb[,j]*estLamb[,j])
    estCompK[2,j] <- sqrt(sum(genLamb[,j]^2)*sum(estLamb[,j]^2))
    estCompK[3,j] <- estCompK[1,j]/estCompK[2,j]
    
  }
  
  df10$estK[i] <- mean(estCompK[3,])
  
}

#path 11
df11 <- filter(dfNA, efa == 2 & rot == 2)

for (i in 1:nrow(df11)) {
  
  #set index variables (levels/conditions)
  n <- as.numeric(df11$n[i])
  nEta <- as.numeric(df11$nEta[i])
  pEta <- as.numeric(df11$pEta[i]) 
  rEta <- as.numeric(df11$rEta[i])
  l <- as.numeric(df11$l[i])
  d <- as.numeric(df11$d[i])
  
  #read in generating loadings
  genLamb <- read.table(paste0('./path 11/',K,'/objects/generating loadings (drop)/','nEta',nEta,'pEta',pEta,'rEta',rEta,'l',l,'n',n,'d',d,'.dat'),
                        header = FALSE,
                        sep = '\t')
  dropCol <- which(colSums(genLamb != 0) == 0)
  genLamb <- genLamb[,-dropCol]
  
  #read in estimated loadings    
  estLamb <- read.table(paste0('./path 11/',K,'/objects/estimated loadings (drop)/','nEta',nEta,'pEta',pEta,'rEta',rEta,'l',l,'n',n,'d',d,'.dat'),
                        header = FALSE,
                        sep = '\t')
  estLamb <- estLamb[,-dropCol]
  
  #compute congruence coefficient (K)
  estCompK <- matrix(NA,
                     nrow = 3,
                     ncol = nEta-1)
  
  for (j in 1:(nEta-1)) {
    
    estCompK[1,j] <- sum(genLamb[,j]*estLamb[,j])
    estCompK[2,j] <- sqrt(sum(genLamb[,j]^2)*sum(estLamb[,j]^2))
    estCompK[3,j] <- estCompK[1,j]/estCompK[2,j]
    
  }
  
  df11$estK[i] <- mean(estCompK[3,])
  
}

#path 12
df12 <- filter(dfNA, efa == 2 & rot == 3)

for (i in 1:nrow(df12)) {
  
  #set index variables (levels/conditions)
  n <- as.numeric(df12$n[i])
  nEta <- as.numeric(df12$nEta[i])
  pEta <- as.numeric(df12$pEta[i]) 
  rEta <- as.numeric(df12$rEta[i])
  l <- as.numeric(df12$l[i])
  d <- as.numeric(df12$d[i])
  
  #read in generating loadings
  genLamb <- read.table(paste0('./path 12/',K,'/objects/generating loadings (drop)/','nEta',nEta,'pEta',pEta,'rEta',rEta,'l',l,'n',n,'d',d,'.dat'),
                        header = FALSE,
                        sep = '\t')
  dropCol <- which(colSums(genLamb != 0) == 0)
  genLamb <- genLamb[,-dropCol]
  
  #read in estimated loadings    
  estLamb <- read.table(paste0('./path 12/',K,'/objects/estimated loadings (drop)/','nEta',nEta,'pEta',pEta,'rEta',rEta,'l',l,'n',n,'d',d,'.dat'),
                        header = FALSE,
                        sep = '\t')
  estLamb <- estLamb[,-dropCol]
  
  #compute congruence coefficient (K)
  estCompK <- matrix(NA,
                     nrow = 3,
                     ncol = nEta-1)
  
  for (j in 1:(nEta-1)) {
    
    estCompK[1,j] <- sum(genLamb[,j]*estLamb[,j])
    estCompK[2,j] <- sqrt(sum(genLamb[,j]^2)*sum(estLamb[,j]^2))
    estCompK[3,j] <- estCompK[1,j]/estCompK[2,j]
    
  }
  
  df12$estK[i] <- mean(estCompK[3,])
  
}

dfNew <- bind_rows(df1,df7,df8,df9,df10,df11,df12)

dfNew <- dfNew %>% arrange(n,nEta,pEta,rEta,l,efa,rot,length,d)

#write out new data set
write.table(dfNew,
            file = 'C:/Users/Patrick Manapat/Dropbox (ASU)/Research/Dissertation/Paper 3/Simulation/Estimation & Recovery/recovery/anovaDataFixed.dat',
            col.names = TRUE,
            row.names = FALSE,
            sep = '\t')

#dfNA %>% group_by(length) %>% summarize(n())
