library(dplyr)

K <- "1K"

df <- read.table(paste0('./results/anovaDataFixed',K,'.dat'),
                 header = TRUE,
                 sep = "\t")

# path 7 ################################################################################################################################################################

df1 <- filter(df, efa == 1 & rot == 1 & length == 2) %>%
       mutate(biasMean = NA,
              RMSE = NA)

for (i in 1:nrow(df1)) {
  
  if(df1$nEta[i] > 1) {
    
    #set index variables (levels/conditions)
    n <- as.numeric(df1$n[i])
    nEta <- as.numeric(df1$nEta[i])
    pEta <- as.numeric(df1$pEta[i]) 
    rEta <- as.numeric(df1$rEta[i])
    l <- as.numeric(df1$l[i])
    d <- as.numeric(df1$d[i])
    
    #read in estimated factor correlations
    matPhi <- read.table(paste0('./path 7/',K,'/factor correlations (drop)/','nEta',nEta,'pEta',pEta,'rEta',rEta,'l',l,'n',n,'d',d,'.dat'),
                         header = FALSE,
                         sep = '\t')
    
    estPhi <- t(matPhi)[lower.tri(t(matPhi))]
    
    #generated factor correlations
    genPhi <- rep(as.numeric(df1$rEta[i]), 
                  times = length(estPhi))
    
    #compute bias
    bias <- estPhi - genPhi
    df1$biasMean[i] <- mean(bias)
    
    #RMSE of parameter estimates
    df1$RMSE[i] <- sqrt(mean(bias^2))
      
  }
  
}

# path 8 ################################################################################################################################################################

df2 <- filter(df, efa == 1 & rot == 2 & length == 2) %>%
  mutate(biasMean = NA,
         RMSE = NA)

for (i in 1:nrow(df2)) {
  
  if(df2$nEta[i] > 1) {
    
    #set index variables (levels/conditions)
    n <- as.numeric(df2$n[i])
    nEta <- as.numeric(df2$nEta[i])
    pEta <- as.numeric(df2$pEta[i]) 
    rEta <- as.numeric(df2$rEta[i])
    l <- as.numeric(df2$l[i])
    d <- as.numeric(df2$d[i])
    
    #read in estimated factor correlations
    matPhi <- read.table(paste0('./path 8/',K,'/factor correlations (drop)/','nEta',nEta,'pEta',pEta,'rEta',rEta,'l',l,'n',n,'d',d,'.dat'),
                         header = FALSE,
                         sep = '\t')
    
    estPhi <- t(matPhi)[lower.tri(t(matPhi))]
    
    #generated factor correlations
    genPhi <- rep(as.numeric(df2$rEta[i]), 
                  times = length(estPhi))
    
    #compute bias
    bias <- estPhi - genPhi
    df2$biasMean[i] <- mean(bias)
    
    #RMSE of parameter estimates
    df2$RMSE[i] <- sqrt(mean(bias^2))
    
  }
  
}

# path 9 ################################################################################################################################################################

df3 <- filter(df, efa == 1 & rot == 3 & length == 2) %>%
  mutate(biasMean = NA,
         RMSE = NA)

for (i in 1:nrow(df3)) {
  
  if(df3$nEta[i] == 3) {
    
    #compute bias
    bias <- rep(-as.numeric(df3$rEta[i]), 
                times = 3)
    df3$biasMean[i] <- mean(bias)
    
    #RMSE of parameter estimates
    df3$RMSE[i] <- sqrt(mean(bias^2))
    
  }
  
  if(df3$nEta[i] == 5) {
    
    #compute bias
    bias <- rep(-as.numeric(df3$rEta[i]), 
                times = 10)
    df3$biasMean[i] <- mean(bias)
    
    #RMSE of parameter estimates
    df3$RMSE[i] <- sqrt(mean(bias^2))
    
  }
  
}

# path 10 ################################################################################################################################################################

df4 <- filter(df, efa == 2 & rot == 1 & length == 2) %>%
  mutate(biasMean = NA,
         RMSE = NA)

for (i in 1:nrow(df4)) {
  
  if(df4$nEta[i] > 1) {
    
    #set index variables (levels/conditions)
    n <- as.numeric(df4$n[i])
    nEta <- as.numeric(df4$nEta[i])
    pEta <- as.numeric(df4$pEta[i]) 
    rEta <- as.numeric(df4$rEta[i])
    l <- as.numeric(df4$l[i])
    d <- as.numeric(df4$d[i])
    
    #read in estimated factor correlations
    matPhi <- read.table(paste0('./path 10/',K,'/factor correlations (drop)/','nEta',nEta,'pEta',pEta,'rEta',rEta,'l',l,'n',n,'d',d,'.dat'),
                         header = FALSE,
                         sep = '\t')
    
    estPhi <- t(matPhi)[lower.tri(t(matPhi))]
    
    #generated factor correlations
    genPhi <- rep(as.numeric(df4$rEta[i]), 
                  times = length(estPhi))
    
    #compute bias
    bias <- estPhi - genPhi
    df4$biasMean[i] <- mean(bias)
    
    #RMSE of parameter estimates
    df4$RMSE[i] <- sqrt(mean(bias^2))
    
  }
  
}

# path 11 ################################################################################################################################################################

df5 <- filter(df, efa == 2 & rot == 2 & length == 2) %>%
  mutate(biasMean = NA,
         RMSE = NA)

for (i in 1:nrow(df5)) {
  
  if(df5$nEta[i] > 1) {
    
    #set index variables (levels/conditions)
    n <- as.numeric(df5$n[i])
    nEta <- as.numeric(df5$nEta[i])
    pEta <- as.numeric(df5$pEta[i]) 
    rEta <- as.numeric(df5$rEta[i])
    l <- as.numeric(df5$l[i])
    d <- as.numeric(df5$d[i])
    
    #read in estimated factor correlations
    matPhi <- read.table(paste0('./path 11/',K,'/factor correlations (drop)/','nEta',nEta,'pEta',pEta,'rEta',rEta,'l',l,'n',n,'d',d,'.dat'),
                         header = FALSE,
                         sep = '\t')
    
    estPhi <- t(matPhi)[lower.tri(t(matPhi))]
    
    #generated factor correlations
    genPhi <- rep(as.numeric(df5$rEta[i]), 
                  times = length(estPhi))
    
    #compute bias
    bias <- estPhi - genPhi
    df5$biasMean[i] <- mean(bias)
    
    #RMSE of parameter estimates
    df5$RMSE[i] <- sqrt(mean(bias^2))
    
  }
  
}

# path 12 ################################################################################################################################################################

df6 <- filter(df, efa == 2 & rot == 3 & length == 2) %>%
  mutate(biasMean = NA,
         RMSE = NA)

for (i in 1:nrow(df6)) {
  
  if(df6$nEta[i] == 3) {
    
    #compute bias
    bias <- rep(-as.numeric(df6$rEta[i]), 
                times = 3)
    df6$biasMean[i] <- mean(bias)
    
    #RMSE of parameter estimates
    df6$RMSE[i] <- sqrt(mean(bias^2))
    
  }
  
  if(df6$nEta[i] == 5) {
    
    #compute bias
    bias <- rep(-as.numeric(df6$rEta[i]), 
                times = 10)
    df6$biasMean[i] <- mean(bias)
    
    #RMSE of parameter estimates
    df6$RMSE[i] <- sqrt(mean(bias^2))
    
  }
  
}

# combine data sets #####################################################################################################################################################

dfNew <- bind_rows(df1,df2,df3,df4,df5,df6)

dfNew <- dfNew %>% arrange(n,nEta,pEta,rEta,l,efa,rot,length,d)

#write out new data set
write.table(dfNew,
            file = paste0('./results/phi/phiShort',K,'.dat'),
            col.names = TRUE,
            row.names = FALSE,
            sep = '\t')
