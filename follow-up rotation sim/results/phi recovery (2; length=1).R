library(dplyr)

df <- read.table('./results/anovaData (mini).dat',
                 header = TRUE,
                 sep = "\t")

# path 1 ################################################################################################################################################################

df1 <- filter(df, efa == 1 & rot == 1 & length == 1) %>%
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
    matPhi <- read.table(paste0('./path 1 (mini)/factor correlations/','nEta',nEta,'pEta',pEta,'rEta',rEta,'l',l,'n',n,'d',d,'.dat'),
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

# path 2 ################################################################################################################################################################

df2 <- filter(df, efa == 1 & rot == 2 & length == 1) %>%
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
    matPhi <- read.table(paste0('./path 2 (mini)/factor correlations/','nEta',nEta,'pEta',pEta,'rEta',rEta,'l',l,'n',n,'d',d,'.dat'),
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

# path 3 ################################################################################################################################################################

df3 <- filter(df, efa == 1 & rot == 3 & length == 1) %>%
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

# combine data sets #####################################################################################################################################################

dfNew <- bind_rows(df1,df2,df3)

dfNew <- dfNew %>% arrange(n,nEta,pEta,rEta,l,efa,rot,length,d)

#write out new data set
write.table(dfNew,
            file = paste0('./results/phi/phiFull',K,'.dat'),
            col.names = TRUE,
            row.names = FALSE,
            sep = '\t')
