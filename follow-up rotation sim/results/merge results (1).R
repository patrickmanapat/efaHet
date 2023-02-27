library(dplyr)

minrep <- 1
maxrep <- 5000

# Path 1 #################################################################

#list files in a folder
p1 <- list.files(path = paste0('./path 1/',K,'/recovery/'),
                 pattern = ".dat",
                 full.names = TRUE)

#row bind files in same folder (stack data sets)
p1DF <- do.call("rbind", lapply(p1, read.table, header = TRUE))
  p1DF <- p1DF[order(p1DF$nEta, p1DF$pEta, p1DF$rEta, p1DF$l, p1DF$n),] #put df in ascending order of factors
  row.names(p1DF) <- 1:nrow(p1DF) #put row numbers back in ascending order
  
#create variables for analysis choices
p1DF$efa <- 1
p1DF$rot <- 1
p1DF$length <- 1
p1DF$mean25 <- NA
p1DF$max25 <- NA
p1DF$d <- rep(minrep:maxrep, times = 168)

#change rEta for one factor to NA
p1DF <- within(p1DF, rEta[nEta == 1] <- 'None')

# Path 2 #################################################################
  
#list files in a folder
p2 <- list.files(path = paste0('./path 2/',K,'/recovery/'),
                 pattern = ".dat",
                 full.names = TRUE)

#row bind files in same folder (stack data sets)
p2DF <- do.call("rbind", lapply(p2, read.table, header = TRUE))
  p2DF <- p2DF[order(p2DF$nEta, p2DF$pEta, p2DF$rEta, p2DF$l, p2DF$n),] #put df in ascending order of factors
  row.names(p2DF) <- 1:nrow(p2DF) #put row numbers back in ascending order

#create variables for analysis choices
p2DF$efa <- 1
p2DF$rot <- 2
p2DF$length <- 1
p2DF$mean25 <- NA
p2DF$max25 <- NA
p2DF$d <- rep(minrep:maxrep, times = 168)

#change rEta for one factor to NA
p2DF <- within(p2DF, rEta[nEta == 1] <- 'None')
  
# Path 3 #################################################################

#list files in a folder
p3 <- list.files(path = paste0('./path 3/',K,'/recovery/'),
                 pattern = ".dat",
                 full.names = TRUE)

#row bind files in same folder (stack data sets)
p3DF <- do.call("rbind", lapply(p3, read.table, header = TRUE))
  p3DF <- p3DF[order(p3DF$nEta, p3DF$pEta, p3DF$rEta, p3DF$l, p3DF$n),] #put df in ascending order of factors
  row.names(p3DF) <- 1:nrow(p3DF) #put row numbers back in ascending order

#create variables for analysis choices
p3DF$efa <- 1
p3DF$rot <- 3
p3DF$length <- 1
p3DF$mean25 <- NA
p3DF$max25 <- NA
p3DF$d <- rep(minrep:maxrep, times = 168)

#change rEta for one factor to NA
p3DF <- within(p3DF, rEta[nEta == 1] <- 'None')
  
# Path 4 #################################################################
  
#list files in a folder
p4 <- list.files(path = paste0('./path 4/',K,'/recovery/'),
                 pattern = ".dat",
                 full.names = TRUE)

#row bind files in same folder (stack data sets)
p4DF <- do.call("rbind", lapply(p4, read.table, header = TRUE))
  p4DF <- p4DF[order(p4DF$nEta, p4DF$pEta, p4DF$rEta, p4DF$l, p4DF$n),] #put df in ascending order of factors
  row.names(p4DF) <- 1:nrow(p4DF) #put row numbers back in ascending order

#create variables for analysis choices
p4DF$efa <- 2
p4DF$rot <- 1
p4DF$length <- 1
p4DF$mean25 <- NA
p4DF$max25 <- NA
p4DF$d <- rep(minrep:maxrep, times = 168)

#change rEta for one factor to NA
p4DF <- within(p4DF, rEta[nEta == 1] <- 'None')
  
# Path 5 #################################################################
  
#list files in a folder
p5 <- list.files(path = paste0('./path 5/',K,'/recovery/'),
                 pattern = ".dat",
                 full.names = TRUE)

#row bind files in same folder (stack data sets)
p5DF <- do.call("rbind", lapply(p5, read.table, header = TRUE))
  p5DF <- p5DF[order(p5DF$nEta, p5DF$pEta, p5DF$rEta, p5DF$l, p5DF$n),] #put df in ascending order of factors
  row.names(p5DF) <- 1:nrow(p5DF) #put row numbers back in ascending order

#create variables for analysis choices
p5DF$efa <- 2
p5DF$rot <- 2
p5DF$length <- 1
p5DF$mean25 <- NA
p5DF$max25 <- NA
p5DF$d <- rep(minrep:maxrep, times = 168)

#change rEta for one factor to NA
p5DF <- within(p5DF, rEta[nEta == 1] <- 'None')
  
# Path 6 #################################################################

#list files in a folder
p6 <- list.files(path = paste0('./path 6/',K,'/recovery/'),
                 pattern = ".dat",
                 full.names = TRUE)

#row bind files in same folder (stack data sets)
p6DF <- do.call("rbind", lapply(p6, read.table, header = TRUE))
  p6DF <- p6DF[order(p6DF$nEta, p6DF$pEta, p6DF$rEta, p6DF$l, p6DF$n),] #put df in ascending order of factors
  row.names(p6DF) <- 1:nrow(p6DF) #put row numbers back in ascending order

#create variables for analysis choices
p6DF$efa <- 2
p6DF$rot <- 3
p6DF$length <- 1
p6DF$mean25 <- NA
p6DF$max25 <- NA
p6DF$d <- rep(minrep:maxrep, times = 168)

#change rEta for one factor to NA
p6DF <- within(p6DF, rEta[nEta == 1] <- 'None')

# Path 7 #################################################################

#list files in a folder
p7 <- list.files(path = paste0('./path 7/',K,'/recovery/'),
                 pattern = ".dat",
                 full.names = TRUE)

#row bind files in same folder (stack data sets)
p7DF <- do.call("rbind", lapply(p7, read.table, header = TRUE))
p7DF <- p7DF[order(p7DF$nEta, p7DF$pEta, p7DF$rEta, p7DF$l, p7DF$n),] #put df in ascending order of factors
row.names(p7DF) <- 1:nrow(p7DF) #put row numbers back in ascending order

#create variables for analysis choices
p7DF$efa <- 1
p7DF$rot <- 1
p7DF$length <- 2

#move mean25 and max25 to end of data frame
p7DF <- p7DF %>% relocate(mean25, .after = last_col()) %>% relocate(max25, .after = last_col())

p7DF$d <- rep(minrep:maxrep, times = 168)

#change rEta for one factor to NA
p7DF <- within(p7DF, rEta[nEta == 1] <- 'None')

# Path 8 #################################################################

#list files in a folder
p8 <- list.files(path = paste0('./path 8/',K,'/recovery/'),
                 pattern = ".dat",
                 full.names = TRUE)

#row bind files in same folder (stack data sets)
p8DF <- do.call("rbind", lapply(p8, read.table, header = TRUE))
p8DF <- p8DF[order(p8DF$nEta, p8DF$pEta, p8DF$rEta, p8DF$l, p8DF$n),] #put df in ascending order of factors
row.names(p8DF) <- 1:nrow(p8DF) #put row numbers back in ascending order

#create variables for analysis choices
p8DF$efa <- 1
p8DF$rot <- 2
p8DF$length <- 2

#move mean25 and max25 to end of data frame
p8DF <- p8DF %>% relocate(mean25, .after = last_col()) %>% relocate(max25, .after = last_col())

p8DF$d <- rep(minrep:maxrep, times = 168)

#change rEta for one factor to NA
p8DF <- within(p8DF, rEta[nEta == 1] <- 'None')

# Path 9 #################################################################

#list files in a folder
p9 <- list.files(path = paste0('./path 9/',K,'/recovery/'),
                 pattern = ".dat",
                 full.names = TRUE)

#row bind files in same folder (stack data sets)
p9DF <- do.call("rbind", lapply(p9, read.table, header = TRUE))
p9DF <- p9DF[order(p9DF$nEta, p9DF$pEta, p9DF$rEta, p9DF$l, p9DF$n),] #put df in ascending order of factors
row.names(p9DF) <- 1:nrow(p9DF) #put row numbers back in ascending order

#create variables for analysis choices
p9DF$efa <- 1
p9DF$rot <- 3
p9DF$length <- 2

#move mean25 and max25 to end of data frame
p9DF <- p9DF %>% relocate(mean25, .after = last_col()) %>% relocate(max25, .after = last_col())

p9DF$d <- rep(minrep:maxrep, times = 168)

#change rEta for one factor to NA
p9DF <- within(p9DF, rEta[nEta == 1] <- 'None')

# Path 10 ################################################################

#list files in a folder
p10 <- list.files(path = paste0('./path 10/',K,'/recovery/'),
                 pattern = ".dat",
                 full.names = TRUE)

#row bind files in same folder (stack data sets)
p10DF <- do.call("rbind", lapply(p10, read.table, header = TRUE))
p10DF <- p10DF[order(p10DF$nEta, p10DF$pEta, p10DF$rEta, p10DF$l, p10DF$n),] #put df in ascending order of factors
row.names(p10DF) <- 1:nrow(p10DF) #put row numbers back in ascending order

#create variables for analysis choices
p10DF$efa <- 2
p10DF$rot <- 1
p10DF$length <- 2

#move mean25 and max25 to end of data frame
p10DF <- p10DF %>% relocate(mean25, .after = last_col()) %>% relocate(max25, .after = last_col())

p10DF$d <- rep(minrep:maxrep, times = 168)

#change rEta for one factor to NA
p10DF <- within(p10DF, rEta[nEta == 1] <- 'None')

# Path 11 ################################################################

#list files in a folder
p11 <- list.files(path = paste0('./path 11/',K,'/recovery/'),
                 pattern = ".dat",
                 full.names = TRUE)

#row bind files in same folder (stack data sets)
p11DF <- do.call("rbind", lapply(p11, read.table, header = TRUE))
p11DF <- p11DF[order(p11DF$nEta, p11DF$pEta, p11DF$rEta, p11DF$l, p11DF$n),] #put df in ascending order of factors
row.names(p11DF) <- 1:nrow(p11DF) #put row numbers back in ascending order

#create variables for analysis choices
p11DF$efa <- 2
p11DF$rot <- 2
p11DF$length <- 2

#move mean25 and max25 to end of data frame
p11DF <- p11DF %>% relocate(mean25, .after = last_col()) %>% relocate(max25, .after = last_col())

p11DF$d <- rep(minrep:maxrep, times = 168)

#change rEta for one factor to NA
p11DF <- within(p11DF, rEta[nEta == 1] <- 'None')

# Path 12 ################################################################

#list files in a folder
p12 <- list.files(path = paste0('./path 12/',K,'/recovery/'),
                 pattern = ".dat",
                 full.names = TRUE)

#row bind files in same folder (stack data sets)
p12DF <- do.call("rbind", lapply(p12, read.table, header = TRUE))
p12DF <- p12DF[order(p12DF$nEta, p12DF$pEta, p12DF$rEta, p12DF$l, p12DF$n),] #put df in ascending order of factors
row.names(p12DF) <- 1:nrow(p12DF) #put row numbers back in ascending order

#create variables for analysis choices
p12DF$efa <- 2
p12DF$rot <- 3
p12DF$length <- 2

#move mean25 and max25 to end of data frame
p12DF <- p12DF %>% relocate(mean25, .after = last_col()) %>% relocate(max25, .after = last_col())

p12DF$d <- rep(minrep:maxrep, times = 168)

#change rEta for one factor to NA
p12DF <- within(p12DF, rEta[nEta == 1] <- 'None')

# Merge All Path Data Sets ###############################################

allDF <- rbind(p1DF,p2DF,p3DF,p4DF,p5DF,p6DF,p7DF,p8DF,p9DF,p10DF,p11DF,p12DF)

# Subset Data Where estK is NA ###############################################

naDF <- filter(allDF, is.na(allDF$estK) | is.na(allDF$targK))

write.table(naDF, paste0('./na/naData',K,'.dat'),
            row.names = FALSE,
            col.names = TRUE,
            sep = "\t")

nNA <- naDF %>%
            group_by(n) %>%
            summarize(Freq=n())

write.table(nNA, paste0('./na/naSummary',K,'.dat'),
            row.names = FALSE,
            col.names = TRUE,
            sep = "\t")

nEtaNA <- naDF %>%
               group_by(nEta) %>%
               summarize(Freq=n())

cat('\n', 
    file = paste0('./na/naSummary',K,'.dat'), 
    append = TRUE)

write.table(nEtaNA, 
            file = paste0('./na/naSummary',K,'.dat'), 
            append = TRUE, 
            row.names = FALSE, 
            col.names = TRUE,
            sep = "\t")

pEtaNA <- naDF %>%
               group_by(pEta) %>%
               summarize(Freq=n())

cat('\n', 
    file = paste0('./na/naSummary',K,'.dat'), 
    append = TRUE)

write.table(pEtaNA, 
            file = paste0('./na/naSummary',K,'.dat'), 
            append = TRUE, 
            row.names = FALSE, 
            col.names = TRUE,
            sep = "\t")

rEtaNA <- naDF %>%
               group_by(rEta) %>%
               summarize(Freq=n())

cat('\n', 
    file = paste0('./na/naSummary',K,'.dat'), 
    append = TRUE)

write.table(rEtaNA, 
            file = paste0('./na/naSummary',K,'.dat'), 
            append = TRUE, 
            row.names = FALSE, 
            col.names = TRUE,
            sep = "\t")

lNA <- naDF %>%
            group_by(l) %>%
            summarize(Freq=n())

cat('\n', 
    file = paste0('./na/naSummary',K,'.dat'), 
    append = TRUE)

write.table(lNA, 
            file = paste0('./na/naSummary',K,'.dat'), 
            append = TRUE, 
            row.names = FALSE, 
            col.names = TRUE,
            sep = "\t")

efaNA <- naDF %>%
              group_by(efa) %>%
              summarize(Freq=n())

cat('\n', 
    file = paste0('./na/naSummary',K,'.dat'), 
    append = TRUE)

write.table(efaNA, 
            file = paste0('./na/naSummary',K,'.dat'), 
            append = TRUE, 
            row.names = FALSE, 
            col.names = TRUE,
            sep = "\t")

rotNA <- naDF %>%
              group_by(rot) %>%
              summarize(Freq=n())

cat('\n', 
    file = paste0('./na/naSummary',K,'.dat'), 
    append = TRUE)

write.table(rotNA, 
            file = paste0('./na/naSummary',K,'.dat'), 
            append = TRUE, 
            row.names = FALSE, 
            col.names = TRUE,
            sep = "\t")

lengthNA <- naDF %>%
                 group_by(length) %>%
                 summarize(Freq=n())

cat('\n', 
    file = paste0('./na/naSummary',K,'.dat'), 
    append = TRUE)

write.table(lengthNA, 
            file = paste0('./na/naSummary',K,'.dat'), 
            append = TRUE, 
            row.names = FALSE, 
            col.names = TRUE,
            sep = "\t")

# Write Out ANOVA Data Set ###############################################

write.table(allDF, paste0('./results/anovaData',K,'.dat'),
            row.names = FALSE,
            col.names = TRUE,
            sep = "\t")
