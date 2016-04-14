# Average 
AverageList <- function(dfTemp,year){
# Sum or average the ones with the same names

cnames <- unique(names(dfTemp))
totnames <- names(dfTemp)
is.nan.data.frame <- function(x)
  do.call(cbind, lapply(x, is.nan))



dfTempNew <- list()

for (i in 1:length(cnames)){
  tmpidx <- which(cnames[i] == totnames)
  tmplist<- dfTemp[tmpidx]
  
  if (length(tmplist[[1]]) > 1){
  nYears <- length(year)
  
  dftmp <- matrix(NA,nYears,7)
   
  F0tmp <- matrix(NA,nYears,length(tmpidx))
  Catchtmp <- matrix(NA,nYears,length(tmpidx))
  Landingstmp <- matrix(NA,nYears,length(tmpidx))
  Biomasstmp <- matrix(NA,nYears,length(tmpidx))
  SSBtmp <- matrix(NA,nYears,length(tmpidx))
  Rtmp <- matrix(NA,nYears,length(tmpidx))
  
  for (j in 1:length(tmpidx)){
    
    if (length(tmplist[[j]]$year) < nYears){
      dfNew <- data.frame(F0 = NA, Catch=NA,Biomass=NA,Landings = NA,SSB = NA, year = year, R = NA)
      idxEnd <- which(tmplist[[j]]$year[length(tmplist[[j]]$year)] == year)
      idxStart <- which(tmplist[[j]]$year[1] == year)
      dfNew$F0[idxStart:idxEnd] <- tmplist[[j]]$F0
      dfNew$Catch[idxStart:idxEnd] <- tmplist[[j]]$Catch
      dfNew$Landings[idxStart:idxEnd] <- tmplist[[j]]$Landings
      dfNew$Biomass[idxStart:idxEnd] <- tmplist[[j]]$Biomass
      dfNew$SSB[idxStart:idxEnd] <- tmplist[[j]]$SSB
      dfNew$R[idxStart:idxEnd] <- tmplist[[j]]$R
      
      # Replace just the biomass and SSB
      #dfNew$Biomass[1:(idxStart)] <- tmplist[[j]]$Biomass[idxStart]
      #dfNew$SSB[1:(idxStart)] <- tmplist[[j]]$SSB[idxStart]
      
      tmplist[[j]] <- dfNew
    }
    
    
    F0tmp[,j] <- tmplist[[j]]$F0
    Catchtmp[,j] <- tmplist[[j]]$Catch
    Landingstmp[,j] <- tmplist[[j]]$Landings
    Biomasstmp[,j] <- tmplist[[j]]$Biomass
    SSBtmp[,j] <- tmplist[[j]]$SSB
    Rtmp[,j] <- tmplist[[j]]$R
  }
dftmp[,1] <- rowMeans(F0tmp, na.rm=T)
dftmp[,2] <- rowSums(Catchtmp, na.rm=T)
dftmp[,3] <- rowSums(Landingstmp, na.rm=T)
dftmp[,4] <- rowSums(Biomasstmp, na.rm=T)
dftmp[,5] <- rowSums(SSBtmp, na.rm=T)
dftmp[,6] <- tmplist[[1]]$year
dftmp[,7] <- rowSums(Rtmp, na.rm = T)

dftmp <- as.data.frame(dftmp)
dftmp[is.nan(dftmp)] <- NA
dftmp[dftmp == 0] <- NA

names(dftmp) <- names(dfTemp[[i]])

dfTempNew[[i]] <- dftmp
}else{
  dfTempNew[[i]] = NA
}

}
names(dfTempNew) <- cnames
return(dfTempNew)
}

