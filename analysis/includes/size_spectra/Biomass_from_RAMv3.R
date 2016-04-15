Biomass_from_RAMv3 = function(LMEnumber, dir){
  # Takes input LMEnumber and dir, where dir is the working directory and LMENumber is the number of the large marine ecosystem in the ram database
  
  # Get the biomass distribution currently (or from some particular year, maybe mean?) 
  # Let's make a script where we test the "Northeast U.S. Continental Shelf" - that way I can compare with the LeMans model 
  
  
  ## Indirect effects of fishing marine ecosystems 
  require('RODBC')
  #require(rfishbase)
  #require(RCurl)
  #require(XML)
  #require(grid)
  #require(gridExtra)
  #require(ggplot2)
  #require(rgdal)
  
  # load The ram database into R 

  #db <- "C:/Users/s102045/Dropbox/Ocean life/Size based models in R/RLSADBv3.accdb"
  #db_old <- "C:/Users/s102045/Downloads/RAM-Legacy-v1.0.accdb" # Old RAM database has LME information 
  
  if(.Platform$OS.type == 'unix'){
    require(Hmisc)
    #need mdb - install instructions under linux are here
    old_db <- Hmisc::mdb.get(paste0(db,"RAM-Legacy-v1.0.accdb"))
    attach(old_db)
    full_db <- Hmisc::mdb.get(paste0(db,"RLSADBv3.accdb"))
    attach(full_db)
    
    ts_values <- timeseries_values_views
    ts_units <- timeseries_units_views
  } else {
    db <- paste(dir,"RLSADBv3.accdb", sep="/")
    db_old <- paste(dir,"RAM-Legacy-v1.0.accdb", sep = "/") # Old RAM database has LME information 
    
    # load the new RAM legacy database 
    
    con2 <- odbcConnectAccess(db) #  EUREKA! 
    con_old <- odbcConnectAccess2007(db_old)
    
    # Get the table headers 
    info <- sqlTables(con2) # Find tables in the ramdatabase 
    info <- info$TABLE_NAME
    
    
    
    assessment <- sqlFetch(con2,'assessment')
    lmes <- sqlFetch(con_old,'lmerefs')
    lmetostocks <- sqlFetch(con_old,'lmetostocks')
    area <- sqlFetch(con2,'area')
    # To get the values out of each table use the "sqlFetch" command
    
    
    ts_values <- sqlFetch(con2,'timeseries_values_views')
   
    bioparams <- sqlFetch(con2,'bioparams')
    ts_units <- sqlFetch(con2,'timeseries_units_views')
    stock <- sqlFetch(con2,'stock')
    rp <- sqlFetch(con2,'bioparams_values_views')
    tsmetrics <- sqlFetch(con2,'tsmetrics')
    #fishnames <- unique(fishnames) # Only find the unique names 
    timeseries <- sqlFetch(con2,'timeseries')
    
    
    
    lmenum <- LMEnumber # LME number, used to find the data from the specific lme 
  }
    df_ref <- read.table(paste0(db,'Reference_points_RAMv3.csv'))
    
    # Just remove the secondary things 
    lmetostocks <- lmetostocks[lmetostocks$STOCKTOLMERELATION == "primary",]
    lmetostocks$STOCKID <- as.character(lmetostocks$STOCKID)
    lmetostocks$STOCKID[113] <- "WHITNS-VIId" # disagreement between RAM versions
    stock$commonname <- as.character(stock$commonname)
  stock$commonname[1] <- "Acadian Redfish" # Fix with capital letter
  # fix difference in capital letters
  stock$commonname[which(stock$commonname == "Yellowtail flounder")] <- "Yellowtail Flounder" 
  # Move down to unique stock ID's 
  
  
  LMEidx <- which(lmetostocks$LME == lmenum)
  fishnamesLME <- lmetostocks$STOCKID[LMEidx]
  
  if (length(fishnamesLME) == 0){
    stop('No species from this LME recorded')
  }
  
  
  AreaCode <- matrix(NA,length(fishnamesLME))
  for (i in 1:length(fishnamesLME)){
    idxtmp <- stock[which(as.character(fishnamesLME[i]) == as.character(stock$stockid)),]
    if (length(idxtmp$areaid) > 0){
    AreaCode[i] <- unique(as.character(idxtmp$areaid))
        if (length(unique(as.character(idxtmp$areaid))) > 1){
            print('More than one entry for stock')
            print(i)
                    }
    }
  }
  
  # Find remaining species in the specified areas 
  Areas <- unique(AreaCode)
  
 # ts_values[19522:19553,]$TC <- NA # Catch is in NUMBERS
#  ts_values[19522:19553,]$SSB <- ts_values[19522:19553,]$SSB*2 # Only females in asssessment
  
  
  SpeciesUnique <- list()
  for (i in 1:length(Areas)){
    idxtmp <- stock[which(as.character(stock$areaid) == Areas[i]),]
    SpeciesUnique[[i]] <- idxtmp$stockid
    
    
  }
  
  SpeciesUnique <- unique(unlist(SpeciesUnique,recursive = FALSE))
  fishnamesLME <- SpeciesUnique
  
  FishIdx <- getFishIndices(ts_values,fishnamesLME,singleR = FALSE)
  
  FishIdx <- FishIdx[[1]]
  
  
  
  df <- data.frame(F0 = matrix(NA,length(fishnamesLME)),
                   Catch = matrix(NA,length(fishnamesLME)),
                   Landings =matrix(NA,length(fishnamesLME)),
                   Biomass = matrix(NA,length(fishnamesLME)),
                   SSB = matrix(NA,length(fishnamesLME)),
                   M = matrix(NA,length(fishnamesLME)),
                   Fmsy = matrix(NA,length(fishnamesLME)),
                   wInf = matrix(NA,length(fishnamesLME)),
                   k = matrix(NA,length(fishnamesLME)),
                   cnames = matrix(NA,length(fishnamesLME)),
                   t0 = matrix(NA,length(fishnamesLME)),
                   lMat = matrix(NA,length(fishnamesLME)),
                   aMat = matrix(NA,length(fishnamesLME)),
                   fnames = fishnamesLME)
  # Throw info from reference points into the new df 
  
  OrgNames <- list()
  
  for (i in 1:length(fishnamesLME)){
    #  idxtmp <- grep(fishnames[i],as.character(ts_values[,1]))
    #  species[i] <-list(ts_values[idxtmp,])
    tmpidx <- which(stock$stockid == as.character(fishnamesLME[i]))
    if (length(tmpidx) > 0){
      OrgNames[i] <- as.character(stock$commonname[tmpidx])
      }else{
      OrgNames[i] <- NA}
  }
  #RPd
  
  Units <- data.frame(matrix(NA,length(fishnamesLME),8)) # Why 8??
  names(Units) <- names(ts_units[4:11])
  
  df$cnames <- as.character(OrgNames)
  # Get new info 
  for (i in 1:length(fishnamesLME)){
    
    tempdf <- df_ref[which(as.character(df_ref$fnames) == as.character(fishnamesLME[i])),]
    
    if (dim(tempdf)[1] > 0){
      df$M[i] <- tempdf$M
      df$Fmsy[i] <- tempdf$Fmsy
      df$wInf[i] <- tempdf$wInf
      df$k[i] <- tempdf$k
      df$t0[i] <- tempdf$t0
      df$aMat[i] <- tempdf$aMat
     # df$wMat[i] <- tempdf$wMat
      #df$cnames[i] <- as.character(tempdf$CNames)
    }
    
    # Now get the average B, catch, F and SSB
    if (is.na(FishIdx[i,1])){
      df$F0[i] <- NA
    }else{
      tmpmat <- ts_values[FishIdx[i,1]:FishIdx[i,2],]
      
      # Let's try 1992 to 2002 
      idxMin <- which(tmpmat$year == 1992)
      if (length(idxMin)==0){
        idxMin <- which.min(tmpmat$year)
      }
      
      
      idxMax <- which(tmpmat$year == 2002)
      if (length(idxMax) == 0){  
        idxMax <- which.max(tmpmat$year[tmpmat$TC > 0])-1 # Usually some unavailable information from most recent year
      }
      
      
      df$F0[i] = mean(tmpmat$F[idxMin:idxMax], na.rm = T)
      df$Catch[i] = mean(tmpmat$TC[idxMin:idxMax], na.rm = T)
      df$Landings[i] <- mean(tmpmat$TL[idxMin:idxMax], na.rm = T)
      df$Biomass[i] = mean(tmpmat$TB[idxMin:idxMax], na.rm = T)
      df$SSB[i] = mean(tmpmat$SSB[idxMin:idxMax], na.rm = T)
      
      # Save the units to make sure there's consistency 
      nm <- tmpmat[1,1]
      un <- as.matrix(ts_units[as.character(ts_units$assessid) == as.character(nm),4:11])
      Units[i,] <- un
      
    }
  }  
  odbcCloseAll() # close off connections in the end 
  
  return(list(df,Units))
}
