getTemporalInfo = function(LMEnumber,tStart,tEnd){
  
  
  # Get the biomass distribution currently (or from some particular year, maybe mean?) 
  # Let's make a script where we test the "Northeast U.S. Continental Shelf" - that way I can compare with the LeMans model 
  
  
  ## Indirect effects of fishing marine ecosystems 
  library('RODBC')
  library(rfishbase)
  library(RCurl)
  library(XML)
  library(grid)
  library(gridExtra)
  library(ggplot2)
  library(rgdal)
  
  
  
  source("getFishIndices.R")
  # load The ram database into R 
  
  #db <- "C:/Users/s102045/Dropbox/Ocean life/Size based models in R/RLSADBv3.accdb"
  #db <- "C:/Users/Nis/Documents/My Dropbox/Ocean life/Size based models in R/RLSADBv3.accdb"
  db <- paste(getwd(),"RLSADBv3.accdb", sep = '/')
  
  db_old <- "C:/Users/s102045/Downloads/RAM-Legacy-v1.0.accdb" # Old RAM database has LME information 
  db_old <- paste(getwd(),"RAM-Legacy-v1.0.accdb", sep = '/')
  
  #db_old <- "C:/Users/Nis/Documents/My Dropbox/RAM database/RAM-Legacy-v1.0.accdb"
  # load the new RAM legacy database 
  
  con2 <- odbcConnectAccess2007(db) #  EUREKA! 
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
  
  # Load cusk total catches 
  cusk <- read.table('cuskLandings.txt')
  cusk <- cusk[,-c(2:17)]
  
  ts_values$TC[6470:6506] <- cusk[,2]
  ts_values[19522:19553,]$TC <- ts_values[19522:19553,]$TC/1000 # Catch is in NUMBERS
  ts_values[19522:19553,]$SSB <- ts_values[19522:19553,]$SSB*2 # Only females in asssessment
  ts_values$F[5010] <- 0.97
  
  bioparams <- sqlFetch(con2,'bioparams')
  ts_units <- sqlFetch(con2,'timeseries_units_views')
  stock <- sqlFetch(con2,'stock')
  rp <- sqlFetch(con2,'bioparams_values_views')
  tsmetrics <- sqlFetch(con2,'tsmetrics')
  #fishnames <- unique(fishnames) # Only find the unique names 
  timeseries <- sqlFetch(con2,'timeseries')
  
  
  
  lmenum <- LMEnumber # LME number, used to find the data from the specific lme 
  
  df_ref <- read.table('Reference_points_RAMv3.csv')
  
  # Just remove the secondary things 
  lmetostocks <- lmetostocks[lmetostocks$STOCKTOLMERELATION == "primary",]
  lmetostocks$STOCKID <- as.character(lmetostocks$STOCKID)
  lmetostocks$STOCKID[113] <- "WHITNS-VIId" # disagreement between RAM versions
  
  stock$commonname <- as.character(stock$commonname)
  stock$commonname[1] <- "Acadian Redfish" # Fix with capital letter
  # fix difference in capital letters
  stock$commonname[which(stock$commonname == "Yellowtail flounder")] <- "Yellowtail Flounder" 
  
  
  # Move down to unique stock ID's 
  
  
  LMEidx <- which(lmetostocks$LME_NUMBER == lmenum)
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
  
  SpeciesUnique <- list()
  for (i in 1:length(Areas)){
    idxtmp <- stock[which(as.character(stock$areaid) == Areas[i]),]
    SpeciesUnique[[i]] <- idxtmp$stockid
    
    
  }
  SpeciesUnique <- unlist(SpeciesUnique,recursive = FALSE)
  fishnamesLME <- SpeciesUnique
  
  FishIdx <- getFishIndices(ts_values,fishnamesLME,singleR = FALSE)
  
  FishIdx <- FishIdx[[1]]
  
  
  # Throw info from reference points into the new df 
  
  OrgNames <- list()
  
  for (i in 1:length(fishnamesLME)){
    #  idxtmp <- grep(fishnames[i],as.character(ts_values[,1]))
    #  species[i] <-list(ts_values[idxtmp,])
    tmpidx <- which(stock$stockid == as.character(fishnamesLME[i]))
    if (length(tmpidx) > 0){
      OrgNames[i] <- as.character(stock$commonname[tmpidx])}
    else{
      OrgNames[i] <- NA}
  }
  #RPd
  
  Units <- data.frame(matrix(NA,length(fishnamesLME),8)) # Why 8??
  names(Units) <- names(ts_units[4:11])
  
  df <- list()
  # Get new info 
  for (i in 1:length(fishnamesLME)){
    
   
      if (is.na(FishIdx[i,1]) == 0){
      tmpmat <- ts_values[FishIdx[i,1]:FishIdx[i,2],]
      
      # Let's try 1992 to 2002 
      idxMin <- which(tmpmat$year == tStart)
      if (length(idxMin)==0){
        idxMin <- which.min(tmpmat$year)
      }
      
      
      idxMax <- which(tmpmat$year == tEnd)
      if (length(idxMax) == 0){  
        idxMax <- which.max(tmpmat$year[tmpmat$TC > 0])-1 # Usually some unavailable information from most recent year
      }
      if (length(idxMax)==0){
        idxMax <- which.max(tmpmat$year)
      }
      
      
      df[[i]]= data.frame(F0 = tmpmat$F[idxMin:idxMax], Catch = tmpmat$TC[idxMin:idxMax],
        Landings = tmpmat$TL[idxMin:idxMax],Biomass = tmpmat$TB[idxMin:idxMax],
        SSB = tmpmat$SSB[idxMin:idxMax],year = tmpmat$year[idxMin:idxMax], R = tmpmat$R[idxMin:idxMax])
      
      # Save the units to make sure there's consistency 
      }
      else{df[[i]] = NA}
  }  
  odbcCloseAll() # close off connections in the end 
  names(df) <- fishnamesLME
  return(df)
}
