RmaxCalc_oneReturn <-  function(state_new,param,S0,Rmax){
  
  source('load_files.R')
  #
  library(stats4)
  require('optimx')
  
  SScalc <- function(data,par){ 
    
    
    state_new <- data$state_new
    S <- data$S
    param <- data$param
    
    
    param$Rmax <- exp(par)
    
    S <- IterateSpectrum(param,S)
    
    tEnd <- param$tEnd/param$dt
    SSB <- calcSSB(param,S,tEnd) # Calculate SSB
    Biomass <- S$Biomass[tEnd,] # Calculate Biomass
    
    SSBio <- matrix(NA,param$nSpecies) # Set the right quantities up
    SSBio[state_new$totFlag == 0] <- SSB[state_new$totFlag == 0]
    SSBio[is.na(SSBio)] <- Biomass[which(is.na(SSBio) == 1)]
    
    # Sum of squares
    SSQ <- sum((log10(state_new$SSBio)-log10(SSBio))^2)
    print(SSQ)
    return(SSQ)
    
  }
  
  
  #MLE <- mle(SScalc, start=list(kappa = 10^(5.2)),fixed=list(state_new = state_new,S0=S,wInf = state_new$wInf))

#   Rmax[1] <- 1e5
#   Rmax[1] <- 150000
#   Rmax[2] <- 25000
#   
#   
  param$F0 <- state_new$F0
  param$tEnd <- 50
  param$fishing <- "Trawl" 
  options(digits = 3)
  param$Rmax <- as.numeric(Rmax)
  
  
  minKappa <- optimx(par = log(param$Rmax), fn = SScalc, 
                    data = list(state_new = state_new, S = S0, param = param),
                    control =list(maxit = 1000),method = "bobyqa")
  Rmaxfit <- minKappa
  
  return(Rmaxfit)  
}