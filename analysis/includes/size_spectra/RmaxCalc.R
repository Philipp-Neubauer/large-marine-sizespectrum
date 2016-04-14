RmaxCalc <- function(state_new,param,S0){
  
source('load_files.R')
    #
library(stats4)

    

SScalc <- function(data,par){ 
      
      
      state_new <- data$state_new
      idx <- data$idx
      S <- data$S
      param <- data$param
      
      
      param$Rmax[idx] <- exp(par)
      param$Rmax <- as.numeric(param$Rmax)
      
      S <- IterateSpectrum(param,S)
      
      tEnd <- param$tEnd/param$dt
      SSB <- calcSSB(param,S,tEnd) # Calculate SSB
      Biomass <- S$Biomass[tEnd,] # Calculate Biomass
      
      SSBio <- matrix(NA,param$nSpecies) # Set the right quantities up
      SSBio[state_new$totFlag == 0] <- SSB[state_new$totFlag == 0]
      SSBio[is.na(SSBio)] <- Biomass[which(is.na(SSBio) == 1)]
      
      # Sum of squares
      SSQ <- ((log10(state_new$SSBio[idx])-log10(SSBio[idx]))^2)
      print(SSQ)
      print(idx)
      return(list(SSQ))
      
    }
    
    
    #MLE <- mle(SScalc, start=list(kappa = 10^(5.2)),fixed=list(state_new = state_new,S0=S,wInf = state_new$wInf))



param$tEnd <- 60
param$fishing <- "Trawl" 
options(digits = 3)
param$wInf <- as.numeric(as.matrix(param$wInf))
Rmax <- matrix(NA,param$nSpecies)

  for(i in 1:param$nSpecies){
      i <- (param$nSpecies+1)-i    
      minKappa <- optim(par = log(param$Rmax[i]), fn = SScalc, 
                        data = list(state_new = state_new, idx = i, S = S0, param = param),
                      control =list(maxit = 100),method = 'L-BFGS-B',lower = 1e-05, upper = 1e20)
      Rmax[i] <- exp(minKappa$par)
      param$Rmax[i] <- Rmax[i]
      sprintf("%s %f is loading", "Species", as.integer(i))
 
  }
    return(Rmax)  
  }