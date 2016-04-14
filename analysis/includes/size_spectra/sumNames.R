# Function to sum up all the same names 

sumNames <- function(state,cnames){
  
  state_new <- as.data.frame(matrix(NA, length(cnames), length(names(state))))
  names(state_new) <- names(state)
  state_new$cnames <- cnames
  
  is.nan.data.frame <- function(x){
    do.call(cbind, lapply(x, is.nan))}
  
  state[is.nan.data.frame(state)] <- NA
  
  
  for (i in 1:length(cnames)){
    tmpidx <- which(cnames[i] == state$cnames)
    state_new$F0[i] <- mean(state$F0[tmpidx], na.rm = T)
    state_new$Catch[i] <- sum(state$Catch[tmpidx], na.rm = T)
    state_new$Landings[i] <- sum(state$Landings[tmpidx], na.rm = T)
    state_new$Biomass[i] <- sum(state$Biomass[tmpidx], na.rm = T)
    state_new$SSB[i] <- sum(state$SSB[tmpidx], na.rm = T)
    state_new$M[i] <-  mean(state$M[tmpidx], na.rm = T)
    state_new$Fmsy[i] <-  mean(state$Fmsy[tmpidx], na.rm = T)
    state_new$k[i] <-   mean(state$k[tmpidx], na.rm = T)
    state_new$t0[i] <- mean(state$t0[tmpidx], na.rm = T)
    state_new$wInf[i] <- mean(state$wInf[tmpidx], na.rm = T)
  }
  
  state_new$Biomass[state_new$Biomass == 0] <- NA
  state_new$Catch[state_new$Catch == 0] <- NA
  state_new$Landings[state_new$Landings == 0] <- NA
  state_new$SSB[state_new$SSB == 0] <- NA
  
 return(state_new) 
  
}
