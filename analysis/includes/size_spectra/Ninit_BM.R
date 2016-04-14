Ninit_BM <- function(param,S,obsBio,state){
  
# Calculate the start biomass with number distribution from SF 
tEnd <- param$tEnd/param$dt

Ninit <- S$N[tEnd, , ]

BMinit <- S$Biomass[tEnd,]
SSBinit <- as.numeric(calcSSB(param,S,tEnd))

SSBio <- SSBinit
SSBio[state$totFlag == 1] <- BMinit[state$totFlag == 1]

initFactor <- obsBio/SSBio

for (i in 1:param$nSpecies){
S$N[tEnd,i,] <- S$N[tEnd,i,]*initFactor[i]
S$Biomass[tEnd,i] <- sum(S$N[tEnd,i,]*S$w*S$dw)
}

return(S)
}