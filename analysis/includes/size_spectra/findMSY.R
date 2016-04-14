# Find MSY for each species 

findMSY <- function(param,S0,Species,state){
  
if (is.na(state) == 1){
  state <- data.frame(F0 = as.numeric(0.4*matrix(1,param$nSpecies))) 
}

param$tEnd <- 40

F0 <- seq(0.01,2,length.out = 50) # start with F = 0.1

Yield <- matrix(0,length(F0))
for (i in 1:length(F0)){
  
  if (i == 1){
    S0 <- S0
  }else{
    S0 <- SF
  }
  param$F0 <- state$F0
  param$F0[Species] <- F0[i]

  SF <- IterateSpectrum(param,S0)
  Yieldtmp <- YieldCalc(param,SF)
  Yield[i] <- Yieldtmp[Species]
  
  if (i > 1){
    if (Yield[i] < Yield[i-1]){
      break
  }}
  print(F0[i])
}

Fmsy <- F0[which.max(Yield)]
  return(Fmsy)
}