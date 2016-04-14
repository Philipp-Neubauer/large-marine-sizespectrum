
plotVonBGrowth  <- function (param,S,data){
  require(ggplot2)
  require(deSolve)
  require(signal)
  require(rfishbase)
  
  t <- param$tEnd
  w <- S$w
  
  iPlot <- round(S$nSave*t/param$tEnd);
  idxRange <- floor(iPlot/2):iPlot;
  
  tRange <- seq(0,50,by= 0.02)
  hbar <- param$alpha*param$f0*param$h-param$ks
  
  
  if (length(hbar) == 1){
    hbar = hbar * ones(1,param$nSpecies);
  }
  state <- param$w0
  # Ode solver 
  deriv <-function(t, state, parameters) {
    with(as.list(c(state, parameters)),{
      # rate of change
      dX <- interp1(w,g,state)
      # return the rate of change
      list(dX)
    }) # end with(as.list ...
  }
  
  
  
  
  wCalc <- matrix(NA,param$nSpecies,length(tRange))
  for (i in 1:param$nSpecies){
    g = colMeans(S$g[idxRange,i,])
    
    SolveG = ode(y = param$w0,times = tRange, func = deriv, parms = c(g,w))
    wCalc[i,] <- SolveG[,2]
  }
  
  
  # Solve the von Bertalanffy growth 
  
  dfVonB <- data[which(is.na(data$k) == 0),]
  dfVonB$Linf <- (dfVonB$wInf/0.01)^(1/3)
  
  dfSS <- wCalc[which(is.na(data$k) == 0),]
  

  # Adjust time for birth at w = 0:
  t = tRange - 1/(hbar[1]*(param$n-1)) * param$w0^(1-param$n)
  tNew <- 1:50
  VonB <- matrix(0,length(dfVonB$k),length(tNew))
  for (i in 1:length(dfVonB$k)){
    VonB[i,] <- dfVonB$Linf[i]*(1-exp(-dfVonB$k[i]*(tNew-dfVonB$t0[i])))  
  }
  
  deNom <- matrix(0,param$nSpecies,length(tNew))
  tScaled <- matrix(0,param$nSpecies,length(t))
  for (i in 1:param$nSpecies){
  deNom[i,] <- tNew/as.numeric((param$h[i]*(param$alphaMature*param$wInf[i])^(1/4)))
  tScaled[i,] <- t/as.numeric((param$h[i]*(param$alphaMature*param$wInf[i])^(1/4)))
  }
  
  par(ask=TRUE)
  p <- list()
  for (i in 1:length(dfVonB$k)){
    
  print(ggplot(data =data.frame(growth = VonB[i,],time = tNew), aes(x = time, y = growth))+
      geom_line(linetype = 2, size = 2)+
          theme_classic()+geom_line(data = data.frame(growth = (dfSS[i,]/0.01)^(1/3), time = t), size = 2)+
          scale_y_continuous(name = "size (cm)")+
        scale_x_continuous(name = 'time (years)')+
      annotate("text",x=30,y=2,label = as.character(dfVonB$cnames[i]), size = 15)+theme(text = element_text(size=20)))
  p[[i]] <-ggplot(data =data.frame(growth = VonB[i,],time = tNew), aes(x = time, y = growth))+
    geom_line(linetype = 2, size = 2)+
    theme_classic()+geom_line(data = data.frame(growth = (dfSS[i,]/0.01)^(1/3), time = t), size = 2)+
    scale_y_continuous(name = "size (cm)")+
    scale_x_continuous(name = 'time (years)')+
    annotate("text",x=30,y=2,label = as.character(dfVonB$cnames[i]), size = 15)+theme(text = element_text(size=20))
#   plot(t,VonB[i,], type = 'l')
#   lines(t,(dfSS[i,]/0.01)^(1/3))
#     
    
  }
  par(ask = FALSE)
  
  #width = seq(0.5,2,length.out = param$nSpecies)
    
  return(p)
}