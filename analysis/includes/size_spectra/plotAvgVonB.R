
plotAvgVonB  <- function (param,S,data){
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
  
  dfVonB$t0[is.na(dfVonB$t0 )] <- 0
  # Adjust time for birth at w = 0:
  t = tRange - 1/(hbar[1]*(param$n-1)) * param$w0^(1-param$n)
  tNew <- 1:50
  VonB <- matrix(0,length(dfVonB$k),length(tNew))
  for (i in 1:length(dfVonB$k)){
    VonB[i,] <- dfVonB$Linf[i]*(1-exp(-dfVonB$k[i]*(tNew-dfVonB$t0[i])))  
  }
  
  VonBweight <- 0.01*VonB^3
  df <- data.frame(t(VonBweight),tNew)

  names(df) <- data$cnames
  
#   
#   V <- ggplot(data=df,aes(x=tNew/(param$h[1]*param$wInf[1]^(1-param$n)),
#                           y=df[,1]/data$wInf[1]))+geom_line()+theme_classic()
#   
#   for (i in 2:10){
#   V <- V + geom_line(data=data.frame(x=tNew/(param$h[i]*param$wInf[i]^(1-param$n)),y=df[,i]/data$wInf[i]),aes(y=y) )  
#     
#   }
#   V
  # Plot the averaged von B from the system at hand
  tmat <-t(param$h*(param$alphaMature*param$wInf)^(1-param$n))
  df <-data.frame(t=t/tmat[1],weight = wCalc[1,]/param$wInf[1])
  
  p <- ggplot(data =df,aes(x = t,y = weight))+geom_line()+theme_classic()+
    scale_x_continuous(expand=c(0,0), name = 't/(h*m^(1/4))', limit = c(0,2))+
    scale_y_continuous(expand = c(0,0), name = expression(paste("w/",'W'[infinity])))+
                         theme(text = element_text(size=9))
    
  for (i in 2:param$nSpecies){
    df1 <-data.frame(t=t/tmat[i],weight = wCalc[i,]/param$wInf[i])
    p <- p+geom_line(data=df1,aes(x=t,y=weight))
    }
  p
  wmat <- param$h*(param$alphaMature*dfVonB$wInf)^(1-param$n)
  
  for (i in 1:dim(VonB)[1]){
    df2 <- data.frame(t=tNew/wmat[i],weight = 0.01*VonB[i,]^3/(0.01*dfVonB$Linf[i]^3))
    p <- p+geom_line(data=df2,aes(x=t,y=weight), color = 'gray', linetype = 1)
  }

p
return(p)
}