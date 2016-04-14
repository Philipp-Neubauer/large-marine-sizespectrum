# Plot the temporal evolution of biomasses in the system (check for species dying)

plotRpR <- function(param,S){
  
  tEnd <- param$tEnd/param$dt
  wInf <- t(param$wInf)
  Rp <- S$Rp[tEnd,]
  R <- S$R[tEnd,]
  
  dfnew <- data.frame(wInf = wInf, Rp = Rp, R = R)
  p <- ggplot(data = dfnew, aes(x = wInf, y = Rp/R))+geom_line()+
    theme_classic()+scale_y_log10('Rp/R', breaks = c(1,5,10,100,1000))+
    scale_x_log10('Asymptotic weight (g)', breaks = 10^seq(log10(100),log10(1e6),length.out =5))+
    geom_point()+
    geom_hline(yintercept = 1)
  
  
  p
  return(p)
}
