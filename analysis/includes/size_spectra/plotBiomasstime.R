# Plot the temporal evolution of biomasses in the system (check for species dying)

plotBiomasstime <- function(param,S){
  
  time <- S$t
  Biomass <- S$Biomass
  
  dfnew <- data.frame(time = time, Biomass = Biomass[,1])
  p <- ggplot(data = dfnew, aes(x = time, y = Biomass))+geom_line()+
    theme_bw()+scale_y_log10('Biomass')+scale_x_continuous('time (years)')+theme(text = element_text(size=9))
  
  if (param$nSpecies> 1){
    for (i in 2:param$nSpecies){
    p <- p + geom_line(data = data.frame(x=time,y = Biomass[,i]),aes(x=x,y=y))
  }
}
p
return(p)
}

