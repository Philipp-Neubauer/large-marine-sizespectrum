plotTemporalBiomass <- function(df_timeCalc,df_timeCalcBM,dfTemp,years,parOld,BiomassStart,S92){
  source('multiplot.R')

  nSpecies <- parOld$nSpecies
  
  nm <- names(df_timeCalc)
  p <- list()
  
  # For publication show every second year 
#  xtick <- as.character(c(1980,'',1982,'',1984,'',1986,'',1988,'',1990,'',
#                         1992,'',1994,'',1996,'',1998,'',2000,
#                          '',2002,'',2004,'',2006,'',2008,'',2010))
  xtick <- years
  
 
#  YtiKbarents <- data.frame(Capelin = c(1e5,1e6),Redfish = c(1e4,1e5),
#                            Haddock = c(1e5,5e5),Pollock = c(1e5,1e6),
#                            halibut = c(1e4,1e5),cod = c(1e5,1e6))
#   
for (i in 1:nSpecies){
  
  yStart <- which(df_timeCalc$years ==min(dfTemp[[i]]$year,na.rm = TRUE))
  yEnd <- which(df_timeCalc$years ==max(dfTemp[[i]]$year,na.rm=TRUE))
  
  if (dfTemp[[i]]$flag[1] == 0){
  df <- data.frame(Bobs = NA,Bsim = NA, years = years)
  df$Bobs[yStart] <- BiomassStart[i]
  df$Bobs[(yStart+1):(yEnd+1)] <- dfTemp[[i]]$SSB
  
  df$Bsim[2:length(years)] = df_timeCalc[,i]
  df$Bsim[1] <- calcSSB(parOld,S92,length(S92$t))[i]
  
  
  }else{
  df <- data.frame(Bobs = NA,Bsim = NA, years = years)
  df$Bobs[yStart] <- BiomassStart[i]
  df$Bobs[(yStart+1):(yEnd+1)] <- dfTemp[[i]]$Biomass
  
  df$Bsim[2:length(years)] = df_timeCalcBM[,i]
  df$Bsim[1] <-S92$Biomass[length(S92$t),i]
  }
  y_limit_min <- min(c(df$Bsim,df$Bobs),na.rm=T)
  y_limit_min <- y_limit_min * 0.5
  
  y_limit_max <- max(c(df$Bsim,df$Bobs),na.rm=T)
  y_limit_max <- y_limit_max * 2
  
  
  # Calculate the correct ticks 
  # Get the correct tick 
  ytickmin <- signif(y_limit_min,1)
  ytickmax <- signif(y_limit_max,1)
  
  ytick <- c(ytickmin,ytickmax)
  ytick <- seq(log10(ytickmin),log10(ytickmax),length.out = 3)
  ytick[2] <- signif(ytick[2],1)
  ytick <- 10^ytick
  
  y_limit_min <- ytick[1]
  y_limit_max <- ytick[3]

  p[[i]] <- ggplot(data=df, aes(x = years,y=Bsim))+geom_line()+theme_classic()+
  geom_line(aes(y=Bobs),linetype = 2)+scale_x_continuous('year',breaks = years,labels = xtick)+
    scale_y_log10('Biomass')+# , breaks = YtiKbarents[,i])+
    coord_cartesian(ylim = c(y_limit_min,y_limit_max))+
    labs(title = nm[i],size = 9)+
    geom_rect(aes(xmin=1992,xmax=2002,ymin = 0,ymax =1e10), alpha = 0.005, color='gray', linetype = 0)+
    theme(axis.text.x = element_text(angle = 90, hjust = 1))+
    theme(text = element_text(size=7),
          text=element_text(family="Times"))#+
    #geom_line(data=state[i,],aes(x=1991:2002,y=seq(SSBio[1],length.out =12)), color='red')

}




return(p)
}


