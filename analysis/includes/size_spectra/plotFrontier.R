plotFrontier <- function(BMall,param,SF_curr,LMEnumber,year,cuts){

# BMall = list that contains Biomass, yield, SSB and Rent
# param = parameters for the fitted model
# SF_curr = model run for fitted model
# LMEnumber = Number of large marine ecosystem plotted 
# Direction = Direction of the ecological velocity arrow
# year = time period considered for plotting
# cuts = fleet division vector, 1 = small, 2 = medium, 3 = large
  
source('yieldplot.R')

BM <- BMall[[1]]
Yield <- BMall[[2]]
SSB <- BMall[[3]]
#Rent <- BMall[[4]]
#Yield <- Rent # REMOVE FOR YIELD

  
BM0 <- BM[1,]
SSB0 <- SSB[1,]

BM_E <- matrix(NA,ninit,param$nSpecies)
SSB_E <- matrix(NA,ninit,param$nSpecies)
for (i in 1:ninit){
  BM_E[i,] <- as.numeric(BM[i,])/as.numeric(BM0)
  SSB_E[i,] <- SSB[i,]/SSB0
}
x0 <- 0.2
x <- 100

ESSB_cal <- 1-x^(-(SSB_E-x0))
# Weigh by the biomass 
Emean <- matrix(NA,ninit)
for (i in 1:ninit){
  Emean[i] <- mean(ESSB_cal[i,])
}



#plotBiomasstime(param,SF_curr)
#compareBiomass(param,SF_curr,state)

SSBcurrent <- as.numeric(calcSSB(param,SF_curr,param$tEnd/param$dt))
YieldCurrent <- as.numeric(YieldCalc(param,SF_curr))
RentCurrent <- as.numeric(calcRent(param,SF_curr)[1])
#YieldCurrent <- RentCurrent # REMOVE FOR YIELD
                          
Ecurrent <- (SSBcurrent)/SSB0
Ecurrent <- 1-x^(-(Ecurrent-x0))
Ecurrent <-mean(Ecurrent)

df <- data.frame(exploit = Emean, Yield = as.numeric(Yield))
dfCurrent <- data.frame(YieldCurrent=sum(YieldCurrent), Ecurrent = Ecurrent)

p2 <- ggplot(data=df, aes(x = Yield, y = exploit))+#geom_point(color = 'gray', alpha = 1)+
  theme_classic()+
  geom_point(data=dfCurrent,aes(x=sum(YieldCurrent),y=Ecurrent), 
                             color = 'black', size = 2, shape =3)
  #geom_point(aes(x=Yield[1],y=exploit[1]), color = 'blue', size = 5, shape = 18)+
  

# Economic frontier for the Baltic Sea
p2 <- ggplot(data=df, aes(x = Yield, y = exploit))+geom_point(color = 'gray', alpha = 0.00)+
  theme_classic()+
  geom_point(data=dfCurrent,aes(x=sum(YieldCurrent),y=Ecurrent), 
             color = 'black', size = 2, shape =3)+
  #geom_point(aes(x=Yield[1],y=exploit[1]), color = 'blue', size = 5, shape = 18)+
  theme(text = element_text(size=9))

parFmsy <- param
parFmsy$F0 <- Fmsy
SF_MSY <- IterateSpectrum(parFmsy,SF_curr)
SSB_MSY <- as.numeric(calcSSB(parFmsy,SF_MSY,parFmsy$tEnd/parFmsy$dt))
Yield_Fmsy <- as.numeric(YieldCalc(parFmsy,SF_MSY))
Rent_Fmsy <- as.numeric(calcRent(parFmsy,SF_MSY)[1])
#Yield_Fmsy <- Rent_Fmsy # REMOVE FOR YIELD

EFmsy <- (SSB_MSY)/SSB0
EFmsy <- 1-x^(-(EFmsy-x0))
EFmsy <-mean(EFmsy)
dfMsy <- data.frame(exploit = EFmsy,Yield = sum(Yield_Fmsy))

# Plot the biomass relative to unfished under Fmsy 
dfMsyplot <- data.frame(SSBMSY= SSB_MSY,SSB0=SSB0,wInf = t(param$wInf),EFmsy = EFmsy) 
dfCurrentplot <- data.frame(SSBcurrent=SSBcurrent,SSB0 = SSB0, wInf = t(param$wInf))

ant1 <- paste("e['MSY'] ==",signif(EFmsy[1], digits =2))
ant2 <-  paste("e['Calibration'] ==",signif(Ecurrent, digits =2))
pFmsy <- ggplot(dfMsyplot, aes(y= SSBMSY/SSB0,x=wInf))+geom_line()+geom_point()+
  geom_line(data=dfCurrentplot, aes(y=SSBcurrent/SSB0,x=wInf), linetype = 3)+
  geom_point(data=dfCurrentplot, aes(y=SSBcurrent/SSB0,x=wInf))+
  theme_classic()+
  scale_x_log10(name = "Asymptotic weight (g)")+
  scale_y_log10(name = expression(paste('SSB'['MSY'],' /SSB'[0])))+
  geom_hline(yintercept = 1)+
  geom_hline(yintercept = 0.2, linetype = 2)+
  annotate('text',x = param$wInf[param$nSpecies]/2, y = 0.8, 
           label = ant1, parse=TRUE, size = 9/2.8)+
  annotate('text',x = param$wInf[param$nSpecies]/2, y = 0.7, 
           label = ant2, parse=TRUE, size = 9/2.8)+ # 2.8 is the ridiculous conversion from theme to text size
  theme(text = element_text(size=9))
# Add the frontier point with the same yield as Fmsy 


pFmsy
p2 <- p2+geom_point(data=dfMsy, aes(x=Yield,y=exploit), color = 'gray', size = 4)
# Preferences 
pref <- high(df$Yield) * high(df$exploit)
df$idx <- 1:(dim(df)[1])
df_frontier <- psel(df,pref)
p2 <- p2 + #geom_point(data= df_frontier, aes(x = Yield, y = exploit), color = 'red')+
  geom_line(data= df_frontier, aes(x = Yield, y = exploit), color = 'red', size = 2)

Y_Fmsy <- sum(dfMsy$Yield)
ix <- which.min((Y_Fmsy - df_frontier$Yield)^2) # Index in the frontier with approx same yield 

df_frontMsy <- data.frame(wInf = state$wInf, SSB0 = SSB0, SSBfrontMsy = SSB[df_frontier$idx[ix],])

pFmsy <- pFmsy + geom_point(data=df_frontMsy, aes(x = wInf, y = SSBfrontMsy/SSB0), colour = 'gray')
# Plot the observed E-Y in the frontier plot 
state$SSYield <- state$Landings
state$SSYield[which(is.na(state$SSYield) == 1)] = state$Catch[which(is.na(state$SSYield) == 1)]

#
# See how the Barents Sea captures temporal biomasses and yields 

dfTemp <- getTemporalInfo(LMEnumber = LMEnumber,tStart = year[1],tEnd = year[length(year)])

# Remove the ones that are not in the calibration 


dfTemp <- AverageList(dfTemp,year)

rmidx <- logical(length(dfTemp))
for (i in  1:length(dfTemp)){
  idxtmp <- which(names(dfTemp)[i] == state$cnames)
  if (length(idxtmp) > 0){
    rmidx[i] <- TRUE
  }
  
}
dfTemp <- dfTemp[rmidx]
dfTemp <- SortDFtemp(state,dfTemp)

# Make the SSBiomass (either SSB or Biomass)
BiomassStart <- matrix(NA,length(dfTemp))
for (i in 1:param$nSpecies){
  dfTemp[[i]]$SSBiomass<- dfTemp[[i]]$SSB  
  dfTemp[[i]]$flag <- 0
  if (all(is.na(dfTemp[[i]]$SSBiomass))){
    dfTemp[[i]]$flag[is.na(dfTemp[[i]]$SSBiomass)] <- 1
    dfTemp[[i]]$SSBiomass[is.na(dfTemp[[i]]$SSBiomass)] <- dfTemp[[i]]$Biomass 
  }
  BMtemp <- dfTemp[[i]]$SSBiomass[dfTemp[[i]]$flag == 0]
  BiomassStart[i] <- BMtemp[1]
  # Correct with exploitation index if all F are 0 
  if (all(is.na(dfTemp[[i]]$F0)) == 1){
    dfTemp[[i]]$F0 <- dfTemp[[i]]$Landings/dfTemp[[i]]$Biomass
    dfTemp[[i]]$F0 <- dfTemp[[i]]$Catch/dfTemp[[i]]$Biomass
  }
  if(is.na(BiomassStart[i] == 1)){
    BiomassStart[i] <-dfTemp[[i]]$SSBiomass[which(is.na(dfTemp[[i]]$SSBiomass) == 0)][1]
  }
  
  if (length(dfTemp[[i]]$year) == length(year)){
  dfTemp[[i]] <- dfTemp[[i]][2:length(dfTemp[[i]]$year),]
  }# Remove the biomass used to initate the simulation
  

}


BiomassStart <- as.numeric(BiomassStart)


F1 <- NA
SSB <- NA
for (i in 1:length(dfTemp)){
  tmp <- as.data.frame(dfTemp[[i]])
  
  
  F1[i] <- tmp[1,1]
  SSB[i] <- tmp$SSBiomass[1]
  
}

naidx <- which(is.na(F1) == 1)
if (length(naidx) > 0){

  for (i in 1:length(naidx)){
    tmp <- dfTemp[[naidx[i]]]
F1[naidx[i]] <- tmp$F0[which(is.na(tmp$F0) == 0)][1]

if (is.na(F1[naidx[i]]) ==1){
Ftmp <- tmp$Landings/tmp$Biomass
F1[naidx[i]]<- Ftmp[which(is.na(Ftmp) == 0)][1]
    }
  }
      
      
      
if (LMEnumber == 7){
  dfTemp[12]$Weakfish$F0 <- dfTemp[12]$Weakfish$Catch/dfTemp[12]$Weakfish$SSB
  F1[12] <- 0.784
  #F1[18] <- 0.0597
  #BiomassStart[18] <- 18900 # Hack the cusk that decreases > 2 fold within one year (no info on F)
}

}
# 

    
year <- year[2:length(year)]   
param$F0 <- F1
parOld <- param
S92 <- Ninit_BM(param,SF_curr,BiomassStart,state)

S92_SSB <- calcSSB(param,S92,param$tEnd/param$dt)
#plot(S92_SSB)
#points(SSB,col='red')


SF <- S92
SSB <- matrix(NA,length(year),param$nSpecies)
Biomass <- matrix(NA,length(year),param$nSpecies)
Yieldtime <- matrix(NA,length(year),param$nSpecies)
RentTime <- matrix(NA,length(year))
M2Save <- matrix(NA,length(year), param$nSpecies)

FeqSave <- data.frame("y1980" = matrix(NA,param$nSpecies),"y2010" =matrix(NA,param$nSpecies))

F0.save <- matrix(NA,length(year),param$nSpecies)

for (i in 1:(length(year))){
  Fy <- matrix(NA,1,param$nSpecies)
  for (j in 1:length(dfTemp)){
    tmp <- dfTemp[[j]]
    yidx <- which(year[i] == tmp[,6])
    if (length(yidx > 0)){
      Fy[j] <- tmp[,1][yidx]  
    }else{
      Fy[j] <- F1[j]}
    if(is.na(Fy[j]) == 1){
      Fy[j] <- F1[j]
      #Fy[j] <- Fmsy[j]
      
      # Find the nearest fishing mortality
      
      
    }
  }
  Fy <- as.numeric(Fy)
  if(i == 1){param$tEnd = 1
  #FeqSave[,1] <- Fy
  }else{
  param$tEnd <- 1
  }
  
  F0.save[i,] <- Fy
  param$dt <- 0.1
  param$F0 <- Fy
  SF <- IterateSpectrum(param,SF) 
  M2Save[,i] <- mean(SF$M2[5,,])
  SSB[i,] <- t(calcSSB(param,SF,param$tEnd/param$dt))
  Biomass[i,] <- SF$Biomass[param$tEnd/param$dt,]
  Yieldtime[i,] <- YieldCalc(param,SF)
  RentTime[i] <- as.numeric(calcRent(param,SF)[1])
  
  if (i == length(year)){
    FeqSave[,2] <- Fy
  }
  
}
df_timeCalc <- as.data.frame(SSB)
df_timeCalcBM <- as.data.frame(Biomass)
df_Yieldtime <- as.data.frame(Yieldtime)
names(df_timeCalc) <- state$cnames
df_timeCalc$years <- year[1:length(year)]  
names(df_Yieldtime) <- state$cnames
df_Yieldtime$years <- year[1:length(year)]  

meanTmp <- NA
for (i in 1:length(dfTemp)){
  meanTmp[i] <- mean(dfTemp[[i]]$SSBiomass, na.rm = T)
}
# Calculate the E-value from year to year 
df_Etime <- matrix(NA,length(df_Yieldtime$years),param$nSpecies)

for (i in 1:param$nSpecies){
  df_Etime[,i] <- t(df_timeCalc[,i]/SSB0[i])
}
Emean = rowMeans(1-x^(-(df_Etime-x0)))


 df_E <- data.frame(Emean = Emean,
                    Yield = rowSums(df_Yieldtime[,1:param$nSpecies]), years =year[1:length(year)])  

 #df_E <- data.frame(Emean = Emean,
  #                Yield = RentTime, years =year[1:length(year)])  # Uncomment below for yield


Yield_vec = seq(min(df_E$Yield),max(df_E$Yield),length.out = 2)
E_vec <- seq(min(df_E$Emean),max(df_E$Emean), length.out = 2)

# Create a dataframe to make an arrow in the plot

p2 <- p2+geom_path(data=df_E,aes(x=Yield,y = Emean, label =years), color = 'black')+
  geom_text(data=df_E[c(1,length(df_E$years)),],aes(x=Yield,y = Emean, label =years),
            hjust=0, vjust=0, size = 2, color = 'black')+
  geom_point(data=df_E,aes(x=Yield,y = Emean, label =years), color = 'black')+
  theme(legend.position="none")



year <- c(min(year)-1,year)

p1 <- plotTemporalBiomass(df_timeCalc,df_timeCalcBM,dfTemp,year,parOld,BiomassStart,S92)

p3 <- yieldplot(dfTemp,year)

source('ArrowCalc.R')
dfArrowExp <- ArrowCalc(df_E)


source('eqArrow.R')
FeqSave[,1] <- colMeans(F0.save[1:5,])

dfArrowExp <- list(eqArrow(parOld,SF_curr,FeqSave,SSB0,x,x0))
# Calculate the arrow with eq values from 1980 to 2010 




# Calculate the frontier distance
p2 <- p2+geom_segment(data=dfArrowExp[[1]],aes(x=yStart,xend=yEnd,y=eStart,yend=eEnd),
                      arrow=arrow(length = unit(0.5,'cm')), size = 2.5,alpha = 1, color='grey')+
  scale_x_continuous(expand = c(0.01,0.01), name = "")+
  scale_y_continuous(expand = c(0.01,0.01), name = 'I')+coord_cartesian(ylim = c(0,1))
# p2 <- p2+geom_segment(data=dfArrowExp[[2]],aes(x=yStart,xend=yEnd,y=eStart,yend=eEnd),
#                       arrow=arrow(length = unit(0.5,'cm')), size = 2.5,alpha = 0.5, color='grey')
# p2 <- p2+geom_segment(data=dfArrowExp[[3]],aes(x=yStart,xend=yEnd,y=eStart,yend=eEnd),
#                       arrow=arrow(length = unit(0.5,'cm')), size = 2.5,alpha = 0.6, color='grey')

# Average point
#dfArrowExp<- data.frame(Yield = dfArrow$Yield/max(df_frontier$Yield),e = dfArrow$E)

source('CalcFrontierExploit.R')
expPlots <- CalcFrontierExploit(BMall,param,SF_curr,cuts)
source('CalcFrontierExploit_vs2.R')  
expPlots2 <- CalcFrontierExploit_vs2(BMall,param,SF_curr,cuts,parOld,FeqSave)
l2010 <- list(parOld,FeqSave)
return(list(p1,p2,p3,dfArrowExp,expPlots,pFmsy,df_Yieldtime,expPlots2,l2010, M2Save))
}
