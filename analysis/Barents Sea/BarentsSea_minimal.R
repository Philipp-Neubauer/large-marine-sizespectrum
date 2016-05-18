# Size based calibration to large marine ecosystems
# Calibrate to the Barents

files <- paste0('../includes/size_spectra/', dir('../includes/size_spectra'))
sapply(files, source)

db <- '../includes/RAM/'
# Barents sea

state <- Biomass_from_RAMv3(lmenum =20, db)
units <- state[[2]]
state <- state[[1]]

RAMnames <- state$fnames


# Sum all the same species 

cnames <- as.character(unique(state$cnames))

is.nan.data.frame <- function(x){
  do.call(cbind, lapply(x, is.nan))}

state[is.nan.data.frame(state)] <- NA
state$Biomass[state$Biomass == 0] = NA
state$SSB[state$SSB == 0] = NA



totBioidx <- which(is.na(state$Biomass) == 1)

state$SSBio <- state$SSB
state$totFlag <- as.numeric(is.na(state$SSBio))
state$SSBio[totBioidx] <- state$Biomass[totBioidx]

# Lme specific corrections 
state <- state[-6,] # Deepwater redfish has no info 
state$F0[6] <- state$Landings[6]/state$Biomass[6]
# Lets set up the size based model 

# Make parameters ready for the size based model 
# Sort the df by wInf 
state <- state[with(state, order(wInf)), ]

# Remove those that do not have any biomass denoted. 
#

rmidx <- which(is.na(state$SSBio) == 0)


state <- state[rmidx,] 
state$F0[which(is.na(state$F0) == 1)] <- state$Catch[which(is.na(state$F0) == 1)]/
  state$SSBio[which(is.na(state$F0) == 1)]

h <- 15

param <- baseparameters(state$wInf, kappa = 5e5, h = h) # Kappa estimated from LME 
param$tEnd <- 50
param$dt <- 0.1
S <- NA
param$F0 <- state$F0
param$fishing <- "Trawl"
ptm <- proc.time()

S <- IterateSpectrum(param,S)

proc.time()-ptm

pdf('plots/testSpectrum.pdf')
plotSpectrum(param,S)
dev.off()

plotBiomasstime(param,S)
plotBiomass(param,S)
# Plot the results
# Make sure the Rp/Rmax isn't too high
plotRpRmax(param,S)

# Calculate species specific values of maximum uptake (h)
KappaPP <- 5e5
h <- 3*state$k/(0.6*state$wInf^(-1/3))
 h[is.na(h)] <- mean(h,na.rm=T)
 cV <- 1/mean(h,na.rm=T)
# #param$h <- mean(param$h, na.rm = T)
# #h <- 20*matrix(1,40)
# param <- baseparameters(state$wInf,KappaPP,h = h)
# param$v <- param$h*cV
# param$fishing <- "Trawl"
# param$F0 <- state$F0
# param$tEnd <- 120
SF <- IterateSpectrum(param,S=S)
#
 plotBiomass(param,SF)
 plotBiomasstime(param,SF)
 plotRpR(param,SF)
#
# #kappaNew <- kappaCalc(state,S = SF,param = param,Rmax = param$Rmax) # Rerun for new calibration
 kappaNew <- 1391729
#
# # Plot the relative spectrum
#
plotSpectrum(param,SF)
#
#
# h <- 3*state$k/(0.6*state$wInf^(-1/3))
# h[is.na(h)] <- mean(h,na.rm=T)
# cV <- 1/mean(h,na.rm=T)
#
# param <- baseparameters(state$wInf,kappaNew,h = h)
# param$v <- param$h*cV
# param$F0 <- state$F0
# param$fishing <- "Trawl"
# param$tEnd <- 120
# param$nF <- 0.1
#
# SF <- IterateSpectrum(param,S=S)
# plotBiomasstime(param,SF)
# compareBiomass(param,SF,state)
# #Rmax_new <- RmaxCalc(state,param,SF)
#
# #Rmax <-Rmax_new
#
# # #
# # Rmax_news <- RmaxCalc_oneReturn(state,param,SF,param$Rmax)
# # Rmax <- as.numeric(exp(Rmax_new[1:param$nSpecies]))
# # write.table(Rmax,'Rmax_Barents.csv')
# #
#
# # Rmax <- exp(Rmax_new$par)
 Rmax <- as.numeric(as.matrix(read.table('../../Size based models LMEs/Rmax_Barents.csv')))
# param$Rmax <- Rmax
# param$eRepro <- matrix(0.05,param$nSpecies)
# param$eRepro[2] <- 0.05 # Very slow growing fish
# param$tEnd <- 120
# SF <- IterateSpectrum(param,S=SF)
# plotBiomasstime(param,SF)
#
# compareBiomass(param,SF,state)
# #Calculate Fmsy
# # Fmsy <- matrix(NA,param$nSpecies)
# # #
# # for (i in 1:param$nSpecies){
# #   Fmsy[i] <- findMSY(param = param,SF,i,state)
# #
# # }
#  Fmsy <- as.numeric(as.matrix((read.table('Fmsy_Barents.csv'))))
# #
#  cuts <- c(1,1,2,2,3,3)
# p <- calibrationPlots(state,param,SF,Fmsy,cuts)
#
#
#
# # Use the Barents to calculate the yield - efficiency frontier trade-off
# mvec <- seq(0,2,length.out = 15)
# p <- expand.grid(mvec,mvec,mvec) # Make 3 fleets
# ninit <- dim(p)[1]
# # Yield <- matrix(NA,ninit)
# # YieldSpecies <- matrix(NA,ninit,param$nSpecies)
# # SSB <- matrix(NA,ninit,param$nSpecies)
# # BM <- matrix(NA,ninit,param$nSpecies)
# #
# # for (i in 1:ninit){
# #
# #   param <- baseparameters(state$wInf,kappaNew,h = h)
# #   param$v <- param$h*cV
# #
# #   param$F0[1:2] <- as.numeric(p[i,1]*Fmsy[1:2])
# #   param$F0[3:4] <-as.numeric(p[i,2]*Fmsy[3:4])
# #   param$F0[5:6] <- as.numeric(p[i,3]*Fmsy[5:6])
# #
# #   param$fishing <- "Trawl"
# #   param$tEnd <- 50
# #   param$eRepro <- matrix(0.05,param$nSpecies)
# #   param$Rmax <- Rmax
# #   param$nF <- 0.1
# #
# #
# #   SF <- IterateSpectrum(param,SF)
# #   Yield[i] <- sum(YieldCalc(param,SF))
# #   YieldSpecies[i,] <- as.numeric(YieldCalc(param,SF))
# #   BM[i,] <- SF$Biomass[param$tEnd/param$dt,]
# #   SSB[i,] <- calcSSB(param,SF,param$tEnd/param$dt)
# #   print(i)
# #
# # }
# # #Save the calculation
# #  write.table(BM,'BarentsBiomassEF.csv')
# #  write.table(Yield,'BarentsYieldEF.csv')
# #  write.table(SSB,'BarentsSSBEF.csv')
# #  write.table(YieldSpecies,'BarentsYieldSpecies.csv')
# #
# #
#
# BM <- (as.matrix((read.table('BarentsBiomassEF.csv'))))
# Yield <- (as.matrix(read.table('BarentsYieldEF.csv')))
# SSB <- (as.matrix((read.table('BarentsSSBEF.csv'))))
# Rent <- NA
# YieldSpecies <- (as.matrix(read.table('BarentsYieldSpecies.csv')))
# source('plotFrontier.R')
# cuts <- c(1,1,2,2,3,3)
# # Make a list of calculated things
# BMall <- list(BM,Yield,SSB,Rent,YieldSpecies)
# plots <- plotFrontier(BMall,param,SF,20, year = 1980:2010,cuts)
#
# #
# #
# #  cairo_ps(file = "C:/Users/s102045/Dropbox/Ocean life/Size based models LMEs/Figures/Temporal figures/Barentstemp.eps",
#  #          width = 3.5*2, height = 3.5*2)
# # # # # #
#     pTemp <- multiplot(plotlist = plots[[1]], cols = 2)
#     dev.off()
# # # #
# # #  dev.off()
# # # # #
# #  cairo_ps(file = "C:/Users/s102045/Dropbox/Ocean life/Size based models LMEs/Figures/BarentsEF.eps",
# #           width = 3.3, height =1)
#  pBarents <- plots[[2]]+ theme(plot.margin = unit(c(0.2,0.1,-0.4,0), "cm"))
#   pBarents
# #
# #  dev.off()# ArrowBarents2 <- plots[[4]]
# # dev.off()
# # #
# #
# # #
# #  cairo_ps(file = "C:/Users/s102045/Dropbox/Ocean life/Size based models LMEs/Figures/BarentsCalibration.eps",
# #           width = 3.4*2,height =3.4*2.5)
# # # #
#    pCal <- calibrationPlots(state,param,SF,Fmsy,cuts)
# #  dev.off()
# # pdf(file = "C:/Users/s102045/Dropbox/Ocean life/Size based models LMEs/Figures/BarentsSpectrum.pdf",
# #     width = 3.5, height = 2)
# # pCal[[3]]
# # dev.off()
#
#
# # #
# #   cairo_ps(file =  "C:/Users/s102045/Dropbox/Ocean life/Size based models LMEs/Figures/ExpPlots/BarentsExP.eps",
# #           width = 3.30,height =2)
# #    cairo_ps(file =  paste(dir,"Figures/ExpPlots/BarentsExP.eps",sep = '/'),
# #             width = 3.30,height =2)
# #       plots[[8]][[1]]
# #   dev.off()
# # # #
# #
# # cairo_ps(file = "C:/Users/s102045/Dropbox/Ocean life/Size based models LMEs/Figures/BarentsRelSpec.eps",
# #          width = 3.4,height =2.5)
# # # #
# # plots[[6]]
# # dev.off()
# #
# # cairo_ps(file = "C:/Users/s102045/Dropbox/Ocean life/Size based models LMEs/Figures/Sheldon spectra/Barents_Sheldon.eps",
# #            width = 2.5,height =1.5)
# # #  # #
# # p <-  pCal[[3]][2]
# # p <- p[[1]] + coord_cartesian(ylim = c(1e1,1e7))
# # p
# # dev.off()
# #
#
# # Split the calibration plots
#    pCal <- calibrationPlots(state,param,SF,Fmsy,cuts)
#
# pCal[[1]] <- pCal[[1]]+annotate("text",x=35,y=7e5,label='A',size = 3)
#
# pCal[[2]] <- pCal[[2]]+annotate("text",x=4.4,y=5.8,label='C',size = 3)
#
# pCal[[4]] <- pCal[[4]]+annotate("text",x=0.15,y=1,label='B',size = 3)
# pCal[[5]] <- pCal[[5]]+annotate("text",x=4.2,y=5.8,label='D',size = 3)
#
#
# # cairo_pdf(file = "C:/Users/s102045/Dropbox/Ocean life/PhD latex fil/BarentsSplit.pdf",
# #                width = 4,height =3)
# #
# # multiplot(plotlist = list(pCal[[1]],pCal[[2]],pCal[[4]],pCal[[5]]), cols = 2)
# # dev.off()
#
# pCal <- calibrationPlots(state,param,SF,Fmsy,cuts)
# pCal[[6]] <- pCal[[6]]+annotate("text",x=1,y=3,label='A',size = 3)
#
# pCal[[7]] <- pCal[[7]]+annotate("text",x=0.7,y=0.7,label='C',size = 3)
#
# pCal[[8]] <- pCal[[8]]+annotate("text",x=35,y=30,label='E',size = 3)
# pCal[[10]] <- pCal[[10]]+annotate("text",x=5,y=8e6,label='B',size = 3)
# pCal[[11]] <- pCal[[11]]+annotate("text",x=0.001,y=0.9,label='D',size = 3)
# pCal[[12]] <- pCal[[12]]+annotate("text",x=35,y=0.98,label='F',size = 3)
#
# # cairo_pdf(file = "C:/Users/s102045/Dropbox/Ocean life/PhD latex fil/BarentsSplit2.pdf",
# #           width = 4,height =3.5)
# #
# # multiplot(plotlist = pCal[c(6,7,8,10,11,12)], cols = 2)
# # dev.off()
# # source('plotTraitFrontier.R')
# # sensBar <- plotTraitFrontier(BMall,param,SF,cuts)
# # sensBar
#
#
#
#
# # Potential code to fix plot labels
#
# # Consider using " scale_x_log10(breaks = trans_breaks("log10", function(x) 10^x),
# # labels = trans_format("log10", math_format(10^.x))" '
# # to get powers of 10 on the log-axis instead of E-notation