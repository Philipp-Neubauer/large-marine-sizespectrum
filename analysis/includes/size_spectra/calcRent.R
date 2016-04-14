calcRent <- function(param,S){

# parameters from Andersen et al 2015 (or closely)
a <- 2.5021225e4 # Extremely sensitive to this parameter
b <- 0.4 # Exponent that determines strength of species specific cost
c <- 0.41 # Exponent that determines weight specific price 
ap = 0.073 # factor for weight specific price 

w <- makegrid(param)
dw <- w[[2]]
w <- w[[1]]

N <- S$N[length(S$t),,]
# Backcalculate to spectrum 

Cost <-matrix(0,param$nSpecies)
profit <- matrix(0,param$nSpecies)


for (i in 1:param$nSpecies){
profit[i] <- sum(N[i,]*S$Fin[i,]*w*ap*w^c*dw)
Cost[i] <- a*param$F0[i]*param$wInf[i]^b
  
}

Rent <- sum(profit)-sum(Cost)  
df <- data.frame(profit = profit, Cost = Cost, wInf = t(param$wInf))

p <- ggplot(df, aes(x = wInf,y = profit/Cost))+theme_classic()+
  scale_x_log10(expand = c(0,0), name = 'Asymptotic size (g)', limit = c(min(param$wInf)/2,max(param$wInf)*2))+
  scale_y_log10(expand = c(0,0), name = 'Profit/Cost', limit = c(0.1,10))+
  geom_point(size = 5)+geom_line(size = 1)+geom_hline(yintercept = 1, alpha = 0.5)

p
return(list(Rent,profit,Cost,p))
}