plotSpectrum <- function(param,S){
  
require(ggplot2)
w <- S$w
wInf <- param$wInf

Ntot <-  S$Ntot[param$tEnd/param$dt,w < max(param$wInf)*2]
NPP <- S$nPP[param$tEnd/param$dt,,]


p1 <- ggplot(data = data.frame(x = S$w[w < max(param$wInf)*2], y = Ntot), aes(x=x,y=y)) + geom_line(size = 1)+
  scale_y_log10(bquote('N (# '*g^-1*')'),expand = c(0,0), limits = c(min(Ntot[Ntot > 0]),max(NPP[S$wPP > 1e-2])))+
  scale_x_log10('weight (g)',expand = c(0,0),limits  = c(1e-2,max(param$wInf)))+theme_classic()+
  theme(text = element_text(size=9))

df_sheldon <- data.frame(x = S$w[w < max(param$wInf)*2], y = Ntot*S$w[w < max(param$wInf)*2]^2)
p2 <- ggplot(data = df_sheldon, 
             aes(x=x,y=y)) + geom_line(size = 1)+
  scale_y_log10('Sheldon spectrum (g)',expand = c(0,0), limits = c(min(df_sheldon$y[df_sheldon$y>0]),
                                                                     max(df_sheldon$y[df_sheldon$y > 1e-2]*10)))+
  scale_x_log10('weight (g)',expand = c(0,0),limits  = c(1e-2,max(param$wInf)))+theme_classic()+
  theme(text = element_text(size=9))

N <- S$N[param$tEnd/param$dt,,]

if (param$nSpecies> 0){
for (i in 1:param$nSpecies){
  
}
}else{
  p1 <- p1 + geom_line(data = data.frame(x = S$w[w <param$wInf[i]], 
                                         y = N[w<wInf[i]]), size = 0.5, alpha = 1/2, size = 0.5)
  p2 <- p2+ geom_line(data = data.frame(x = S$w[w <param$wInf[i]], 
                                        y = N[w<wInf[i]]*w[w<wInf[i]]^2), size = 0.5, alpha = 1/2, size = 0.5)
}
p1
p1 <- p1+ geom_line(data = data.frame(x = S$wPP, y = NPP), size = 1, linetype = 2)

p2 <- p2 +geom_line(data = data.frame(x = S$wPP, y = NPP*S$wPP^2), size = 1, linetype = 2)
return(p2)

}