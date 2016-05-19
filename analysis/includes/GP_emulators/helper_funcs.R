run_SS_sim <- function(inputs){
  this.setup <- set_trait_model(no_sp = no_sp, 
                                r_pp = inputs$r_pp,
                                kappa = 0.01,
                                h=inputs$h,
                                beta = 300,
                                sigma= inputs$sigma,
                                q=0.9,#inputs$q,
                                p=0.75,#inputs$p,
                                n=2/3,#inputs$n,
                                f_0 = 0.8,
                                min_w_inf = min_w_inf, 
                                max_w_inf = max_w_inf,
                                knife_edge_size = knife_edges)
  
  this.sim <- project(this.setup, t_max = 100, dt=0.2, effort = 0.4)
  this.biom <- colMeans(getBiomass(this.sim)[75:100,])
  this.biom
}

run_nis_spec <- function(x, state) {
  
  param = baseparameters(state$wInf,
                         x$kappa,
                         h = as.numeric(x[grepl('\\h',names(x))]))
  cV <- 1/mean(param$h,na.rm=T)
  param$v <- param$h*cV
  param$F0 <- state$F0
  param$fishing <- "Trawl"
  param$tEnd <- 120
  param$nF <- 0.1
  param$Rmax =  as.numeric(x[grepl('Rmax',names(x))])
  SF <- IterateSpectrum(param, S=NA)
  tEnd <- param$tEnd/param$dt
  SSB <- calcSSB(param,SF,tEnd)
  Biomass <- SF$Biomass[tEnd,]
  
  SSBio <- matrix(NA,param$nSpecies)
  SSBio[state$totFlag == 0] <- SSB[state$totFlag == 0]
  SSBio[is.na(SSBio)] <- Biomass[which(is.na(SSBio) == 1)]
  SSBio
}

run_SS_sim_comp <- function(inputs){
  this.setup <- set_trait_model(no_sp = no_sp, 
                                r_pp = inputs$r_pp,
                                kappa = inputs$kappa,
                                h=unlist(inputs[grepl('h.',colnames(inputs))]),
                                k0 = unlist(inputs[grepl('k0.',colnames(inputs))]),
                                beta = inputs$beta,
                                sigma= inputs$sigma,
                                q=0.9,#inputs$q,
                                p=0.75,#inputs$p,
                                n=2/3,#inputs$n,
                                f_0 = 0.8,
                                min_w_inf = min_w_inf, 
                                max_w_inf = max_w_inf,
                                knife_edge_size = knife_edges,
                                gear_names = sprintf('ke %d',1:no_sp))
  
  this.sim <- project(this.setup, t_max = 50, dt=0.2,effort = effort)
  this.biom <- colMeans(getBiomass(this.sim)[25:50,])
  this.biom
}

lseq <- function(mins,maxs,l) 10^(seq(log10(mins),log10(maxs),l=l))

jitter_preds <- function(preds, keep=NULL){
  pred_pre_list <- expand.grid(preds)
  pred_pre_list_jitter <- sapply(1:ncol(pred_pre_list),function(x) {
    x_old <- pred_pre_list[,x]
    new_x <- x_old+rnorm(length(x_old),0,sd(x_old)/5)
    if (!is.null(keep)){
      sapply(new_x, function(z) min(max(z,keep[x,2]),keep[x,3]))} else {
      new_x[new_x<0] <- x_old[new_x<0]
      new_x
      }
  })
 # print(head(pred_pre_list_jitter))
  colnames(pred_pre_list_jitter) <- c('h','r_pp','sigma')
  pred_pre_list_jitter
}

repmat = function(X,m,n){
  ##R equivalent of repmat (matlab)
  X<- as.matrix(X)
  mx = dim(X)[1]
  nx = dim(X)[2]
  matrix(t(matrix(X,mx,nx*n)),mx*m,nx*n,byrow=T)
}

input_reg <- function(inputs) {
  out <- cbind(1,inputs,inputs^2,inputs^3)
  rownames(out) <- sprintf('input %d',1:nrow(as.matrix(inputs)))
  as.matrix(out)
}


output_reg <- function(outputs) {
  out <- cbind(1,outputs,outputs^3)
  rownames(out) <- sprintf('outputs %d',1:nrow(as.matrix(outputs)))
  as.matrix(out)
}


require(emulator)

input_var <- function(inputs,finputs=NULL,scales) {
  cmat <- corr.matrix(inputs,yold=finputs,scales=scales)
  if(!is.null(finputs)) cmat <- t(cmat)
  cmat
}


define_OPE <- function(outscale,inscales, inData, outGrid, outData, opt=F){
  
  input_var <- function(inputs,finputs=NULL) {
    cmat <- corr.matrix(inputs,yold=finputs,scales=inscales)
    if(!is.null(finputs)) cmat <- t(cmat)
    cmat
  }
  
  output_var <- corr.matrix(as.matrix(outGrid),scales=outscale)
  output_regs <- output_reg(outGrid)
  
  vr <- ncol(input_reg(inData))
  vs <- ncol(output_regs)
  ms <- rep(0, vr * vs)
  V <- diag(1^2, vr * vs)
  a <- 1
  d <- 1
  NIG <- list(m = ms, V = V, a = a, d = d)
  
  OPE <- initOPE(gr = input_reg, 
                 kappar = input_var, 
                 Gs = output_regs , 
                 Ws = output_var, 
                 NIG = NIG)
  
  OPE <- adjustOPE(OPE, R = inData, Y = outData)
  if(opt) {
    return(marlikOPE(OPE))
  } else {
    OPE
  }
  
}

m <- as.matrix
stdise <- function(x) apply(x,2,function(y) (y-mean(y))/sd(y))
stdise_wrt <- function(x,z) sapply(1:ncol(x), function(y) (x[,y]-mean(z[,y]))/sd(z[,y]))

#rev_stdise <- function(y) apply(y,2,function(x) x*sd(x)+mean(x))
rev_stdise_wrt <- function(y,z) sapply(1:ncol(y), function(x) y[,x]*sd(z[,x])+mean(z[,x]))

emulate <- function(ins,emulator,split=F,rev_std_out=T,rev_std_data=NULL){
  
  this.pred <- predictOPE(emulator, Rp = as.matrix(ins), type='EV')
  if (rev_std_out) this.pred$mu <- rev_stdise_wrt(this.pred$mu,rev_std_data)
  dim(this.pred$Sigma) <- rep(length(this.pred$mu),2)
  out <- data.frame(mu=as.vector(this.pred$mu),
                    Sigma=diag(this.pred$Sigma),
                    ind = 1:nrow(this.pred$mu))
  out <- arrange(out,ind)
  if (split) out <- split(out,as.factor(out$ind))
 
  out
}

reshape_preds <- function(pred){
  reshape2:::melt.list(pred,level = 1,id.vars = "ind") %>% 
    mutate(L1=as.numeric(L1),
           Ind = paste(L1,ind,sep='-')) %>% 
    group_by(Ind,variable) %>% 
    mutate(Species = 1:n()) %>% 
    ungroup() %>%
    select(L1,ind,Species,Ind,variable, value) %>% 
    arrange(L1,ind,Species) %>%
    tidyr::spread(variable,value) %>% 
    select(L1,ind,Species,Ind,mu,Sigma) %>% 
    arrange(L1,ind,Species)
}

get_impl <- function(prd, bioms, Vobs, Vdiscr) {
  
  prd %>% 
    group_by(L1,ind) %>% 
    summarise(impl = t(bioms-mu)%*%solve(Vobs+diag(Sigma))%*%(bioms-mu)) %>%
    ungroup %>%
    select(impl) %>%
    unlist() %>% as.vector()
  
}

#' Test marginal response of the emulator given a testset
test_em <- function(testset, 
                    emulator, 
                    variable) {
  
  test_data<- mclapply(split(testset,1:nrow(testset)),run_SS_sim,mc.cores = 8)
  
  testdata <- stdise_wrt(log(do.call('rbind',test_data)),log(sim_data))
  testdata <- reshape2::melt(testdata)
  colnames(testdata) <- c(variable,'Species','Biomass')
  testdata[variable] <- testset[variable]
  
  testset_st <- stdise_wrt(testset,m(pred_pre_list))
  
  testdata_em <- emulate(testset_st,
                         emulator,
                         split=T,
                         rev_std_out=F,
                         rev_std_data = log(sim_data))
  
  testdata_em <- do.call('rbind',testdata_em)
  
  colnames(testdata_em)[3] <- variable
  testdata_em[variable] <- rep(testset[[variable]],each=10)
  testdata_em$Species <- 1:10
  
  test <- inner_join(testdata_em,testdata)
  test$Species <- factor(test$Species)
  
  ggplot(test,aes_string(x=variable,y='Biomass')) +
    geom_line(aes(linetype='Biomass',col=Species)) +
    geom_ribbon(aes(y=mu,ymin=mu-sqrt(Sigma),ymax=mu+sqrt(Sigma),fill=Species),alpha=0.2)+
    geom_line(aes(y=mu,linetype='Predictions',col=Species)) + 
    facet_wrap(~Species) + 
    scale_linetype_discrete('') +
    theme_classic()
}

#' Test marginal response of the emulator given a testset
test_em_comp <- function(testset, 
                    emulator, 
                    variable,mult) {
  
  test_data<- mclapply(split(testset,1:nrow(testset)),run_SS_sim_comp,mc.cores = 4)
  
  testdata <- stdise_wrt(log(do.call('rbind',test_data))[,sp],log(sim_data))
  testdata <- reshape2::melt(testdata)
  colnames(testdata) <- c(variable,'Species','Biomass')
  testdata[variable] <- mult
  
  testset_st <- stdise_wrt(testset,m(pred_pre_list))
  
  testdata_em <- emulate(testset_st,
                         OPE_opt,
                         split=T,
                         rev_std_out=F,
                         rev_std_data = log(sim_data))
  
  testdata_em <- do.call('rbind',testdata_em)
  
  colnames(testdata_em)[3] <- variable
  testdata_em[variable] <- rep(mult,each=no_sp)
  testdata_em$Species <- 1:no_sp
  
  test <- inner_join(testdata_em,testdata)
  test$Species <- factor(test$Species)
  
  ggplot(test,aes_string(x=variable,y='Biomass')) +
    geom_line(aes(linetype='Biomass',col=Species)) +
    geom_ribbon(aes(y=mu,ymin=mu-sqrt(Sigma),ymax=mu+sqrt(Sigma),fill=Species),alpha=0.2)+
    geom_line(aes(y=mu,linetype='Predictions',col=Species)) + 
    facet_wrap(~Species) + 
    scale_linetype_discrete('') +
    theme_classic()
}