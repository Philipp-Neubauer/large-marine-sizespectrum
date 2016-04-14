SortDFtemp <- function(state,dfTemp){
 
  rmidx <- matrix(NA,length(dfTemp))
  # Remove the ones we dont use 
  for (i in 1:length(dfTemp)){
    
    idxtmp <- which(names(dfTemp)[i] == state$cnames)
    
    if ((length(idxtmp)) > 0){
      rmidx[i] <- idxtmp}
    else{
      rmidx[i] <- NA}
  }

  # Remove the ones without information
  rmidx <- as.numeric(rmidx)
  No_idx <- which(is.na(rmidx) == 1)
  if (length(No_idx)> 0){
  dfTemp <- dfTemp[-No_idx]
  }

  
  
  
  Sortidx <- NA
  for (i in 1:length(dfTemp)){
    idxtmp <- which(state$cnames[i] == names(dfTemp))
    Sortidx[i] <- idxtmp
  }
  
  dfTemp <- dfTemp[Sortidx]
  
return(dfTemp)
  
}