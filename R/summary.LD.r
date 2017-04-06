# summary for LD objects

summary.LDdf <- function(object,cores=1,...){
  multiLapply <- function(x,y,...,cores=cores){
    if(.Platform$OS.type == "windows" & cores>1){
      cl <- makeCluster(min(cores, detectCores()))
      registerDoParallel(cl)
      parLapply(cl,x,y,...)
      stopCluster(cl)
    } else {
      mclapply(x,y,...,mc.cores=cores)
    }
  }
  LDdf <- object
  # help functions
  noM <- function(x) length(unique(x$marker1))+1
  avgr2 <- function(x) mean(x$r2,na.rm=TRUE)
  minr2 <- function(x) min(x$r2,na.rm=TRUE)
  maxr2 <- function(x) max(x$r2,na.rm=TRUE)
  pr <- function(x) mean(x$r2>0.2,na.rm=TRUE)
  maxd <- function(x) max(x$dist,na.rm=TRUE)
  ret <- data.frame(noM=unlist(multiLapply(LDdf,noM, cores=cores)),
                    avgr2= unlist(multiLapply(LDdf,avgr2, cores=cores)),
                    minr2= unlist(multiLapply(LDdf,minr2, cores=cores)),
                    maxr2= unlist(multiLapply(LDdf,maxr2, cores=cores)),
                    Pr02=unlist(multiLapply(LDdf,pr, cores=cores)),
                    averDist = unlist(multiLapply(LDdf,avgr2, cores=cores)),
                    maxDist = unlist(multiLapply(LDdf,maxd, cores=cores)))
     return(ret)
}

summary.LDmat <- function(object, cores=1,...){
  multiLapply <- function(x,y,...,cores=cores){
    if(.Platform$OS.type == "windows" & cores>1){
      cl <- makeCluster(min(cores, detectCores()))
      registerDoParallel(cl)
      parLapply(cl,x,y,...)
      stopCluster(cl)
    } else {
      mclapply(x,y,...,mc.cores=cores)
    }
  }
  LDmat <- object
  # help functions
  avgr2 <- function(x) mean(x[upper.tri(x)],na.rm=TRUE)
  minr2 <- function(x) min(x[upper.tri(x)],na.rm=TRUE)
  maxr2 <- function(x) max(x[upper.tri(x)],na.rm=TRUE)
  pr <- function(x) mean(x[upper.tri(x)]>0.2,na.rm=TRUE)
  ret <- data.frame(noM=unlist(multiLapply(LDmat$LD,ncol)),
                    avgr2= unlist(multiLapply(LDmat$LD,avgr2, cores=cores)),
                    minr2= unlist(multiLapply(LDmat$LD,minr2, cores=cores)),
                    maxr2= unlist(multiLapply(LDmat$LD,maxr2, cores=cores)),
                    Pr02=unlist(multiLapply(LDmat$LD,pr, cores=cores)),
                    averDist = unlist(multiLapply(LDmat$distance,avgr2, cores=cores)),
                    maxDist = unlist(multiLapply(LDmat$distance,max, cores=cores)))
  return(ret)
}
