LDDist <- function(LDdf,chr=NULL,type="p",breaks=NULL,file=NULL,n=NULL,...){


    if(class(LDdf)!="LDdf") stop("'LDdf' must be of class 'LDdf'")
    if(type=="nls" & is.null(n)) stop("number of obeservations must be specified using argument 'n'")
    if(is.null(chr)) lg <- (1:length(LDdf))[!as.logical(lapply(LDdf,is.null))]
    else lg <- chr
  

    # function for fit according to Hill and Weir (1988)
    smooth.fit <- function(overallDist,overallr2,n){
      # nls estimate
      nonlinearoverall <- nls(overallr2 ~ ((10 + p*overallDist)) / ((2+p*overallDist) * (11 + p*overallDist) ) *
      ( 1 + ( (3+ p*overallDist) * (12 + 12 * p + p^2*overallDist^2)) / ( n*(2+p*overallDist) * (11 + p*overallDist))),
      start=list(p=1))

      # estimated value for p
      p <- coef(nonlinearoverall)
      x <- NA
      fitcurve <- function(x,p,n) {
        ((10 + p*x)) / ((2+p*x) * (11 + p*x) ) *
        ( 1 + ( (3+ p*x) * (12 + 12 * p + p^2*x^2)) / ( n*(2+p*x) * (11 + p*x)))
      }
      # fit curve to data
      curve(fitcurve(x,p=p,n=n), from=min(overallDist), to = max(overallDist), add=TRUE,col=2,lwd=2,...)
    }
    
    # use LD from input arguement
    ret <- LDdf

    if(!is.null(file)) pdf(file)
    # compute distances within each linkage group
    for (i in lg){
    
       # create plots
       # scatterplot
       if(type=="p") plot(r2~dist,data=ret[[i]],main=names(ret)[[i]],...)

       # scatterplot with nls curve
       if(type=="nls"){
               plot(r2~dist,data=ret[[i]],main=paste("Linkage Group",i),...) 
               smooth.fit(ret[[i]][,4],ret[[i]][,3],n=n)
       }

       # stacked histogramm
       if(type=="bars"){
          # use default breaks
          if(is.null(breaks)){
             breaks.dist <- seq(from=min(ret[[i]]$dist),to=max(ret[[i]]$dist),length=6)
             breaks.r2 <- seq(from=1,to=0,by=-0.2) 
          }
          # use user- specified breaks
          else{
             breaks.dist <- breaks$dist
             breaks.r2 <- breaks$r2
          }
          cut.dist <- cut(ret[[i]]$dist,breaks=breaks.dist)
          cut.r2 <- cut(ret[[i]]$r2,breaks=breaks.r2)
          
          # create matrix with relative frequencies
          tab.abs <- table(cut.r2,cut.dist)
          colSum <- matrix(rep(colSums(tab.abs),nrow(tab.abs)),nrow=nrow(tab.abs),byrow=TRUE)
          
          # baplot out of frequency matrix
          barplot((tab.abs/colSum)[nrow(tab.abs):1,],col=grey(1:nrow(tab.abs)/nrow(tab.abs))[1:nrow(tab.abs)],space=c(.2),main=names(ret)[[i]],xlim=c(0,ncol(tab.abs)+2.8),ylab="Fraction of SNP pairs",...)
          legend(ncol(tab.abs)+1.2,0.95,fill=grey(1:nrow(tab.abs)/nrow(tab.abs))[nrow(tab.abs):1],legend=levels(cut.r2),title="LD (r2)",cex=1)
         }
        }
     # close graphic device
     if(!is.null(file)) dev.off()
}
