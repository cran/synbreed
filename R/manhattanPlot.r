manhattanPlot <- function(b,gpData=NULL,colored=FALSE,add=FALSE,pch=19,ylab=NULL,colP=NULL,...){
  if(is.null(gpData)) plot(b,...)
  else{
    if(is.null(gpData$map)) stop("missing map in gpData object ",substitute(gpData))
    if (class(b) == "gpMod") b <- b$markerEffects
    b <- b[!(is.na(gpData$map$pos) | is.na(gpData$map$chr))]
    gpData$map <- gpData$map[!(is.na(gpData$map$pos) | is.na(gpData$map$chr)),]
    if(is.null(colP)){
      if(colored) colP <- rainbow(6)
      else colP <- rep(c(grey(0.3),grey(0.7)),times=length(unique(gpData$map$chr)))
      colP=colP[(as.numeric(gpData$map$chr)-1)%%6+1]
    }
    chrs <- cumsum(tapply(gpData$map$pos, gpData$map$chr, max))
    namChrs <- names(chrs)
    chrs <- c(0,chrs[1:(length(chrs)-1)])
    names(chrs) <- namChrs
    chr <- as.numeric(chrs[gpData$map$chr]) + as.numeric(gpData$map$pos)+ as.numeric(as.factor(gpData$map$chr))*0.01
    if(!add){
      plot(chr,b,col=colP,type="p",axes=FALSE,pch=pch,ylab=ylab,...)
      axis(side=1,at=c(chr[!duplicated(gpData$map$chr)],max(chr,na.rm=TRUE)),labels=NA, cex=.9, lwd.ticks=1,...)
      axis(side=1,at=chr[!duplicated(gpData$map$chr)]+diff(c(chr[!duplicated(gpData$map$chr)],max(chr,na.rm=TRUE)))/2,
           tick=FALSE,labels=unique(gpData$map$chr),hadj=0, padj=0,...)
      axis(side=2,...)
      box()
      } else points(chr, b, col=colP[(as.numeric(gpData$map$chr)-1)%%6+1], pch=pch,...)
  }
}
