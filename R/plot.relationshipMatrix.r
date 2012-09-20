# heatmap for relationshipMatrix objects

plot.relationshipMatrix <- function(x,...){

    relMat <- x
    class(relMat) <- "matrix"
    size <- nrow(relMat)
    color <- c("#FFF7EC","#FEE8C8","#FDD49E","#FDBB84","#FC8D59","#EF6548","#D7301F","#B30000","#7F0000")
 
    if (size < 35){
      levelplot(relMat,axes=FALSE,col.regions=color,cuts=8,xlab="",ylab="",scales=list(cex=1-(size-20)/50,rot=c(40,0),abbreviate=TRUE,minlength=5),ylim=c(size+1,0),...)
    } else {
      levelplot(relMat,axes=FALSE,col.regions=color,cuts=8,xlab="",ylab="",scales=list(at=c(1,size/2,size),labels=c(1,"...",paste(size)),tck=0),ylim=c(size+1,0),...)  
    }                                   
}
