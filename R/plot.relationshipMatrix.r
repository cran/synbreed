# heatmap for relationshipMatrix objects

#plot.relationshipMatrix <- function(x,limits,...){
plot.relationshipMatrix <- function(x,...){

    relMat <- x
    class(relMat) <- "matrix"
    size <- nrow(relMat)
    color <- c("#FFFFFF", "#FFF7EC", "#FEE8C8", "#FDD49E", "#FDBB84", "#FDAB68", "#FC8D5A",
               "#FB7C41", "#EF6548", "#ED5031", "#D73018", "#CC0000", "#B30000", "#A52A2A",
               "#990000", "#8B2323", "#7F0000", "#660000", "#4B1B06", "#341304", "#000000")
    Min <- min(relMat, na.rm=TRUE)
    Max <- max(relMat, na.rm=TRUE)
    if(Min == Max){
        levelbreaks <- c(Min-.1, Max+.1)
        color="#000000"
    } else {
        quantiles <- round(quantile(relMat, probs=c(.01, .99), na.rm=TRUE), digits=1)
        if(quantiles[1] == quantiles[2]){
            levelbreaks <- c(Min, quantiles[1], Max)
            color <- c("#FFFFFF", "#000000")
        } else if(quantiles[2] - quantiles[1] <.1){
            levelbreaks <- c(Min, quantiles[1],  quantiles[2], Max)
            color <- color[c(1,11,21)]
        } else {
            rangbreaks <- round(sum(abs(quantiles))/19*2, digits=1)*.5
            if(rangbreaks == 0){
                levelbreaks <- seq(from=quantiles[1], length.out=20, by=.05)
            } else {
                levelbreaks <- sort(c(seq(mean(quantiles), quantiles[1], -rangbreaks), seq(mean(quantiles)+rangbreaks, quantiles[2], rangbreaks)))
            }
            while(length(color)-1 > length(levelbreaks)){
                levelbreaks <- c(levelbreaks, max(levelbreaks)+rangbreaks)
            }
            levelbreaks <- levelbreaks[levelbreaks < Max]
            if(length(levelbreaks) >20) levelbreaks <- levelbreaks[1:20]
            levelbreaks <- c(Min, levelbreaks, Max)
            if(length(levelbreaks) <= length(color)) color <- color[c(1:(length(levelbreaks)-1), length(color))]
        }
    }
    if (size < 35){
      levelplot(relMat,axes=FALSE,
                col.regions=color,xlab="",ylab="",
                scales=list(cex=1-(size-20)/50,rot=c(40,0),abbreviate=TRUE,minlength=5),
                at=levelbreaks,
                ylim=c(size+.05,0.5),
                pretty=TRUE,...)
    } else {
      levelplot(relMat,axes=FALSE,
                col.regions=color,xlab="",ylab="",
                scales=list(at=c(1,size/2,size),labels=c(1,"...",paste(size)),tck=0),
                at=levelbreaks,
                ylim=c(size+.05,0.5),
                pretty=TRUE)#,...)
    }
}
