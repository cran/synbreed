# heatmap for relationshipMatrix objects

plot.relationshipMatrix <- function(x,y=NULL,levelbreaks=NULL,axes=TRUE,cols=NULL,...){
  oldPar <- par(no.readonly = TRUE)
  plotRelMatS <- function(x,levelbreaks=levelbreaks,axes=axes,cols=cols,...){
    relMat <- x[, nrow(x):1]
    class(relMat) <- "matrix"
    size <- nrow(relMat)
    if(is.null(cols))
      col <- c("#ffffff", "#ffe4c8", "#ffdab4", "#ffd0a0", "#ffc182", "#ffb76e", "#ffad5a",
               "#ffa346", "#ff9932", "#ff8f1e", "#ff850a", "#e17100", "#cd6700", "#b95d00",
               "#a55300", "#914900", "#7d3f00", "#5f3000", "#4b2600", "#371c00", "#000000")
    Min <- min(relMat, na.rm=TRUE)
    Max <- max(relMat, na.rm=TRUE)
    if(is.null(levelbreaks)){
      if(Min == Max){
        levelbreaks <- c(Min-.1, Max+.1)
        col="#000000"
      } else {
        quantiles <- round(quantile(relMat, probs=c(.01, .99), na.rm=TRUE), digits=1)
        if(quantiles[1] == quantiles[2]){
          levelbreaks <- c(Min, quantiles[1], Max)
          col <- c("#FFFFFF", "#000000")
        } else if(quantiles[2] - quantiles[1] <.1){
          levelbreaks <- c(Min, quantiles[1],  quantiles[2], Max)
          col <- col[c(1,11,21)]
        } else {
          rangbreaks <- round(sum(abs(quantiles))/19*2, digits=1)*.5
          if(rangbreaks == 0){
            levelbreaks <- seq(from=quantiles[1], length.out=20, by=.05)
          } else {
            levelbreaks <- sort(c(seq(mean(quantiles), quantiles[1], -rangbreaks), seq(mean(quantiles)+rangbreaks, quantiles[2], rangbreaks)))
          }
          while(length(col)-1 > length(levelbreaks)){
            levelbreaks <- c(levelbreaks, max(levelbreaks)+rangbreaks)
          }
          levelbreaks <- levelbreaks[levelbreaks < Max]
          if(length(levelbreaks) >20) levelbreaks <- levelbreaks[1:20]
          levelbreaks <- unique(c(Min, levelbreaks, Max))
        }
      }
    }
    if(is.list(levelbreaks)) levelbreaks <- levelbreaks[[1]]
    if(is.list(cols)) col <- cols[[1]]
    if(length(levelbreaks) <= length(col))
      col <- col[round(seq(1, length(col), length.out=length(levelbreaks)-1))]
    par(mar = c(1.5, 1, 2, 2.8) + 0.1)
    layout(matrix(c(2,1), ncol = 2), widths = c(0.83, 0.17))
    image(x=c(-.5,.5), y=levelbreaks, z=t(matrix((levelbreaks[-length(levelbreaks)]+levelbreaks[-1])/2)), col=col, axes=FALSE, xlab="", ylab="")
    box()
    axis(side = 4, las = 1)
    par(mar=c(5,4,4,1)+.1)
    image(relMat,col=col, xaxt="n", yaxt="n")
    box()
    if(axes){
      (size <- nrow(relMat))
      if (size < 35){
        axis(1, at=seq(0, 1, length.out=size), labels=FALSE); axis(2, at=seq(0, 1, length.out=size), labels=FALSE)
        text(x=seq(0, 1, length.out=size)-.4/size,labels=colnames(x), srt=40,xpd=TRUE, y=par()$usr[3]-0.075*(par()$usr[4]-par()$usr[3]))
        text(y=seq(1, 0, length.out=size)-.4/size,labels=rownames(x), srt=40,xpd=TRUE, x=par()$usr[3]-0.075*(par()$usr[4]-par()$usr[3]))
      } else {
        axis(1, at=c(0,1), labels=FALSE); axis(2, at=c(0,1), labels=FALSE)
        text(x=c(0,1)-.4/size,labels=colnames(x)[c(1,size)], srt=40,xpd=TRUE, y=par()$usr[3]-0.075*(par()$usr[4]-par()$usr[3]))
        text(y=c(1,0)-.4/size,labels=rownames(x)[c(1,size)], srt=40,xpd=TRUE, x=par()$usr[3]-0.075*(par()$usr[4]-par()$usr[3]))
        text(y=.5-.4/size,labels="...", srt=90,xpd=TRUE, x=par()$usr[3]-0.075*(par()$usr[4]-par()$usr[3]))
        text(x=.5-.4/size,labels="...", srt=0,xpd=TRUE, y=par()$usr[3]-0.075*(par()$usr[4]-par()$usr[3]))
      }
    }
    par(oldPar)
  }
  plotRelMatD <- function(x,y,levelbreaks,axes,cols,...){
    if(nrow(x)>nrow(y)) relMat1 <- relMat2 <- x[, ncol(x):1]*NA else  relMat1 <- relMat2 <- y[, ncol(y):1]*NA
    x[upper.tri(x, diag=TRUE)] <- NA
    y[lower.tri(y, diag=TRUE)] <- NA
    relMat1[rownames(y), colnames(y)] <- y
    relMat2[rownames(x), colnames(x)] <- x
    class(relMat1) <- class(relMat2) <- "matrix"
    size <- nrow(relMat1)
    col1 <- c("#ffffff", "#ffe4c8", "#ffdab4", "#ffd0a0", "#ffc182", "#ffb76e", "#ffad5a",
              "#ffa346", "#ff9932", "#ff8f1e", "#ff850a", "#e17100", "#cd6700", "#b95d00",
              "#a55300", "#914900", "#7d3f00", "#5f3000", "#4b2600", "#371c00", "#000000")
    col2 <- c("#ffffff", "#c8e4ff", "#b4daff", "#a0d0ff", "#82c1ff", "#6eb7ff", "#5aadff",
              "#46a3ff", "#3299ff", "#1e8fff", "#0a85ff", "#0071e1", "#0067cd", "#005db9",
              "#0053a5", "#004991", "#003f7d", "#00305f", "#00264b", "#001c37", "#000000")
    Min1 <- min(relMat1, na.rm=TRUE); Min2 <- min(relMat2, na.rm=TRUE);
    Max1 <- max(relMat1, na.rm=TRUE); Max2 <- max(relMat2, na.rm=TRUE)
    levelbreaks1 <- levelbreaks2 <- levelbreaks
    if(is.list(levelbreaks)) {
      levelbreaks1 <- levelbreaks[[1]]
      levelbreaks2 <- levelbreaks[[2]]
    }
    if(is.list(cols)){
      col1 <- cols[[1]]
      col2 <- cols[[2]]
    }
    if(is.null(levelbreaks1)){
      if(Min1 == Max1){
        levelbreaks1 <- c(Min1-.1, Max1+.1)
        col1="#000000"
      } else {
        quantiles <- round(quantile(relMat1, probs=c(.01, .99), na.rm=TRUE), digits=1)
        if(quantiles[1] == quantiles[2]){
          levelbreaks1 <- c(Min1, quantiles[1], Max1)
          col1 <- c("#FFFFFF", "#000000")
        } else if(quantiles[2] - quantiles[1] <.1){
          levelbreaks1 <- c(Min1, quantiles[1],  quantiles[2], Max1)
          col1 <- col1[c(1,11,21)]
        } else {
          rangbreaks1 <- round(sum(abs(quantiles))/19*2, digits=1)*.5
          if(rangbreaks1 == 0){
            levelbreaks1 <- seq(from=quantiles[1], length.out=20, by=.05)
          } else {
            levelbreaks1 <- sort(c(seq(mean(quantiles), quantiles[1], -rangbreaks1), seq(mean(quantiles)+rangbreaks1, quantiles[2], rangbreaks1)))
          }
          while(length(col1)-1 > length(levelbreaks1)){
            levelbreaks1 <- c(levelbreaks1, max(levelbreaks1)+rangbreaks1)
          }
          levelbreaks1 <- levelbreaks1[levelbreaks1 < Max1]
          if(length(levelbreaks1) >20) levelbreaks1 <- levelbreaks1[1:20]
          levelbreaks1 <- unique(c(Min1, levelbreaks1, Max1))
        }
      }
    }
    if(is.null(levelbreaks2)){
      if(Min2 == Max2){
        levelbreaks2 <- c(Min2-.1, Max2+.1)
        col1="#000000"
      } else {
        quantiles <- round(quantile(relMat2, probs=c(.01, .99), na.rm=TRUE), digits=1)
        if(quantiles[1] == quantiles[2]){
          levelbreaks2 <- c(Min2, quantiles[1], Max2)
          col1 <- c("#FFFFFF", "#000000")
        } else if(quantiles[2] - quantiles[1] <.1){
          levelbreaks2 <- c(Min2, quantiles[1],  quantiles[2], Max2)
          col1 <- col1[c(1,11,21)]
        } else {
          rangbreaks2 <- round(sum(abs(quantiles))/19*2, digits=1)*.5
          if(rangbreaks2 == 0){
            levelbreaks2 <- seq(from=quantiles[1], length.out=20, by=.05)
          } else {
            levelbreaks2 <- sort(c(seq(mean(quantiles), quantiles[1], -rangbreaks2), seq(mean(quantiles)+rangbreaks2, quantiles[2], rangbreaks2)))
          }
          while(length(col1)-1 > length(levelbreaks2)){
            levelbreaks2 <- c(levelbreaks2, max(levelbreaks2)+rangbreaks2)
          }
          levelbreaks2 <- levelbreaks2[levelbreaks2 < Max2]
          if(length(levelbreaks2) >20) levelbreaks2 <- levelbreaks2[1:20]
          levelbreaks2 <- unique(c(Min2, levelbreaks2, Max2))
        }
      }
    }
    if(length(levelbreaks1) <= length(col1))
      col1 <- col1[round(seq(1, length(col1), length.out=length(levelbreaks1)-1))]
    if(length(levelbreaks2) <= length(col2))
      col2 <- col2[round(seq(1, length(col2), length.out=length(levelbreaks2)-1))]
    par(mar = c(2.5, 1.5, 1.5, 2.8) + 0.1)
    layout(matrix(c(3,3,1,2), ncol = 2), widths = c(0.87, 0.13))
    image(x=c(-.5,.5), y=levelbreaks2, z=t(matrix((levelbreaks2[-length(levelbreaks2)]+levelbreaks2[-1])/2)), col=col2, axes=FALSE, xlab="", ylab="")
    box()
    axis(side = 4, las = 1)
    image(x=c(-.5,.5), y=levelbreaks1, z=t(matrix((levelbreaks1[-length(levelbreaks1)]+levelbreaks1[-1])/2)), col=col1, axes=FALSE, xlab="", ylab="")
    box()
    axis(side = 4, las = 1)
    par(mar=c(5,4,4,2)+.1)
    image(relMat1,col=col1, xaxt="n", yaxt="n")
    image(relMat2,col=col2, add=TRUE, xaxt="n", yaxt="n")
    box()
    if(axes){
      (size <- nrow(relMat1))
      if (size < 35){
        axis(1, at=seq(0, 1, length.out=size), labels=FALSE); axis(2, at=seq(0, 1, length.out=size), labels=FALSE)
        text(x=seq(0, 1, length.out=size)-.4/size,labels=colnames(relMat1), srt=40,xpd=TRUE, y=par()$usr[3]-0.075*(par()$usr[4]-par()$usr[3]))
        text(y=seq(0, 1, length.out=size)-.4/size,labels=rownames(relMat1), srt=40,xpd=TRUE, x=par()$usr[3]-0.075*(par()$usr[4]-par()$usr[3]))
      } else {
        axis(1, at=c(0,1), labels=FALSE); axis(2, at=c(0,1), labels=FALSE)
        text(x=c(0-1.2^(-size),.99),labels=colnames(relMat1)[c(size,1)], srt=40,xpd=TRUE, y=par()$usr[3]-0.05*(par()$usr[4]-par()$usr[3]))
        text(y=c(1-1.2^(-size),.01),labels=rownames(relMat1)[c(1,size)], srt=40,xpd=TRUE, x=par()$usr[3]-0.05*(par()$usr[4]-par()$usr[3]))
        text(y=.5-.4/size,labels="...", srt=90,xpd=TRUE, x=par()$usr[3]-0.05*(par()$usr[4]-par()$usr[3]))
        text(x=.5-.4/size,labels="...", srt=0,xpd=TRUE, y=par()$usr[3]-0.05*(par()$usr[4]-par()$usr[3]))
      }
    }
    par(oldPar)
  }

 if(is.null(y)) {
   plotRelMatS(x,levelbreaks,axes,cols,...)
 } else if(class(y)[1] != "relationshipMatrix") {
   plotRelMatS(x,y,levelbreaks,axes,cols,...)
 } else {
   plotRelMatD(x,y,levelbreaks,axes,cols,...)
 }
}
