summary.relationshipMatrix <- function(object,...){
     relMat <- object
     ans <- list(dim=c(nrow=nrow(relMat),ncol=ncol(relMat)),rank=qr(relMat)$rank,range.off.diagonal=c(min=min(relMat[upper.tri(relMat,diag=FALSE)]),max=max(relMat[upper.tri(relMat,diag=FALSE)])),mean.diag=mean(diag(relMat)),mean.off.diag=mean(relMat[lower.tri(relMat,diag=FALSE)]),nUnique=length(unique(relMat[upper.tri(relMat,diag=TRUE)])),diag.val=summary(diag(relMat)))
     class(ans) <- "summary.relationshipMatrix"
     ans
}

print.summary.relationshipMatrix <- function(x,...){
    cat(" dimension                   ",x$dim[1],"x",x$dim[2],"\n")
    cat(" rank                        ",x$rank,"\n")
    cat(" range of off-diagonal values",x$range.off.diagonal[1],"--",x$range.off.diagonal[2],"\n")
    cat(" mean off-diagonal values    ",x$mean.off.diag,"\n")
    cat(" range of diagonal values    ",x$diag.val[1],"--",x$diag.val[6],"\n")
    cat(" mean diagonal values        ",x$mean.diag,"\n")
    cat(" number of unique values     ",x$nUnique,"\n")
}