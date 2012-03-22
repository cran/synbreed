# predictions for objects of class gpMod

predict.gpMod <- function(object,newdata=NULL,...){
  if(class(object)!="gpMod") stop("'object' must be of class 'gpMod'")
  if (is.null(newdata)) prediction <- object$g
  else{
  model <- object$model

  if(model %in% c("BL","BRR")){
      if(!is.null(object$kin)) ("including a polygenic effect is not yet implemented")
      if(class(newdata)!="gpData") stop("object 'newdata' must be of class 'gpData'")
      X <- newdata$geno
      m <- object$m
      mu <- object$fit$mu
      prediction <- mu + X %*% m
  }
  if(model == "BLUP" & !is.null(object$m)){    # if marker effects are available
     if(class(newdata)!="gpData") stop("object 'newdata' must be of class 'gpData'")
     X <- newdata$geno
     m <- object$m
     mu <- object$fit$beta
     prediction <- mu + X %*% m
  }
  if(model == "BLUP"){
#      prediction <- gpData$geno %*% t(gpData$geno[rownames(kin), ]) %*% ginv(kin) %*% genVal[rownames(kin)]
#      prediction <- prediction[!names(prediction) %in% names(genVal)] / mean(prediction[names(genVal)]/genVal)

      G <- object$kin[c(names(object$y),newdata),c(names(object$y),newdata)]
      y <- object$y
      n <- length(y) # size of the training set
      Z <- cbind(diag(n),matrix(data=0,nrow=n,ncol=length(newdata)))
      X <- matrix(1,nrow=n)
      sigma2g <- object$fit$sigma[1]
      sigma2  <- object$fit$sigma[2]
      #diag(G) <- diag(G) + 0.00001
      GI <- ginv(G)*sigma2/sigma2g    # to ensure a solution for the inverse
      RI <- diag(n)
      sol <- MME(X, Z, GI, RI, y)
      prediction <- sol$b + sol$u[-(1:n)]
      names(prediction) <- newdata
  }

  }
  return(prediction)
}
