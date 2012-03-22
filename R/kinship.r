kin <- function(gpData,ret=c("add","kin","dom","gam","realized","realizedAB","sm","sm-smin"),DH=NULL){

    ret <- match.arg(ret,choices=c("add","kin","dom","gam","realized","realizedAB","sm","sm-smin"),several.ok = FALSE)

    # (1) expected relatedness
    
    if (ret %in% c("add","kin","dom","gam")){

    # check for 'gpData'
    if(any(class(gpData)=="gpData")){
      if (is.null(gpData$pedigree)) stop("no pedigree found")
      else ped <- gpData$pedigree 
    }  
    
    # number of ids
    n <- nrow(ped)
    if(is.null(DH)) DH <- rep(0,n)
    if(!is.null(DH) & (length(DH) != n)) stop("DH must have same length as pedigree")

    # set up extended pedigree
    ID <- rep(seq_along(ped$ID),each=2)
    par1 <- pmatch(ped$Par1,ped$ID,nomatch = 0, duplicates.ok = TRUE)
    par2 <- pmatch(ped$Par2,ped$ID,nomatch = 0, duplicates.ok = TRUE)

    # set up gametic pedigree data.frame
    gamMat <- matrix(data=0,nrow=n*2,ncol=3,byrow=FALSE)
    gamMat[,1] <- ID

    # loop over ID
    for (i in 1:n){
        par1gam <- par1[i]
        par2gam <- par2[i]
        j <- (i-1)*2 + 1
        k <- j + 1
        #  parents of male genome contribution
        if(par1gam > 0){
           gamMat[j,2] <- (par1gam - 1)*2 + 1
           gamMat[j,3] <- (par1gam - 1)*2 + 2
        }
        #  parents of female genome contribution
        if(par2gam > 0){
           gamMat[k,2] <- (par2gam - 1)*2 + 1
           gamMat[k,3] <- (par2gam - 1)*2 + 2
        }
    }  # end of loop over ID

    #  Build Gametic Relationship
    ngam <- 2*n
    DHgam <- rep(DH,each=2)
    G <- diag(ngam)
    dimnames(G) <- list(paste(rep(ped$ID,each=2),rep(1:2,times=n),sep="_"), paste(rep(ped$ID,each=2),rep(1:2,times=n),sep="_"))

    # set inbreed coefficients of DHs on 1
    G[cbind((1:ngam)*DHgam,((1:ngam)+c(1,-1))*DHgam)] <- 1

    # caluclate gametic relationship
    # loop over gamets
    for(i in 1:(ngam-1-DHgam[2*n])){
      ip <- i+1 + (DHgam* rep(c(1,0),ngam))[i]

      for(j in ip:ngam){
          if(gamMat[j,2] > 0) {

            x <- 0.5*(G[i,gamMat[j,2]]+G[i,gamMat[j,3]])

            G[i,j] <- G[j,i] <- x
            }
    }
    } # end of loop over gamets

    # calculate addiditive and dominance relationship
    if(any(c("add","dom","kin") %in% ret)){
      A <- D <- matrix(data=NA,nrow=n,ncol=n)
      dimnames(A) <-  dimnames(D) <- list(ped$ID, ped$ID)

   # set up A and D matrices
   # loop over individuals
      for(i in 1:n){
         ka <- (i-1)*2 + 1
         for(j in i:n){
            kb <- (j-1)*2 + 1
            fab <- 0.25*(G[ka,kb]+G[ka,kb+1]+G[ka+1,kb]+G[ka+1,kb+1])
            A[i,j] <- A[j,i] <- 2*fab
            dab <- (G[ka,kb]*G[ka+1,kb+1] + G[ka+1,kb]*G[ka,kb+1])#*(1-G[ka,ka+1])*(1-G[kb,kb+1])
            # acoount for inbreeding
            # dominance = 0 if Fi=1
            D[i,j] <- D[j,i] <- dab
        }
      } # end of loop over individuals

      #diag(D) <- 1 - (diag(A)-1)
      
    }  # end of if

    # set return matrices
    if(ret == "add") kmat <- A
    if(ret == "dom") kmat <- D
    if(ret == "kin") kmat <- A/2
    if(ret == "gam") kmat <- G
    
    }
    
    # (2) realized relatedness
    
    if (ret == "realized"){ # former method vanRaden
    
        # extract information from arguments
          if(any(class(gpData)=="gpData")){
             if(!gpData$info$codeGeno) stop("use function 'codeGeno' before using 'kin'") 
             marker <- gpData$geno
          }
           else stop("object is not of class 'gpData'")

    # M supposed to be coded with 0,1,2
    M <- marker
    n <- nrow(M)
    p <- ncol(M)
    
    # 2* minor allele frequency as expectation
    maf <- colMeans(M,na.rm=TRUE)
    P <- matrix(rep(maf,each=n),ncol=p)
    
    # compute realized relationship matrix G
    Z <- M - P
    Zq <- tcrossprod(Z)
    U <- Zq/(2*sum(maf/2*(1-maf/2)))
    
    kmat <- U
    }

    if (ret == "realizedAB"){ # based an Astle & Balding (2009)
    
        # extract information from arguments
          if(any(class(gpData)=="gpData")){
             if(!gpData$info$codeGeno) stop("use function 'codeGeno' before using 'kin'") 
             marker <- gpData$geno
          }
           else stop("object is not of class 'gpData'")

        # M supposed to be coded with 0,1,2
        M <- marker
        n <- nrow(M)
        p <- ncol(M)
    
        # 2* minor allele frequency as expectation
        maf <- colMeans(M,na.rm=TRUE)
        pq2 <- 2*maf/2*(1-maf/2)
        # compute realized relationship matrix U
        Z <- sweep(M,2,maf)
        for (i in 1:p){  # loop for standardizing columns by sd
	    Z[,i]<-Z[,i]/sqrt(pq2[i])
        }
        U <- (Z %*% t(Z))/(2*p)
        kmat <- U
    }

    if (ret %in% c("sm","sm-smin")){      # simple matchin coefficient (only for homozygous inbreed lines)
              
          # extract information from arguments
          if(any(class(gpData)=="gpData")){
             if(!gpData$info$codeGeno) stop("use function 'codeGeno' before using 'kin'")
             if(any(gpData$geno == 1)) stop("simple matching coefficient is only for homozygous inbred lines")  
             marker <- gpData$geno
          }
          else stop("object is not of class 'gpData'")
            
          # code marker to -1/0/1 from 0,1,2
          marker <- marker - (max(marker,na.rm=TRUE)-1)
          m <- ncol(marker)
          
          # rogers distance
          d <- 1- (tcrossprod(marker) + m)/(2*m)
          
          #  simple matching coefficient
          if(ret=="sm-smin"){
            s <- 1-d
            smin <- min(s,na.rm=TRUE)
            f <- (s-smin)/(1-smin)
          }
          else f <- 1-d
          
          kmat <- 2*f
          
    }
    

    class(kmat) <- "relationshipMatrix"
    return(kmat)
}
