write.relationshipMatrix <- function(relationshipMatrix,file=NULL,sorting=c("WOMBAT","ASReml"),type=c("ginv","inv","none"),digits=10){
        
        type <- match.arg(type)
        sorting <- match.arg(sorting)

        if(sorting=="WOMBAT" & type!="ginv") stop("'type' must be 'ginv' for WOMBAT")
        
        # pass (inverse) relationship matrix
        if(type=="ginv") rMinv <- ginv(relationshipMatrix)
      	if(type=="inv")  rMinv <- solve(relationshipMatrix)
        if(type=="none") rMinv <- relationshipMatrix
        
        rMinv <- round(rMinv,digits)
        
        # add rownames and colnames
        res <- data.frame(coeff = as.numeric(rMinv),
                        rowname = rep(1:nrow(rMinv), nrow(rMinv)),
                        colname = rep(1:nrow(rMinv), each = nrow(rMinv)),
                        lower = as.logical(lower.tri(rMinv, diag = TRUE)))
                        
      
    
        # only use lower triangle
        res <- res[res$lower == TRUE, c("rowname", "colname", "coeff")]
          
        if (sorting=="ASReml"){    
          res <-  res[order( res$rowname,  res$colname), ] 
        }
        if (sorting=="WOMBAT"){
          res <- res[,c(2,1,3)]
          res <-  res[order( res$colname,  res$rowname), ]  
        }
        
        # write to given file
        if (!is.null(file)) write.table(res, file, sep = " ", quote = FALSE, col.names = FALSE, row.names = FALSE)

        else return(res)


}
