"[.relationshipMatrix" <- function(x,...) {
   y <- NextMethod("[")
   class(y) <- oldClass(x)
   y
}