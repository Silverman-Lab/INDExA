#' @export
setClass("indexa.mc",
         slots=c(
             reads="matrix",
             dirichlet.samples="array",
             log2.norm.samples="array",
             scales="matrix",
             denom="ANY"
         )
)

validIndexaMC <- function(object) {
    ## check 0 dirichlet check dirichlet sum 1
    if(any(abs(apply(object@dirichlet.samples, c(2,3), sum)-1)>(1*10^-10))) {
        return("dirichlet sample cols must sum to 1")
    }
    if(!ncol(object@scales)==ncol(object@reads)) {
        return("ncol scales must be same as ncol reads")
    }
    if(!nrow(object@scales)==dim(object@dirichlet.samples)[3]) {
        return("ncol scales must be same as num mc samples")
    }
    if(!is.null(row.names(object@reads))) {
        if(sum(duplicated(row.names(object@reads)))!=0) {
            return("sequence count data contained duplicate row names")
        }
    }
    if(!is.null(colnames(object@reads))) {
         if(sum(duplicated(colnames(object@reads)))!=0) {
            return("sequence count data contained duplicate col names")
        }      
    }
    if((length(dim(object@dirichlet.samples))!=3)|
       (length(dim(object@log2.norm.samples))!=3)) {
        return(paste0("dirichlet.samples and log2.norm.samples",
                      " must have three dimensions"))
    }
    if(any(dim(object@reads)!=dim(object@dirichlet.samples)[c(1,2)])) {
        return(paste0("dim of reads must match first two dims ",
                      "of dirichlet.samples"))
    }
    if(any(dim(object@reads)!=dim(object@log2.norm.samples)[c(1,2)])) {
        return(paste0("dim of reads must match first two dims ",
                      "of log2.norm.samples"))
    }
    if(any(dim(object@dirichlet.samples)!=dim(object@log2.norm.samples))) {
        return(paste0("dims of dirichlet.samples must dims of",
                      " log2.norm.samples"))
    }
    return(TRUE)
}

setValidity("indexa.mc", validIndexaMC)

