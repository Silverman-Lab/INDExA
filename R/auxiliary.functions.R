#' Generate Random Samples from a Dirichlet Distribution
#'
#' @param n An integer specifying the number of samples to generate.
#' @param alpha A positive numeric vector representing the parameters of the Dirichlet distribution.
#' 
#' @return A matrix with `n` rows and length of `alpha` columns, where each column 
#' is a random sample drawn from the Dirichlet distribution.
#' 
#' @importFrom stats rgamma
#' @export
rdirichlet <- function(n, alpha) {
    if(any(alpha<0)) {
        stop("rdirichlet cannot contain negative elements")
    }
    res <- matrix(0, ncol=n, nrow=length(alpha))
    for(i in 1:n) {
        res[,i] <- sapply(alpha, function(a) rgamma(1, a))    
    }
    res <- apply(res, 2, function(col) col/sum(col))
    return(res)
}

#' Compute the Geometric Mean of a Numeric Vector
#'
#' @param x A numeric vector of positive values.
#' @param na.rm A logical indicating whether to remove all NA from `x`.
#' 
#' @return The geometric mean of the values in `x`.
#' 
#' @export
gm.mean <- function(x, na.rm=TRUE) {
    if(any(x<0)) {
        stop("gm.mean cannot take negative numbers")
    }
    exp(sum(log(x), na.rm=na.rm) / length(x))
}

