#' Compute an \code{indexa.mc} Object using Monte Carlo Sampling
#'
#' @aliases indexa.mc
#' 
#' Generates Monte Carlo samples from the Multinomial-Dirichlet distribution for each sample 
#' in a given matrix of sequence counts. These samples provide estimates of the microbial 
#' or gene composition, which are then used to compute estimates of the scale and absolute 
#' abundance. The resulting \code{indexa.mc} object serves as input for further analyses.
#'
#' @param reads A matrix of measured sequence counts, where each row represents a gene,
#' microbe, or taxon, and each column represents a measured sample.
#' @param mc.samples An integer specifying the number of Monte Carlo samples to generate.
#' DEFAULT: 250
#' @param denom A parameter for specifying the method used to calculate the scale. Options include:
#' \describe{
#'   \item{"none"}{(Default) No normalization. Use this option when specifying an interval over 
#'                 the LFC in scales. For example, if the total microbial load (i.e., scale) is expected
#'                 to increase by up to 2.5-fold between treatment and control conditions, use 
#'                 "none" and specify the interval \code{c(0, log2(2.5))} in the next step.}
#'   \item{"all"}{Equivalent to Centered Log-Ratio (CLR) normalization as in ALDEx2. Use this 
#'                when specifying an interval assumption to account for error in CLR normalization.}
#'   \item{numeric vector}{Specifies row indices in the `reads` matrix (i.e., the genes/microbes) to
#'                         use for normalization. This option allows normalization with housekeeping
#'                         genes, under the assumption these genes or microbes have LFCs of `0`. An
#'                         interval assumption in the next step would account for errors in this
#'                         assumption.}
#' }
#' @param custom.scale A custom function for specifying the scale, intended for advanced use cases. 
#' This function should accept a three-dimensional array (with dimensions corresponding to the number 
#' of rows in \code{reads}, the number of samples, and \code{mc.samples}) and return a matrix. The
#' returned matrix should have `mc.samples` rows and the same number of columns as the `reads` matrix.
#' The returned matrix of the custom function should be log2 values representing the scale in each sample
#' for each monte carlo sample. The output of the custom matrix is then added to the log2 estimated
#' proportions to obtain the log2 absoute abundance values used in the next step. See examples below
#' for how to define `custom.scale`.
#'
#' @return An \code{indexa.mc} object.
#'
#' @examples
#'    # The 'reads' matrix has microbes/genes for the rows
#'    # and samples for the columns. Like:
#'    #
#'    #            sample_1  sample_2  sample_3 ...etc...
#'    # microbe_1         0         3         7
#'    # microbe_2        12        23         0
#'    # microbe_3       212         0        98
#'    # ... etc ...
#' 
#' # Simple sequence count data
#' Y <- rbind(c(4,0,1,9), c(5,1,9,0))
#' colnames(Y) <- c("control_1", "control_2", "treatment_1", "treatment_2")
#' rownames(Y) <- c("microbe_1", "microbe_2")
#'
#' # Example 1: No normalization, log2 transformed proportions
#' indexa.obj <- indexa.mc(Y, mc.samples=1000)
#'
#' # Example 2: CLR normalization (default in ALDEx2)
#' indexa.obj.clr <- indexa.mc(Y, denom="all", mc.samples=1000)
#'
#' # Example 3: Custom scale function (equivalent to CLR normalization)
#' custom.func <- function(dirichlet.samples) {
#'   return(-log2(apply(dirichlet.samples, c(3,2), gm.mean)))
#' }
#' indexa.obj.clr <- indexa.mc(Y, custom.scale=custom.func, mc.samples=1000)
#'
#' # Example 4: Custom scale function (equivalent to denom="none")
#' custom.func <- function(dirichlet.samples) {
#'   return(matrix(0, nrow=dim(dirichlet.samples)[3], ncol=dim(dirichlet.samples)[2]))
#' }
#' indexa.obj.clr <- indexa.mc(Y, custom.scale=custom.func, mc.samples=1000)
#'
#' @export
indexa.mc.function <- function(reads, mc.samples=250, denom="none", custom.scale=NULL) {
    ## Checks
    if(any(rowSums(reads==0)==ncol(reads))) {
        if(is.null(row.names(reads))) {
            all_zero_samples <- which(rowSums(reads==0)==ncol(reads))
        } else {
            all_zero_samples <- row.names(reads)[which(rowSums(reads==0)==ncol(reads))]
        }
        warning(paste("The 'reads' matrix contains rows with all zero values.",
                      "Please remove the following rows for better results:",
                      paste(all_zero_samples, collapse=", ")))
    }
    if(mc.samples<100) {
        warning(paste0("Using ", mc.samples, " Monte Carlo samples. It is ",
                       "recommended to use at least 100 for reliable results."))
    }
    
    ## Dirichlet resampling
    dirichlet.samples <- array(0, dim=c(nrow(reads), ncol(reads), mc.samples))
    colnames(dirichlet.samples) <- colnames(reads)
    row.names(dirichlet.samples) <- row.names(reads)
    for(m in 1:mc.samples) {
        for(col in 1:ncol(reads)) {
            dirichlet.samples[,col,m] <- rdirichlet(1, reads[,col]+0.5)
        }
    } 

    ## Build Scale Matrix
    if(!is.null(custom.scale)) {
        if(!is.function(custom.scale)) {
            stop("The 'custom.scale' parameter must be a function.")
        }
        scales <- custom.scale(dirichlet.samples)
        if((nrow(scales)!=mc.samples)|(ncol(scales)!=ncol(reads))) {
            stop(paste0("'custom.scale' should return a matrix with dimensions",
                        " (", mc.samples, ", ", ncol(reads), "), but returned a matrix",
                        " with dimensions (", dim(scales)[1], ", ", dim(scales)[2], ")."))
        }
    } else if (!is.null(denom)) {
        if(is.numeric(denom)) {
            scales <- -log2(apply(dirichlet.samples, c(3,2),
                                  function(col) gm.mean(col[denom])))
        } else if(length(denom)>1) {
            stop("'denom' must be 'none', 'all', or a numeric vector.")
        } else if (denom=="none") {
            scales <- matrix(0, nrow=mc.samples, ncol=ncol(reads))
        } else if (denom=="clr") {
            scales <- -log2(apply(dirichlet.samples, c(3,2), gm.mean))
        } else {
            stop("'denom' must be 'none', 'all', or a numeric vector.")
        }
    } else {
        stop("'custom.scale' and 'denom' cannot both be NULL.")
    }
    colnames(scales) <- colnames(reads)
    
    ## Peform scale transformation of dirichlet samples
    log2.norm.samples <- array(0, dim=c(nrow(reads), ncol(reads), mc.samples))
    colnames(log2.norm.samples) <- colnames(reads)
    row.names(log2.norm.samples) <- row.names(reads)
    for(col in 1:ncol(reads)) {
        for(m in 1:mc.samples) {
            log2.norm.samples[,col,m] <- log2(dirichlet.samples[,col,m]) + scales[m,col]
        }
    }
    
    return(new("indexa.mc", reads=reads, dirichlet.samples=dirichlet.samples,
               log2.norm.samples=log2.norm.samples, scales=scales, denom=denom))
}

setMethod("getReads", signature(.object="indexa.mc"),
          function(.object) .object@reads)
setMethod("getDirichletSamples", signature(.object="indexa.mc"),
          function(.object) .object@dirichlet.samples)
setMethod("getLog2NormSamples", signature(.object="indexa.mc"),
          function(.object) .object@log2.norm.samples)
setMethod("getScales", signature(.object="indexa.mc"),
          function(.object) .object@scales)
setMethod("getDenom", signature(.object="indexa.mc"),
          function(.object) .object@denom)
setMethod("indexa.mc", signature(reads="matrix"), function(reads, mc.samples=250, denom="none", custom.scale=NULL) indexa.mc.function(reads, mc.samples, denom, custom.scale))

