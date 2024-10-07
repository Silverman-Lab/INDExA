#' Calculate Log-Fold Changes and Perform Interval Null Hypothesis Testing
#'
#' This function applies the CTT, GTT, and CWT tests in the simple case of
#' a comparison of two conditions.
#' 
#'
#' @param indexa.obj An object of class `indexa.mc`, as returned by the \code{indexa.mc} function.
#' @param conds A numeric or logical vector indicating the conditions. For a numeric vector, 
#'              values should be binary (1 for treatment or after condition, 0 for control or 
#'              before condition). For a logical vector, `TRUE` is equivalent to 1, and `FALSE` 
#'              is equivalent to 0.
#' @param epsilon.interval A numeric vector of length 2 specifying the interval bounds of the 
#'                         interval assumption (see examples below). The first element is the
#'                         lower bound, and the second element is the upper bound.
#' @param equal.var Logical; if `TRUE`, the t-tests assume equal variances between the two 
#'                  conditions. If `FALSE`, the Welch approximation is used.
#'
#' @return A `data.frame` containing the following columns:
#' \describe{
#'     \item{\code{lfc.lower.bound}}{Lower bound of the median log fold change
#'                                   across monte carlo samples.}
#'     \item{\code{lfc.lower.bound}}{Upper bound of the median log fold change
#'                                   across monte carlo samples.}
#'     \item{\code{ctt.pval}}{Expected p-value for the composite t-test across all Monte
#'                            Carlo samples for each gene/microbe.}
#'     \item{\code{cwt.pval}}{Expected p-value for the composite Wilcoxon rank test across
#'                            all Monte Carlo samples for each gene/microbe.}
#'     \item{\code{gtt.pval}}{Expected p-value for the generalized t-test across all Monte
#'                            Carlo samples for each gene/microbe.}
#'     \item{\code{ctt.pval.BH.adj}}{Expected Benjamini-Hochberg adjusted p-value for the
#'                                   composite t-test across all Monte Carlo samples for
#'                                   each gene/microbe.}
#'     \item{\code{cwt.pval.BH.adj}}{Expected Benjamini-Hochberg adjusted p-value for the
#'                                   composite Wilcoxon rank test across all Monte Carlo
#'                                   samples for each gene/microbe.}
#'     \item{\code{gtt.pval.BH.adj}}{Expected Benjamini-Hochberg adjusted p-value for the
#'                                   generalized t-test across all Monte Carlo samples for
#'                                   each gene/microbe.}
#' }
#'
#' @details #' This function computes the expected values (means) of p-values for various
#' statistical tests, including the composite t-test, generalized t-test, and 
#' composite Wilcoxon rank test, across all Monte Carlo samples in the provided 
#' `indexa.mc` object. It also computes the lower and upper bounds for the median
#' Log Fold Change (LFC) across all Monte Carlo Samples based on the interval
#' assumption given.
#' 
#' @note The median LFC is computed rather than the mean to accouny for potential
#' outliers in the distribution of LFCs (i.e., effect sizes) calculated across samples.
#' 
#' @examples
#' # Example: Suppose a study compares the absolute abundances of microbes in the colons
#' # of individuals in a treatment group (indicated by `1`) vs. a control group (indicated
#' # by `0`). Suppose we expect total colonic microbial load (i.e., the scale) to be up to
#' # 3 fold higher in the treatment group. We can represent this as the interval assumption
#' # [0, log2(2)] in our analysis. As long as this interval assumption is correct, our
#' # absolute differential abundance analysis will control the false positive rate.
#'
#' # Made-up dataset for colon study
#' Y <- rbind(c(4,0,1,9), c(5,1,9,0))
#' colnames(Y) <- c("control_1", "control_2", "treatment_1", "treatment_2")
#' rownames(Y) <- c("microbe_1", "microbe_2")
#'
#' ## Conditions: 0=control, 1=treatment
#' conds <- c(0, 0, 1, 1)
#'
#' ## Interval assumption (from example above)
#' epsilon.interval <- c(0, log2(3))
#' 
#' # Get the indexa object 
#' indexa.obj <- indexa.mc(Y, mc.samples=1000)
#' 
#' # Perform INDExA analysis
#' indexa.res <- indexa.test(indexa.obj, conds, epsilon.interval)
#' 
#' @author Kyle McGovern
#' 
#' @importFrom stats p.adjust
#' @export
indexa.test <- function(indexa.obj, conds, epsilon.interval, equal.var=T) {
    ## Get Monte Carlo Samples from Object
    mc <- getLog2NormSamples(indexa.obj)
    mc.samples <- dim(mc)[3]
    mc.sample <- mc[,,1]
    epsilon.l <- epsilon.interval[1]
    epsilon.u <- epsilon.interval[2]
    
    ## Conds Checks
    if((!is.vector(conds))|(length(conds)!=ncol(mc.sample))) {
        stop(paste0("'conds' must be a vector with the same length as the",
                    " number of samples in the sequence count reads matrix."))
    }
    if(!all(conds%in%c(1,0))) {
        if(!all(conds%in%c(TRUE,FALSE))) {
            if(!all(conds%in%c("1", "0"))) {
                stop(paste0("The 'conds' parameter must be either a binary",
                     " (1/0) numeric vector or a logical (TRUE/FALSE) vector."))
            }
        }
    }
    conds <- as.numeric(conds)
    
    ## epsilon check
    if (length(epsilon.interval) != 2) {
        stop(paste0("'epsilon.interval' must contain exactly 2 elements."))
    }
    if (epsilon.l > epsilon.u) {
        stop(paste0("The first element of 'epsilon.interval' cannot be",
                    " greater than the second."))
    }
    
    ## Iterate over monte carlo samples, run interval tests
    ctt.lower.pvals <- c()
    ctt.upper.pvals <- c()
    ctt.lower.adj.pvals <- c()
    ctt.upper.adj.pvals <- c()
    
    gtt.lower.pvals <- c()
    gtt.upper.pvals <- c()
    gtt.lower.adj.pvals <- c()
    gtt.upper.adj.pvals <- c()
    
    cwt.lower.pvals <- c()
    cwt.upper.pvals <- c()
    cwt.lower.adj.pvals <- c()
    cwt.upper.adj.pvals <- c()
    lfcs <- c()
    
    for(i in 1:mc.samples) {
        mc.sample <- mc[,,i]

        ## Iterate over each microbe (i.e., row)
        ctt.lower.pval <- c()
        ctt.upper.pval <- c()
        gtt.lower.pval <- c()
        gtt.upper.pval <- c()
        cwt.lower.pval <- c()
        cwt.upper.pval <- c()
        lfc <- c()
        for(j in 1:nrow(mc.sample)) {
            ## The interval assumption \(\theta^\perp \in [epsilon.l, pesilon.u]\)
            ## corresponds to the null hypothesis interval \([-epsilon.u, -epsilon.l]\)
            data <- mc.sample[j,]
            ctt.res <- ttest.interval.fast(data[conds==1], data[conds==0], -epsilon.u,
                                            -epsilon.l, equal.var=equal.var)
            ctt.lower.pval <- c(ctt.lower.pval, assign.pval(ctt.res, -epsilon.u, -epsilon.l)$l)
            ctt.upper.pval <- c(ctt.upper.pval, assign.pval(ctt.res, -epsilon.u, -epsilon.l)$u)
            gtt.res <- generalized.ttest(data[conds==1], data[conds==0], -epsilon.u, -epsilon.l,
                                          equal.var=equal.var)
            gtt.lower.pval <- c(gtt.lower.pval, assign.pval(gtt.res, -epsilon.u, -epsilon.l)$l)
            gtt.upper.pval <- c(gtt.upper.pval, assign.pval(gtt.res, -epsilon.u, -epsilon.l)$u)
            wilcox.res <- wilcox.interval.test(data[conds==1], data[conds==0],
                                               -epsilon.u, -epsilon.l)
            cwt.lower.pval <- c(cwt.lower.pval, wilcox.res$lower.pval)
            cwt.upper.pval <- c(cwt.upper.pval, wilcox.res$upper.pval)
            lfc <- c(lfc, gtt.res$m)
        }
        ctt.lower.pvals <- cbind(ctt.lower.pvals, ctt.lower.pval)
        ctt.upper.pvals <- cbind(ctt.upper.pvals, ctt.upper.pval)
        ctt.lower.adj.pvals <- cbind(ctt.lower.adj.pvals,
                                     p.adjust(ctt.lower.pval, method="BH"))
        ctt.upper.adj.pvals <- cbind(ctt.upper.adj.pvals,
                                     p.adjust(ctt.upper.pval, method="BH"))
        gtt.lower.pvals <- cbind(gtt.lower.pvals, gtt.lower.pval)
        gtt.upper.pvals <- cbind(gtt.upper.pvals, gtt.upper.pval)
        gtt.lower.adj.pvals <- cbind(gtt.lower.adj.pvals,
                                     p.adjust(gtt.lower.pval, method="BH"))
        gtt.upper.adj.pvals <- cbind(gtt.upper.adj.pvals,
                                     p.adjust(gtt.upper.pval, method="BH"))
        cwt.lower.pvals <- cbind(cwt.lower.pvals, cwt.lower.pval)
        cwt.upper.pvals <- cbind(cwt.upper.pvals, cwt.upper.pval)
        cwt.lower.adj.pvals <- cbind(cwt.lower.adj.pvals,
                                     p.adjust(cwt.lower.pval, method="BH"))
        cwt.upper.adj.pvals <- cbind(cwt.upper.adj.pvals,
                                     p.adjust(cwt.upper.pval, method="BH"))
        lfcs <- cbind(lfcs, lfc)
    }

    ## Aggregate results into single matrix for return
    ctt.all.pvals <- cbind(rowMeans(ctt.lower.pvals), rowMeans(ctt.upper.pvals))
    ctt.all.adj.pvals <- cbind(rowMeans(ctt.lower.adj.pvals),
                               rowMeans(ctt.upper.adj.pvals))
    gtt.all.pvals <- cbind(rowMeans(gtt.lower.pvals), rowMeans(gtt.upper.pvals))
    gtt.all.adj.pvals <- cbind(rowMeans(gtt.lower.adj.pvals),
                               rowMeans(gtt.upper.adj.pvals))
    cwt.all.pvals <- cbind(rowMeans(cwt.lower.pvals), rowMeans(cwt.upper.pvals))
    cwt.all.adj.pvals <- cbind(rowMeans(cwt.lower.adj.pvals),
                               rowMeans(cwt.upper.adj.pvals))
    lfc.eff <- apply(lfcs, 1, median)
    lfc.lower <- apply(lfcs, 1, median) + epsilon.l
    lfc.upper <- apply(lfcs, 1, median) + epsilon.u
    m <- cbind(lfc.eff, lfc.lower, lfc.upper,
               apply(ctt.all.pvals, 1, min), 
               apply(cwt.all.pvals, 1, min),
               apply(gtt.all.pvals, 1, min),
               apply(ctt.all.adj.pvals, 1, min), 
               apply(cwt.all.adj.pvals, 1, min),
               apply(gtt.all.adj.pvals, 1, min))
    colnames(m) <- c("median.effect", "median.effect.l", "median.effect.u", "ctt.pval",
                     "cwt.pval", "gtt.pval", "ctt.pval.BH.adj", "cwt.pval.BH.adj",
                     "gtt.pval.BH.adj")
    row.names(m) <- row.names(mc[,,1])
    return(m)
}

