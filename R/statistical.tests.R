#' Composite t-test for an Interval Null Hypothesis Difference in Means
#'
#' Performs a composite t-test to determine if the difference in means between
#' two numeric vectors `x` and `y` is greater than `epsilon_l` and less than `epsilon_u`.
#' The test involves performing two one-sided t-tests: one testing whether the difference
#' in means is less than `epsilon_l` and another testing if it is greater than `epsilon_u`.
#' If the p-values for each of these one-sided tests are `p_l` and `p_u`, respetively, then
#' the p-value for the composite t-test is calculated as `min(2*min(p_l,p_u), 1)`.
#' 
#' @param x A non-empty numeric vector of data values 
#' @param y A non-empty numeric vector of data values
#' @param epsilon_l A single numeric value representing the lower bound of the interval null hypothesis.
#' @param epsilon_u A single numeric value representing the lower bound of the interval null hypothesis.
#'                This must be larger than or equal to `epsilon_l`.
#' @param equal.var A logical value indicating whether to assume equal variances between the two populations.
#'                 If `TRUE`, the function uses the pooled variance; if `FALSE` it uses the Welch approximation.
#'                 The default is `FALSE`
#' @return A list containing the following components:
#' \describe{
#'   \item{p}{The p-value for the test.}
#'   \item{std.err}{The estimated standard error of the difference in means.}
#'   \item{m}{The observed difference in means between `x` and `y`.}
#' }
#' 
#' @importFrom stats pt sd var
#' @export
ttest.interval.fast <- function(x, y, epsilon_l, epsilon_u, equal.var=F) {
    t.res.l <- ttest.fast(x, y, epsilon_l, lower.tail=T, equal.var=equal.var)
    t.res.u <- ttest.fast(x, y, epsilon_u, lower.tail=F, equal.var=equal.var)
    p <- min(2*min(c(t.res.l$p, t.res.u$p)), 1)
    std.err <- ttest.fast(x, y, epsilon_l, lower.tail=T, equal.var=equal.var)$se
    return(list(p=p, std.err=std.err, m=t.res.l$m))
}

#' Composite Wilcoxon Rank Test for Interval Null Hypothesis
#' 
#' Performs a composite Wilcoxon rank test to assess whether the difference between two
#' numeric vectors, `x` and `y`, lies within a specified interval. Specifically, the test
#' determines if the difference is greater than `epsilon_l` and less than `epsilon_u`.
#' 
#' This is done by conducting two one-sided Wilcoxon rank-sum tests: one to test if the 
#' difference is less than `epsilon_l` and another to test if it is greater than `epsilon_u`.
#' The p-values from these tests are denoted as `p_l` and `p_u`, respectively. The p-value 
#' for the composite test is then computed as `min(2 * min(p_l, p_u), 1)`.
#' 
#' @param x A non-empty numeric vector of data values 
#' @param y A non-empty numeric vector of data values
#' @param epsilon_l A single numeric value representing the lower bound of the interval null hypothesis.
#' @param epsilon_u A single numeric value representing the lower bound of the interval null hypothesis.
#'                This must be larger than or equal to `epsilon_l`.
#' @param equal.var A logical value indicating whether to assume equal variances between the two populations.
#'                 If `TRUE`, the function uses the pooled variance; if `FALSE` it uses the Welch approximation.
#'                 The default is `FALSE`
#' @return A list containing the following components:
#' \describe{
#'   \item{p}{The p-value for the test.}
#'   \item{std.err}{The estimated standard error of the difference in means.}
#'   \item{m}{The observed difference in means between `x` and `y`.}
#' }
#' 
#' @importFrom stats wilcox.test
#' @export
wilcox.interval.test <- function(x, y, epsilon_l, epsilon_u) {
    lower.pval <- min(2*wilcox.test(x, y, mu=epsilon_l,
                                    alternative="less", correct=F, paired=F)$p.value, 1)
    upper.pval <- min(2*wilcox.test(x, y, mu=epsilon_u,
                                    alternative="greater", correct=F, paired=F)$p.value, 1)
    return(list(lower.pval=lower.pval, upper.pval=upper.pval))
}

#' Generalized t-test for an Interval Null Hypothesis of Difference in Means
#'
#' Performs a generalized t-test to determine if the difference in means between
#' two numeric vectors `x` and `y` lies within the specified interval
#' \[`epsilon_l`, `epsilon_u`\]. This test is an extension of the standard t-test,
#' incorporating an interval null hypothesis
#' 
#' @param x A non-empty numeric vector of data values 
#' @param y A non-empty numeric vector of data values
#' @param epsilon_l A single numeric value representing the lower bound of the interval null hypothesis.
#' @param epsilon_u A single numeric value representing the lower bound of the interval null hypothesis.
#'                This must be larger than or equal to `epsilon_l`.
#' @param equal.var A logical value indicating whether to assume equal variances between the two populations.
#'                 If `TRUE`, the function uses the pooled variance; if `FALSE` it uses the Welch approximation.
#'                 The default is `FALSE`
#' @return A list containing the following components:
#' \describe{
#'   \item{p}{The p-value for the test.}
#'   \item{std.err}{The estimated standard error of the difference in means.}
#'   \item{m}{The observed difference in means between `x` and `y`.}
#' }
#' 
#' @importFrom stats pt sd var
#' @export
generalized.ttest <- function(x, y, epsilon_l, epsilon_u, equal.var=F) {
    nx <- length(x)
    ny <- length(y)
    vx <- var(x)
    vy <- var(y)
    m <- mean(x) - mean(y)
    if(equal.var==T) {
        df <- nx+ny-2
        sp <- sqrt(vx*(nx-1) + vy*(ny-1)) / sqrt(df)
        se <- sp*sqrt(1/nx+1/ny) 
    } else {
        se <- sqrt(vx/nx+vy/ny)
        df <- ( (vx/nx+vy/ny)**2 ) / ( ((vx/nx)**2 )/(nx-1) + ((vy/ny)**2)/(ny-1) )
    }
    pval <- min(pt((m/se-epsilon_l/se), df=df) + pt((m/se-epsilon_u/se), df=df),
                1-pt((m/se-epsilon_l/se), df=df)+1-pt((m/se-epsilon_u/se), df=df))
    return(list(p=pval, std.err=se, m=m))
}

#' Perform a Composite t-test Using Precalculated Statistics
#'
#' This function performs a composite t-test similar to \code{t.test.interval.fast}, but allows the user 
#' to input precalculated values for the mean difference, standard error, and degrees of freedom.
#' 
#' @param m A numeric value representing the pre-estimated difference in means.
#' @param se A numeric value representing the pre-estimated standard error.
#' @param df A numeric value representing the pre-estimated degrees of freedom.
#' @param m One numeric value: a pre-estimated difference in means
#' @param se One numeric value: a pre-estimated standard error
#' @param df One numeric value: a pre-estimated degrees of freedom
#' 
#' @return A list containing the following components:
#' \describe{
#'   \item{p}{The p-value for the test.}
#'   \item{std.err}{The estimated standard error of the difference in means.}
#'   \item{m}{The observed difference in means between `x` and `y`.}
#' }
#' 
#' @importFrom stats pt
#' @export
ttest.interval.precalc <- function(m, se, df, epsilon_l, epsilon_u) {
    t_stat_l <- (m - epsilon_l) / se
    t_stat_u <- (m - epsilon_u) / se
    p_l <- pt(t_stat_l, df=df, lower.tail=T)
    p_u <- pt(t_stat_u, df=df, lower.tail=F)
    p <- min(2*min(c(p_l, p_u)), 1)
    return(list(p=p, m=m))
}

#' Perform a Generalized t-test Using Precalculated Statistics
#'
#' This function performs a composite t-test similar to \code{generalized.t.test}, but allows the user 
#' to input precalculated values for the mean difference, standard error, and degrees of freedom.
#' 
#' @param m A numeric value representing the pre-estimated difference in means.
#' @param se A numeric value representing the pre-estimated standard error.
#' @param df A numeric value representing the pre-estimated degrees of freedom.
#' @param m One numeric value: a pre-estimated difference in means
#' @param se One numeric value: a pre-estimated standard error
#' @param df One numeric value: a pre-estimated degrees of freedom
#' 
#' @return A list containing the following components:
#' \describe{
#'   \item{p}{The p-value for the test.}
#'   \item{std.err}{The estimated standard error of the difference in means.}
#'   \item{m}{The observed difference in means between `x` and `y`.}
#' }
#' 
#' @importFrom stats pt
#' @export
generalized.ttest.precalc <- function(m, se, df, epsilon_l, epsilon_u) {
    pval <- min(pt((m/se-epsilon_l/se), df=df) + pt((m/se-epsilon_u/se), df=df),
                1-pt((m/se-epsilon_l/se), df=df)+1-pt((m/se-epsilon_u/se), df=df))
    return(list(p=pval, m=m))
}

#' @importFrom stats pt sd var
#' @keywords internal
generalized.ttest.mehring <- function(Z, conds, epsilon_l, epsilon_u) {
    epsilon <- epsilon_u - mean(c(epsilon_l, epsilon_u))
    n <- sum(conds==1)
    m <- sum(conds==0)
    X <- Z[which(conds==1)]
    Y <- Z[which(conds==0)]
    k <- n+m-2
    md <- mean(X) - mean(Y) - mean(c(epsilon_l, epsilon_u))
    sigma <- sqrt(var(X) * (n-1) + var(Y) * (m-1) ) / sqrt( k )
    T_stat <- ( (md) * sqrt(n*m) ) / (sigma * sqrt(n+m))
    gamma_2 <- (epsilon/sigma) * sqrt((n*m)/(n+m))
    pval <- 1 - pt(abs(T_stat)-gamma_2, df=k) + pt(-abs(T_stat)-gamma_2, df=k)
    return(list(m=mean(X)-mean(Y), p=pval))
}

#' @importFrom stats pt sd var
#' @keywords internal
ttest.fast <- function(x, y, mu, lower.tail, equal.var=F) {
    nx <- length(x)
    ny <- length(y)
    vx <- var(x)
    vy <- var(y)
    m <- mean(x) - mean(y) 
    if(equal.var==T) {
        df <- nx+ny-2
        sp <- sqrt(vx*(nx-1) + vy*(ny-1)) / sqrt( df )
        se <- sp*sqrt(1/nx+1/ny) 
    } else {
        se <- sqrt(vx/nx+vy/ny)
        df <- ( (vx/nx+vy/ny)**2 ) / ( ((vx/nx)**2 )/(nx-1) + ((vy/ny)**2)/(ny-1) )
    }
    t_stat <- (m - mu) / se
    return(list(t_stat=t_stat, p=pt(t_stat, df=df, lower.tail=lower.tail), se=se, m=m))
}

#' @keywords internal
assign.pval <- function(t.res, epsilon_l, epsilon_u) {
    if(t.res$m<epsilon_l) {
        return(list(l=t.res$p, u=1))
    } else if(t.res$m>epsilon_u) {
        return(list(l=1, u=t.res$p))
    } else {
        return(list(l=t.res$p, u=t.res$p))
    }
}

