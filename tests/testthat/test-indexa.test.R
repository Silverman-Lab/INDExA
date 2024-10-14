library(ALDEx2)
library(abind)

set.seed(1)
D <- 100
N <- 30
Intercepts <- runif(D, 3, 8)
LFC <- c(rep(0, 0.25*D), runif(0.75*D,-4,4))
conds <- c(rep(0, N/2), rep(1, N/2))
m <- matrix(0, ncol=N, nrow=D)
for(n in 1:N) {
    for(d in 1:D) {
            m[d,n] <- rnorm(1, Intercepts[d]+(LFC[d]*conds[n]), 0.75)
    }
}
true_scale <- mean(log2(colSums(2^m))[((N/2)+1):N]) - mean(log2(colSums(2^m))[1:(N/2)])
Y <- apply(m, 2, function(col){ rmultinom(1, D*25, (2^col)/sum(2^col)) })

clr.obj <- aldex.clr(Y, conds, mc.samples=300, denom="all", useMC=FALSE)
mc <- ALDEx2::getMonteCarloInstances(clr.obj)
dir.d <- getDirichletInstances(clr.obj)
mc.d <- getMonteCarloInstances(clr.obj)
dir_array <- c()
mc_array <- c()
for(i in 1:numMCInstances(clr.obj)) {
    dir_array <- abind(dir_array, sapply(dir.d, function(x) x[, i]), along=3)
    mc_array <- abind(mc_array, sapply(mc.d, function(x) x[, i]), along=3)
}

indexa.obj <- new("indexa.mc", reads=Y, dirichlet.samples=dir_array,
                  log2.norm.samples=mc_array,
                  scales=matrix(0, dim(dir_array)[3], ncol(Y)),
                  denom=NULL)

aldex2_res <- aldex.ttest(clr.obj, verbose=F, paired.test=F)
indexa_res <- indexa.test(indexa.obj, conds, c(0, 0), equal.var=F)
indexa_res_i <- indexa.test(indexa.obj, conds, c(-0.2, 0.2), equal.var=F)
indexa_res_c <- indexa.test(indexa.obj, conds, c(0, true_scale), equal.var=F)

test_that("indexa.test function comparison to ALDEx2 sanity checks", {

    expect_error(
        indexa.test(indexa.obj, conds, c(-1, 0, 1), equal.var=F),
        "'epsilon.interval' must contain exactly 2 elements." 
    )
    
    expect_error(
        indexa.test(indexa.obj, conds, c(1, -1), equal.var=F),
        "The first element of 'epsilon.interval' cannot be greater than the second." 
    )

    expect_error(
        indexa.test(indexa.obj, c(0,1,0,1), c(-1, 1), equal.var=F),
        "'conds' must be a vector with the same length as the number of samples in the sequence count reads matrix." 
    )

    expect_error(
        indexa.test(indexa.obj, cbind(c(0,1), c(1,0)), c(-1, 1), equal.var=F),
        "'conds' must be a vector with the same length as the number of samples in the sequence count reads matrix." 
    )

     expect_error(
        indexa.test(indexa.obj, rep("A", N), c(-1, 1), equal.var=F),
        "parameter must be either a binary" 
    )
    
    expect_equal(
        unname(indexa_res[,"ctt.pval"]),
        aldex2_res[,"we.ep"]
    )
    
    expect_equal(
        unname(indexa_res[,"ctt.pval.BH.adj"]),
        aldex2_res[,"we.eBH"]
    )

    expect_true(
        all(unname(indexa_res_i[,"ctt.pval.BH.adj"])>aldex2_res[,"we.eBH"])
    )

    expect_true(
        indexa_res_c[86,"ctt.pval.BH.adj"]>0.05
    )
    
    expect_true(
        indexa_res_c[70,"ctt.pval.BH.adj"]<=0.05
    )

})

