set.seed(1)

t_test_res_eq <- c()
ctt_res_0_eq <- c()
ctt_res_i_eq <- c()
t_test_int_eq <- c()
for(i in 1:10000) {
    n1 <- sample(seq(4, 30, 2), 1)
    n2 <- sample(seq(4, 30, 2), 1)
    lfc <- runif(1, -3, 3)
    x <- rnorm(n1, 0)
    y <- rnorm(n2, lfc)
    epsilon_u <- runif(1, -1, 2)
    epsilon_l <- epsilon_u + runif(1, -1, 0)
    z <- c(x, y)
    conds <- c(rep(0, n1), rep(1, n2))

    t_test_res_eq <- c(t_test_res_eq, t.test(y, x, var.equal=T)$p.value)
    ctt_res_0_eq <- c(ctt_res_0_eq,
                      ttest.interval.fast(y, x, 0, 0, equal.var=T)$p) 
    ctt_res_i_eq <- c(ctt_res_i_eq,
                      ttest.interval.fast(y, x, epsilon_l, epsilon_u, equal.var=T)$p)
    t_test_int_eq <- c(t_test_int_eq,
                       min(2*min(t.test(y, x, mu=epsilon_l, alternative="less", var.equal=T)$p.value,
                                 t.test(y, x, mu=epsilon_u, alternative="greater", var.equal=T)$p.value),
                           1))
}

t_test_res_neq <- c()
ctt_res_0_neq <- c()
ctt_res_i_neq <- c()
t_test_int_neq <- c()
for(i in 1:10000) {
    n1 <- sample(seq(4, 30, 2), 1)
    n2 <- sample(seq(4, 30, 2), 1)
    s1 <- runif(1, 0.1, 2)
    s2 <- runif(1, 0.1, 2)
    lfc <- runif(1, -3, 3)
    x <- rnorm(n1, 0, s1)
    y <- rnorm(n2, lfc, s2)
    epsilon_u <- runif(1, -1, 2)
    epsilon_l <- epsilon_u + runif(1, -1, 0)
    z <- c(x, y)
    conds <- c(rep(0, n1), rep(1, n2))
    
    t_test_res_neq <- c(t_test_res_neq, t.test(y, x, var.equal=F)$p.value)
    ctt_res_0_neq <- c(ctt_res_0_neq, ttest.interval.fast(y, x, 0, 0, equal.var=F)$p)
    ctt_res_i_neq <- c(ctt_res_i_neq, ttest.interval.fast(y, x, epsilon_l, epsilon_u, equal.var=F)$p)
    t_test_int_neq <- c(t_test_int_neq,
                       min(2*min(t.test(y, x, mu=epsilon_l, alternative="less", var.equal=F)$p.value,
                                 t.test(y, x, mu=epsilon_u, alternative="greater", var.equal=F)$p.value),
                           1))
}

test_that("ttest.interval.fast matches results from t.test function", {
    expect_equal(
        ctt_res_0_eq,
        t_test_res_eq
    )

    expect_equal(
        ctt_res_i_eq,
        t_test_int_eq
    )
    
    expect_equal(
        t_test_res_neq,
        ctt_res_0_neq
    )

    expect_equal(
        ctt_res_i_neq,
        t_test_int_neq
    )
    
})

