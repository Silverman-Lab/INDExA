set.seed(1)
x <- rnorm(10, 0)
y <- rnorm(10, 1)

wilcox_res_1     <- wilcox.test(y, x)
wilcox_int_res_1 <- wilcox.interval.test(y, x, 0, 0)

test_that("wilcox interval function matches results from wilcox.test function", {
    expect_equal(
        round(wilcox_res_1$p.value, 5),
        round(wilcox_int_res_1$upper.pval, 5)
    )

})

