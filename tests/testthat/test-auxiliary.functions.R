test_that("aux function tests", {
    expect_equal(
        gm.mean(c(0.1,0.6,0.3)),
        prod(c(0.1,0.6,0.3))^(1/3)
    )

    expect_true(
        all(apply(rdirichlet(4, c(1e12, 5e12, 4e12)), 2,
                  function(col) abs(col-c(0.1,0.5,0.4)))<0.0001)
    )

    expect_error(
        gm.mean(c(-1,2,4,5)),
        "gm.mean cannot take negative numbers"
    )

    expect_error(
        rdirichlet(10, c(10,20,-2)),
        "rdirichlet cannot contain negative elements"
    )

})
