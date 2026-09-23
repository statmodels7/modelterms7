# A break-point sitting on an observation: which side the observation is read
# on, and which coefficients move a kink.

kink_data <- function() {
  d <- data.frame(x = seq(0, 10, length.out = 41))
  d$y <- 1 + 0.5 * d$x + 2 * pmax(d$x - 6, 0)
  d
}

test_that("an observation on the break-point is read on one side whatever its last digits", {
  d <- kink_data()
  b <- term_build(seg(x), d)
  bp <- b@blueprint
  # the observation at x = 6, reached from either side by the amounts the
  # inner fit locates a mode to
  cols <- lapply(c(-2e-11, 0, 2e-11), function(e) {
    blk <- modelterms7:::.seg_block(bp, d$x, c(0.5, 2, 6 + e), bp$cscale)
    as.matrix(blk$X)[, 3L]
  })
  expect_identical(cols[[1]], cols[[2]])
  expect_identical(cols[[2]], cols[[3]])
  # and the side is the INACTIVE one: the observation at 6 has no entry
  expect_identical(cols[[2]][d$x == 6], 0)
  # an ordinary displacement still moves the column across it
  off <- as.matrix(modelterms7:::.seg_block(bp, d$x, c(0.5, 2, 6 - 1e-3),
                                           bp$cscale)$X)[, 3L]
  expect_identical(off[d$x == 6], -2)
})

test_that("the compiled block reads the tie as its R twin does", {
  d <- kink_data()
  for (e in c(-2e-11, 0, 2e-11, 1e-3)) {
    b <- term_build(seg(x), d)
    bp <- b@blueprint
    cf <- c(0.5, 2, 6 + e)
    a <- modelterms7:::.seg_block(bp, d$x, cf, bp$cscale)
    k <- modelterms7:::.seg_block_cpp(bp, d$x, cf, bp$cscale)
    expect_identical(unname(as.matrix(a$X))[, 3L], unname(k$X)[, 3L],
                     info = format(e))
  }
})

test_that("term_kinks() names the coefficients that move a break-point on an observation", {
  d <- kink_data()
  b <- term_build(seg(x), d)
  expect_identical(term_kinks(b, c(0.5, 2, 6)), 3L)
  expect_identical(term_kinks(b, c(0.5, 2, 6 + 2e-11)), 3L)
  expect_identical(term_kinks(b, c(0.5, 2, 6.1)), integer(0))
  # a term with no kink, and the constructions that are not fitted by a line
  # search on the objective, answer nothing
  expect_identical(term_kinks(term_build(linpar(~ x), d)), integer(0))
  expect_identical(term_kinks(term_build(jump(x), d)), integer(0))
  sm <- term_build(seg(x, smoothed = numericals7::smooth_probit()), d)
  expect_identical(term_kinks(sm, c(0.5, 2, 6)), integer(0))
})

test_that("a developed break-point names its intercept and that group's deviation only", {
  set.seed(1)
  dd <- data.frame(x = runif(90, 0, 10), id = factor(rep(1:3, each = 30)))
  dd$y <- 1 + dd$x + rnorm(90)
  b <- term_build(seg(x, psi ~ random(~1 | id)), dd)
  ix <- b@blueprint$index$psi1
  cf <- b@blueprint$coef
  # group 2's break-point placed on its first observation
  cf[ix[1L]] <- dd$x[dd$id == "2"][1L] - cf[ix[3L]]
  expect_identical(term_kinks(b, cf), ix[c(1L, 3L)])
  cf[ix[1L]] <- cf[ix[1L]] + 1e-3
  expect_identical(term_kinks(b, cf), integer(0))
})
