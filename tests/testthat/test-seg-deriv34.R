# term_block_deriv3() and term_block_deriv4() of a smoothed break-point term
# against a central difference of the order below

seg_d34_data <- function() {
  set.seed(7)
  d <- data.frame(x = sort(runif(60, 0, 10)), h = rbinom(60, 1, 0.4))
  d$y <- 1 + 0.5 * d$x + 2 * pmax(d$x - 5, 0) + rnorm(60, sd = 0.3)
  d
}

seg_d34_check <- function(b, order) {
  m <- length(term_coef_names(b))
  cf <- b@build_info$coef
  set.seed(3)
  dirs <- lapply(1:4, function(i) rnorm(m))
  if (order == 3L) {
    an <- term_block_deriv3(b, coef = cf, v = dirs[[1]], u = dirs[[2]],
                            s = dirs[[3]])
    num <- matrix(numDeriv::jacobian(function(t)
      as.vector(term_block_deriv2(b, coef = cf + t * dirs[[3]], v = dirs[[1]],
                                  u = dirs[[2]])), 0), nrow(an))
    sw <- term_block_deriv3(b, coef = cf, v = dirs[[3]], u = dirs[[1]],
                            s = dirs[[2]])
  } else {
    an <- term_block_deriv4(b, coef = cf, v = dirs[[1]], u = dirs[[2]],
                            s = dirs[[3]], r = dirs[[4]])
    num <- matrix(numDeriv::jacobian(function(t)
      as.vector(term_block_deriv3(b, coef = cf + t * dirs[[4]], v = dirs[[1]],
                                  u = dirs[[2]], s = dirs[[3]])), 0), nrow(an))
    sw <- term_block_deriv4(b, coef = cf, v = dirs[[4]], u = dirs[[3]],
                            s = dirs[[1]], r = dirs[[2]])
  }
  list(err = max(abs(an - num)) / max(abs(num)),
       sym = max(abs(an - sw)) / max(abs(an)), scale = max(abs(an)))
}

test_that("a smoothed seg() block's third and fourth derivatives match", {
  skip_if_not_installed("numDeriv")
  d <- seg_d34_data()
  sm <- numericals7::smooth_probit()
  for (f in list(seg(x, smoothed = sm), seg(x, psi ~ 1 + h, smoothed = sm),
                 seg(x, gamma ~ 1 + h, npsi = 2, smoothed = sm))) {
    b <- term_build(f, d)
    for (o in 3:4) {
      r <- seg_d34_check(b, o)
      expect_gt(r$scale, 0)
      expect_lt(r$err, 1e-7)
      expect_lt(r$sym, 1e-12)
    }
  }
})

test_that("a smoothed jump() and jseg() have the third derivative", {
  skip_if_not_installed("numDeriv")
  d <- seg_d34_data()
  sm <- numericals7::smooth_probit()
  for (f in list(jump(x, smoothed = sm), jseg(x, smoothed = sm))) {
    b <- term_build(f, d)
    r <- seg_d34_check(b, 3L)
    expect_lt(r$err, 1e-7)
    expect_lt(r$sym, 1e-12)
    # the fourth reads the sixth derivative of the smoother
    m <- length(term_coef_names(b))
    expect_error(term_block_deriv4(b, v = rep(1, m), u = rep(1, m),
                                   s = rep(1, m), r = rep(1, m)),
                 "derivative of order 6")
  }
})

test_that("a sharp seg() block has zero third and fourth derivatives", {
  b <- term_build(seg(x), seg_d34_data())
  m <- length(term_coef_names(b))
  expect_true(all(term_block_deriv3(b, v = rep(1, m), u = rep(1, m),
                                    s = rep(1, m)) == 0))
  expect_true(all(term_block_deriv4(b, v = rep(1, m), u = rep(1, m),
                                    s = rep(1, m), r = rep(1, m)) == 0))
})
