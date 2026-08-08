# Segmented, stepmented and both: the working blocks, and the iteration
# that estimates the break-points.

# One step of the algorithm: refresh the block at the current
# coefficients, fit the working linear model, take the new coefficients.
# That loop IS Muggeo's iteration for the continuous case and Fasola's
# for the discontinuous one -- the term supplies the block, nothing here
# knows which construction it is looking at.
seg_iterate <- function(built, y, iters = 40, damp = 1) {
  b <- built@blueprint$coef
  for (it in seq_len(iters)) {
    cur <- term_refresh(built, b)
    X <- term_matrix(cur)
    # for the continuous case the block is a Jacobian, so the fit gives
    # an increment; for the jump it gives the coefficients themselves
    if (cur@kind == "seg") {
      r <- y - term_value(cur)
      step <- qr.solve(crossprod(X) + 1e-10 * diag(ncol(X)), crossprod(X, r))
      b <- b + damp * as.numeric(step)
    } else {
      nb <- qr.solve(crossprod(X) + 1e-10 * diag(ncol(X)), crossprod(X, y))
      b <- as.numeric(b + damp * (nb - b))
    }
  }
  b
}

test_that("seg recovers a change of slope", {
  set.seed(1)
  n <- 300
  dd <- data.frame(x = sort(runif(n, 0, 10)))
  dd$y <- 1 + 0.5 * dd$x + 2 * pmax(dd$x - 6, 0) + rnorm(n, sd = 0.3)

  built <- term_build(seg(x, psi = 4), dd)
  expect_identical(term_coef_names(built), c("seg.lin", "seg.delta1",
                                             "seg.psi1"))
  b <- seg_iterate(built, dd$y - 1)
  expect_equal(seg_psi(built, b), 6, tolerance = 0.15)
  expect_equal(b[2], 2, tolerance = 0.2)   # the change of slope
  expect_equal(b[1], 0.5, tolerance = 0.2) # the slope before it
})

test_that("the seg block is the jacobian of the contribution", {
  set.seed(2)
  dd <- data.frame(x = sort(runif(60, 0, 10)))
  built <- term_build(seg(x, psi = 5), dd)
  b <- c(0.4, 1.7, 5.5)

  f <- function(v) v[1] * dd$x + v[2] * pmax(dd$x - v[3], 0)
  cur <- term_refresh(built, b)
  expect_equal(term_value(cur), f(b), tolerance = 1e-12)

  # A central difference straddling the kink smears it, so at a point
  # within a difference step of the break-point the REFERENCE is wrong
  # and not the value; those rows are dropped rather than loosened.
  ok <- abs(dd$x - b[3]) > 1e-3
  expect_true(sum(ok) > length(ok) - 3)
  expect_equal(unname(term_matrix(cur))[ok, ],
               numDeriv::jacobian(f, b)[ok, ], tolerance = 1e-5)
})

test_that("the jump identity is exact away from the break-point", {
  # kappa * Z + g * W IS the step, exactly, when the weight is evaluated
  # at the same psi the coefficients imply
  set.seed(3)
  dd <- data.frame(x = sort(runif(80, 0, 10)))
  built <- term_build(jump(x, psi = 4), dd)
  kappa <- 2.5
  psi <- 4
  b <- c(0, kappa, -kappa * psi)
  cur <- term_refresh(built, b)
  X <- term_matrix(cur)

  lin <- as.numeric(X %*% b)
  # exact outside the band, where the weight is not capped; inside it the
  # step is deliberately a ramp
  band <- 0.02 * diff(range(dd$x))
  away <- abs(dd$x - psi) > band
  expect_true(sum(away) > 0.8 * length(away))
  expect_equal(lin[away], (kappa * (dd$x > psi))[away], tolerance = 1e-8)
  # and inside the band it is a ramp, not the step
  near <- !away
  if (any(near)) {
    expect_false(isTRUE(all.equal(lin[near], (kappa * (dd$x > psi))[near])))
  }
  # and term_value reports the true step, not the linearization
  expect_equal(term_value(cur), kappa * (dd$x > psi), tolerance = 1e-12)
  expect_equal(seg_psi(cur), psi, tolerance = 1e-12)
})

test_that("jump recovers a discontinuity, from a start well away from it", {
  set.seed(4)
  n <- 400
  dd <- data.frame(x = sort(runif(n, 0, 10)))
  dd$y <- 3 * (dd$x > 6.5) + rnorm(n, sd = 0.4)

  built <- term_build(jump(x, psi = 3, linear = FALSE), dd)
  expect_identical(term_coef_names(built),
                   c("jump.kappa1", "jump.g1"))
  b <- seg_iterate(built, dd$y, iters = 60, damp = 0.5)
  expect_equal(seg_psi(built, b), 6.5, tolerance = 0.2)
  expect_equal(b[1], 3, tolerance = 0.3)   # the size of the jump
})

test_that("jseg recovers a jump and a change of slope at the same point", {
  set.seed(5)
  n <- 500
  dd <- data.frame(x = sort(runif(n, 0, 10)))
  dd$y <- 0.3 * dd$x + 1.5 * pmax(dd$x - 5, 0) + 2 * (dd$x > 5) +
    rnorm(n, sd = 0.3)

  built <- term_build(jseg(x, psi = 3), dd)
  expect_identical(term_coef_names(built),
                   c("jseg.lin", "jseg.delta1", "jseg.kappa1", "jseg.g1"))
  b <- seg_iterate(built, dd$y, iters = 120, damp = 0.2)
  expect_equal(seg_psi(built, b), 5, tolerance = 0.2)
  expect_equal(b[3], 2, tolerance = 0.4)    # the jump
  expect_equal(b[2], 1.5, tolerance = 0.4)  # the change of slope
  expect_equal(b[1], 0.3, tolerance = 0.2)  # the slope before it
})

test_that("a run can end at the limit, and that is the signal it failed", {
  # The objective has local optima in the break-point. Taken too fast the
  # same problem walks out of the data and stops where it is held, with a
  # residual sum of squares tens of times worse -- which is why the
  # segmented literature restarts from several positions, and why
  # multistart() belongs around this iteration.
  set.seed(5)
  n <- 500
  dd <- data.frame(x = sort(runif(n, 0, 10)))
  dd$y <- 0.3 * dd$x + 1.5 * pmax(dd$x - 5, 0) + 2 * (dd$x > 5) +
    rnorm(n, sd = 0.3)
  built <- term_build(jseg(x, psi = 3), dd)

  rss <- function(b) sum((dd$y - term_value(built, b))^2)
  good <- seg_iterate(built, dd$y, iters = 120, damp = 0.2)
  bad <- seg_iterate(built, dd$y, iters = 120, damp = 0.4)

  lim <- as.numeric(stats::quantile(dd$x, c(0.05, 0.95), names = FALSE))
  expect_equal(seg_psi(built, bad), lim[1], tolerance = 1e-8)
  expect_gt(rss(bad), 10 * rss(good))
})

test_that("the compiled block agrees with the R twin", {
  # Same operations in the same order, so the two are bit for bit; a
  # tolerance here would hide a reordering.
  set.seed(11)
  code <- c(seg = 0L, jump = 1L, jseg = 2L)
  for (n in c(37L, 1000L)) {
    xv <- sort(runif(n, 0, 10))
    lim <- as.numeric(stats::quantile(xv, c(0.05, 0.95), names = FALSE))
    floor_w <- 2 * 0.02 * diff(range(xv))
    for (kind in c("seg", "jump", "jseg")) {
      for (npsi in c(1L, 3L)) {
        for (linear in c(TRUE, FALSE)) {
          npar <- (if (linear) 1L else 0L) +
            (if (kind == "jseg") 3L else 2L) * npsi
          cf <- runif(npar, 0.5, 2) * sample(c(-1, 1), npar, TRUE)
          a <- modelterms7:::.seg_block_r(kind, xv, cf, npsi, linear,
                                          floor_w, lim)
          b <- modelterms7:::.seg_block(kind, xv, cf, npsi, linear,
                                        floor_w, lim)
          info <- paste(kind, npsi, linear, n)
          expect_identical(b$X, a$X, info = info)
          expect_identical(b$value, a$value, info = info)
          expect_identical(b$psi, a$psi, info = info)
        }
      }
    }
  }

  # A level of `by` can be empty on a subset, which is how check_term
  # reaches this: taking the address of the first element of a matrix
  # with no rows is out of bounds.
  e <- modelterms7:::.seg_block("jseg", numeric(0), c(1, 1, 2, -4), 1L,
                                TRUE, 0.2, c(0, 10))
  expect_identical(dim(e$X), c(0L, 4L))
  expect_identical(e$value, numeric(0))
})

test_that("several break-points are carried and reported in order", {
  set.seed(6)
  dd <- data.frame(x = sort(runif(200, 0, 10)))
  built <- term_build(seg(x, npsi = 3), dd)
  expect_identical(term_coef_names(built),
                   c("seg.lin", "seg.delta1", "seg.delta2", "seg.delta3",
                     "seg.psi1", "seg.psi2", "seg.psi3"))
  # the default starting positions are the interior quantiles
  expect_equal(seg_psi(built),
               as.numeric(stats::quantile(dd$x, c(0.25, 0.5, 0.75),
                                          names = FALSE)),
               tolerance = 1e-10)

  bj <- term_build(jump(x, npsi = 2), dd)
  expect_identical(term_coef_names(bj),
                   c("jump.lin", "jump.kappa1", "jump.kappa2",
                     "jump.g1", "jump.g2"))
  expect_identical(length(seg_psi(bj)), 2L)
})

test_that("by gives an independent break-point per level", {
  set.seed(7)
  n <- 400
  dd <- data.frame(x = sort(runif(n, 0, 10)),
                   g = factor(sample(c("a", "b"), n, TRUE)))
  built <- term_build(seg(x, by = g), dd)
  expect_identical(term_npar(built), 6L)
  expect_true(all(grepl("^seg\\.(a|b)\\.", term_coef_names(built))))
  expect_identical(length(seg_psi(built)), 2L)

  # a level's columns vanish off its own rows
  X <- term_matrix(built)
  expect_true(all(X[dd$g != "a", 1:3] == 0))
  expect_true(all(X[dd$g != "b", 4:6] == 0))

  res <- check_term(seg(x, by = g), dd, verbose = FALSE)
  expect_true(all(res$status == "OK"),
              info = paste(res$check[res$status != "OK"], collapse = ", "))
})

test_that("a penalty reaches the changes and nothing else", {
  set.seed(8)
  dd <- data.frame(x = sort(runif(120, 0, 10)))
  built <- term_build(seg(x, npsi = 2, penalty = "lasso"), dd)
  pen <- term_penalty(built)
  expect_identical(pen@params, "lambda")
  expect_identical(pen@n_coef, term_npar(built))
  expect_false(term_smooth(built))

  # the map selects the two slope changes out of the five coefficients
  D <- pen@map
  expect_identical(dim(D), c(2L, 5L))
  expect_identical(which(D[1, ] == 1), 2L)
  expect_identical(which(D[2, ] == 1), 3L)

  # ridge is the smooth alternative and reaches the same coefficients
  br <- term_build(seg(x, npsi = 2, penalty = "ridge"), dd)
  expect_true(term_smooth(br))
  expect_identical(dim(term_penalty(br)@map), c(2L, 5L))
  # with no penalty there is none
  expect_null(term_penalty(term_build(seg(x), dd)))
})

test_that("prediction reapplies the break-points and the terms validate", {
  set.seed(9)
  dd <- data.frame(x = sort(runif(150, 0, 10)),
                   g = factor(rep(c("a", "b"), length.out = 150)))
  for (spec in list(seg(x), jump(x), jseg(x), seg(x, npsi = 2),
                    jump(x, linear = FALSE))) {
    res <- check_term(spec, dd, verbose = FALSE)
    expect_true(all(res$status == "OK"),
                info = paste(res$check[res$status != "OK"], collapse = ", "))
  }
  built <- term_build(seg(x, psi = 5), dd)
  nd <- data.frame(x = c(1, 5.5, 9), g = factor(c("a", "b", "a"),
                                                levels = c("a", "b")))
  P <- term_predict(built, nd)
  expect_identical(dim(P), c(3L, 3L))
  expect_equal(P[, 2], pmax(nd$x - 5, 0), ignore_attr = TRUE)
})

test_that("the constructors reject what they cannot honour", {
  expect_error(seg(x, npsi = 0), "at least 1")
  expect_error(seg(x, npsi = 2, psi = 1), "one starting position")
  expect_error(jump(x, linear = NA), "TRUE or FALSE")
  expect_error(jseg(x, label = ""), "non-empty")
  expect_error(seg(x, penalty = "nope"), "should be one of")
  dd <- data.frame(x = rep(1, 10))
  expect_error(term_build(seg(x), dd), "must vary")
  expect_error(seg_psi(seg(x)), "not been built")
})

test_that("the interpreter routes the three and print reports the points", {
  set.seed(10)
  dd <- data.frame(y = rnorm(60), x = sort(runif(60, 0, 10)))
  out <- interpret_formula(y ~ seg(x) + jump(x, label = "j"), dd)
  expect_named(out$terms, c("linpar", "seg(x)", 'jump(x, label = "j")'))
  expect_output(print(out$terms[["seg(x)"]]), "specification")
  expect_output(print(term_build(out$terms[["seg(x)"]], dd)), "break-point")
})
