# The smoothed break-point terms: with an abs_smoother the break-points are
# ordinary parameters, the block is the true Jacobian, and the developments
# the read-off forbade become legal. The default NULL leaves the working
# construction untouched, which the rest of the suite pins.

sm_data <- function(n = 200, seed = 1) {
  set.seed(seed)
  dd <- data.frame(x = sort(runif(n, 0, 10)))
  dd$y <- 1 + 0.3 * dd$x + 2 * (dd$x > 6) + rnorm(n, sd = 0.3)
  dd
}

test_that("the default is NULL and the working construction is untouched", {
  dd <- sm_data()
  b <- term_build(jump(x), dd)
  expect_null(b@spec$smoothed)
  expect_null(b@blueprint$smooth)
  expect_false(term_jacobian_block(jump(x)))
  # the auxiliary coefficient of the identity of Fasola et al. is still there
  expect_true(any(grepl("\\.g1$", term_coef_names(b))))
})

test_that("a smoothed term is a Jacobian block of any kind", {
  dd <- sm_data()
  sm <- penalties7::smooth_probit(h = 0.3)
  for (ctor in list(seg, jump, jseg)) {
    spec <- ctor(x, psi = 4, smoothed = sm)
    expect_true(term_jacobian_block(spec))
    b <- term_build(spec, dd)
    expect_false(is.null(b@blueprint$smooth))
    # no auxiliary coordinates: the break-point is held directly
    expect_true(any(grepl("\\.psi1$", term_coef_names(b))))
    expect_false(any(grepl("\\.g1$", term_coef_names(b))))
    # the block is the derivative of the contribution: every column against
    # a central difference of term_value in that coefficient
    cf <- b@blueprint$coef
    X <- as.matrix(term_matrix(term_refresh(b, cf)))
    for (j in seq_along(cf)) {
      h <- 1e-6 * max(1, abs(cf[j]))
      up <- cf; up[j] <- cf[j] + h
      dn <- cf; dn[j] <- cf[j] - h
      fd <- (term_value(b, up) - term_value(b, dn)) / (2 * h)
      expect_lt(max(abs(X[, j] - fd)), 1e-5 * max(1, max(abs(fd))))
    }
    # and the reported position is the coefficient itself
    expect_equal(as.numeric(seg_psi(b, cf)), cf[b@blueprint$index$psi1])
  }
})

test_that("smoothed developments of a break-point are legal, penalized ones
          included", {
  dd <- sm_data()
  dd$g <- factor(rep(letters[1:4], each = 50))
  sm <- penalties7::smooth_probit(h = 0.3)
  # unsmoothed: the documented refusal stands
  expect_error(term_build(jump(x, psi ~ random(~1 | g)), dd), "penalty")
  # smoothed: it builds, and the sub-term's penalty is declared
  b <- term_build(jump(x, psi ~ random(~1 | g), smoothed = sm), dd)
  expect_true(length(term_penalties(b)) >= 1L)
  # a jseg development of psi, refused outright before, builds too
  b2 <- term_build(jseg(x, psi ~ 0 + g, smoothed = sm), dd)
  expect_equal(length(b2@blueprint$index$psi1), 4L)
})

test_that("c0 is ignored with a message on a smoothed discontinuous term", {
  sm <- penalties7::smooth_probit(h = 0.3)
  expect_message(jump(x, c0 = 0.1, smoothed = sm), "ignored")
  expect_message(jseg(x, c0 = 0.1, smoothed = sm), "ignored")
  expect_silent(jump(x, smoothed = sm))
  expect_error(seg(x, smoothed = "probit"), "abs_smoother")
})

test_that("the width is resolved at build from the covariate's spacing", {
  dd <- sm_data()
  gap <- stats::median(diff(sort(unique(dd$x))))
  b <- term_build(jump(x, smoothed = penalties7::smooth_probit()), dd)
  expect_equal(b@blueprint$smooth$width, gap)
  # the hyperbolic's parameter is a squared length
  b2 <- term_build(jump(x, smoothed = penalties7::smooth_hyperbolic()), dd)
  expect_equal(b2@blueprint$smooth$width, gap^2)
  # a width the smoother holds wins
  b3 <- term_build(jump(x, smoothed = penalties7::smooth_probit(h = 0.5)), dd)
  expect_equal(b3@blueprint$smooth$width, 0.5)
  # and one below the derived floor is rejected, naming the bound
  expect_error(
    term_build(jump(x, smoothed = penalties7::smooth_probit(h = 1e-12)), dd),
    "floor")
})

test_that("a per-group width needs a partition development and gives one
          width per group", {
  dd <- sm_data()
  dd$g <- factor(rep(letters[1:4], each = 50))
  sm <- penalties7::smooth_probit(per_group = TRUE)
  expect_error(term_build(jump(x, smoothed = sm), dd), "per-group")
  b <- term_build(jump(x, psi ~ 0 + g, smoothed = sm), dd)
  expect_equal(length(b@blueprint$smooth$w_group), 4L)
  expect_true(all(b@blueprint$smooth$w_group > 0))
})

test_that("the closed second derivatives agree with a brute-force dX/dbeta", {
  dd <- sm_data(120)
  sm <- penalties7::smooth_probit(h = 0.4)
  for (ctor in list(seg, jump, jseg)) {
    b <- term_build(ctor(x, psi = 5, smoothed = sm), dd)
    cf <- b@blueprint$coef
    m <- length(cf)
    n <- nrow(dd)
    D <- array(0, c(n, m, m))
    for (j in seq_len(m)) {
      h <- 1e-5 * max(1, abs(cf[j]))
      up <- cf; up[j] <- cf[j] + h
      dn <- cf; dn[j] <- cf[j] - h
      D[, , j] <- (as.matrix(term_matrix(term_refresh(b, up))) -
                     as.matrix(term_matrix(term_refresh(b, dn)))) / (2 * h)
    }
    set.seed(3)
    A <- matrix(rnorm(n * m), n, m)
    v <- rnorm(m)
    ct <- term_block_contract(b, coef = cf, A = A)
    ct_ref <- vapply(seq_len(m), function(j) sum(A * D[, , j]), numeric(1))
    expect_lt(max(abs(ct - ct_ref)), 1e-4 * max(1, max(abs(ct_ref))))
    dv <- term_block_deriv(b, coef = cf, v = v)
    dv_ref <- matrix(0, n, m)
    for (j in seq_len(m)) dv_ref <- dv_ref + D[, , j] * v[j]
    expect_lt(max(abs(dv - dv_ref)), 1e-4 * max(1, max(abs(dv_ref))))
    # the two generics are adjoint: <A, dX[v]> = <contract(A), v>
    expect_equal(sum(A * dv), sum(ct * v), tolerance = 1e-10)
  }
})

test_that("a smoothed term predicts and values on new data through its own
          blueprint", {
  dd <- sm_data()
  sm <- penalties7::smooth_probit(h = 0.3)
  b <- term_build(jseg(x, psi = 6, smoothed = sm), dd)
  cf <- b@blueprint$coef
  nd <- dd[seq(1, 200, by = 7), , drop = FALSE]
  expect_equal(term_value(b, cf, newdata = nd),
               term_value(b, cf)[seq(1, 200, by = 7)])
  expect_equal(as.matrix(term_predict(b, nd)),
               as.matrix(term_matrix(b))[seq(1, 200, by = 7), ],
               ignore_attr = TRUE)
})

test_that("relocate and the profile carry over to a smoothed term", {
  dd <- sm_data()
  sm <- penalties7::smooth_probit(h = 0.3)
  b <- term_build(jump(x, psi = 3, smoothed = sm), dd)
  b2 <- seg_relocate(b, 6)
  expect_equal(as.numeric(seg_psi(b2)), 6)
  # the exact profile does not depend on the mollifier: polishing walks to
  # the true position
  b3 <- seg_polish(term_build(jump(x, psi = 2, smoothed = sm), dd), dd$y)
  expect_lt(abs(as.numeric(seg_psi(b3)) - 6), 0.3)
})

test_that("a smoothed seg converges to the sharp answer on a smooth truth", {
  set.seed(4)
  dd <- data.frame(x = sort(runif(300, 0, 10)))
  dd$y <- 1 + 0.5 * dd$x + 2 * pmax(dd$x - 6, 0) + rnorm(300, sd = 0.2)
  fit_gn <- function(spec) {
    b <- term_build(spec, dd)
    cf <- b@blueprint$coef
    for (it in 1:60) {
      b <- term_refresh(b, cf)
      X <- cbind(1, as.matrix(term_matrix(b)))
      r <- dd$y - term_value(b, cf)
      inc <- qr.coef(qr(X), r)
      inc[!is.finite(inc)] <- 0
      cf <- cf + inc[-1L]
      if (max(abs(inc)) < 1e-9) break
    }
    seg_psi(b, cf)
  }
  sharp <- fit_gn(seg(x, psi = 5))
  smoothd <- fit_gn(seg(x, psi = 5,
                        smoothed = penalties7::smooth_probit(h = 0.05)))
  expect_lt(abs(as.numeric(sharp) - as.numeric(smoothd)), 0.05)
})
