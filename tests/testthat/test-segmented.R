# Segmented, stepmented and both: the working blocks, and the iteration
# that estimates the break-points.

# One step of the algorithm: refresh the block at the current
# coefficients, fit the working linear model, take the new coefficients.
# That loop IS Muggeo's iteration for the continuous case and Fasola's
# for the discontinuous one -- the term supplies the block, nothing here
# knows which construction it is looking at.
seg_iterate <- function(built, y, iters = 40, damp = 1, stop_early = TRUE) {
  b <- built@blueprint$coef
  cur <- built
  for (it in seq_len(iters)) {
    # the term is chained rather than rebuilt, because the scaling
    # factor of a discontinuous term is a state of the iteration
    cur <- term_refresh(cur, b)
    X <- term_matrix(cur)
    # for the continuous case the block is a Jacobian, so the fit gives
    # an increment; for the jump it gives the coefficients themselves
    # a QR of the design, not the normal equations: the working block of
    # a discontinuous term is conditioned to eps^-1/2 by construction and
    # crossprod(X) squares that
    if (cur@kind == "seg") {
      b <- b + damp * as.numeric(qr.coef(qr(X), y - term_value(cur)))
    } else {
      b <- as.numeric(qr.coef(qr(X), y))
    }
    if (stop_early && seg_converged(cur)) break
  }
  attr(b, "iters") <- it
  attr(b, "converged") <- seg_converged(cur)
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

test_that("the jump identity is exact at every observation", {
  # kappa * Z + g * W IS the step, with no exceptions: the rescaling
  # moves the observations off the break-point, so the weight is finite
  # everywhere and the identity holds row for row. Capping the weight
  # instead would have made this true only away from the break-point.
  set.seed(3)
  dd <- data.frame(x = sort(runif(80, 0, 10)))
  built <- term_build(jump(x, psi = 4), dd)
  kappa <- 2.5
  psi <- 4
  b <- c(0, kappa, -kappa * psi)
  cur <- term_refresh(built, b)
  X <- term_matrix(cur)

  expect_equal(as.numeric(X %*% b), kappa * (dd$x > psi), tolerance = 1e-12)
  expect_true(all(is.finite(X)))
  # and term_value reports the true step, not the linearization
  expect_equal(term_value(cur), kappa * (dd$x > psi), tolerance = 1e-12)
  expect_equal(seg_psi(cur), psi, tolerance = 1e-12)

  # an observation sitting exactly on the break-point is the case a cap
  # exists for, and the rescaling handles it without one
  d2 <- data.frame(x = c(sort(runif(40, 0, 10)), 4))
  c2 <- term_refresh(term_build(jump(x, psi = 4), d2), b)
  expect_true(all(is.finite(term_matrix(c2))))
  expect_equal(as.numeric(term_matrix(c2) %*% b), kappa * (d2$x > psi),
               tolerance = 1e-12)
})

test_that("jump recovers a discontinuity, from a start well away from it", {
  set.seed(4)
  n <- 400
  dd <- data.frame(x = sort(runif(n, 0, 10)))
  dd$y <- 3 * (dd$x > 6.5) + rnorm(n, sd = 0.4)

  built <- term_build(jump(x, psi = 3, linear = FALSE), dd)
  expect_identical(term_coef_names(built),
                   c("jump.kappa1", "jump.g1"))
  b <- seg_iterate(built, dd$y, iters = 60)
  expect_equal(seg_psi(built, b), 6.5, tolerance = 0.2)
  expect_equal(b[1], 3, tolerance = 0.3)   # the size of the jump
})

test_that("jseg recovers a jump and a change of slope at the same point", {
  set.seed(5)
  n <- 500
  dd <- data.frame(x = sort(runif(n, 0, 10)))
  dd$y <- 0.3 * dd$x + 1.5 * pmax(dd$x - 5, 0) + 2 * (dd$x > 5) +
    rnorm(n, sd = 0.3)

  built <- term_build(jseg(x, psi = 6), dd)
  expect_identical(term_coef_names(built),
                   c("jseg.lin", "jseg.delta1", "jseg.kappa1", "jseg.g1"))
  b <- seg_iterate(built, dd$y, iters = 120)
  expect_equal(seg_psi(built, b), 5, tolerance = 0.2)
  expect_equal(b[3], 2, tolerance = 0.4)    # the jump
  expect_equal(b[2], 1.5, tolerance = 0.4)  # the change of slope
  expect_equal(b[1], 0.3, tolerance = 0.2)  # the slope before it
})

test_that("a distant start reaches a real local optimum, not a failure", {
  # The profile objective of a joint term is riddled with local minima,
  # and from outside the basin the iteration converges to one of them
  # rather than diverging. That is worth telling apart from a broken
  # run, so the test profiles the EXACT objective at the point reached
  # and asserts it really is a minimum of it -- which is the evidence
  # that the algorithm is right and the starting position is not. It is
  # also the argument for multistart() around this iteration.
  set.seed(5)
  n <- 500
  dd <- data.frame(x = sort(runif(n, 0, 10)))
  dd$y <- 0.3 * dd$x + 1.5 * pmax(dd$x - 5, 0) + 2 * (dd$x > 5) +
    rnorm(n, sd = 0.3)

  # the exact profile: psi held fixed, the three linear coefficients out
  prof <- function(psi) {
    Z <- cbind(dd$x, pmax(dd$x - psi, 0), as.numeric(dd$x > psi))
    sum(qr.resid(qr(Z), dd$y)^2)
  }

  near <- seg_iterate(term_build(jseg(x, psi = 6), dd), dd$y, iters = 200)
  far <- seg_iterate(term_build(jseg(x, psi = 3), dd), dd$y, iters = 200)
  p_near <- seg_psi(term_build(jseg(x, psi = 6), dd), near)
  p_far <- seg_psi(term_build(jseg(x, psi = 3), dd), far)

  expect_equal(p_near, 5, tolerance = 0.2)
  expect_gt(abs(p_far - 5), 0.5)
  # the distant run STOPPED, on the scaling schedule's own rule, rather
  # than running out of iterations: what it reached is an optimum of the
  # problem and not the point a budget happened to leave it at
  expect_true(attr(far, "converged"))
  expect_lt(attr(far, "iters"), 200)
  # and it is a worse one, on a plateau of the exact profile: the
  # objective barely moves for a third of a unit either side, which is
  # the shape Fasola et al. display and the reason spurious solutions
  # are so easy to reach
  expect_gt(prof(p_far), 1.5 * prof(p_near))
  swing <- function(p) {
    v <- vapply(p + c(-0.3, -0.1, 0.1, 0.3), prof, numeric(1))
    max(abs(v / prof(p) - 1))
  }
  expect_gt(swing(p_near), 10 * swing(p_far))
})

test_that("the compiled block agrees with the R twin", {
  # The same operations in the same order, so the two agree to one
  # rounding and no more: the tolerance is tight enough that a
  # reordering would fail it. They are NOT bit for bit everywhere --
  # a compiler free to contract `v += d * u` into a fused multiply-add
  # drops the intermediate rounding, which macOS does and this machine
  # does not.
  set.seed(11)
  code <- c(seg = 0L, jump = 1L, jseg = 2L)
  for (n in c(37L, 1000L)) {
    xv <- sort(runif(n, 0, 10))
    lim <- as.numeric(stats::quantile(xv, c(0.05, 0.95), names = FALSE))
    rr <- range(xv)
    cs <- rep(0.05, 3L)
    for (kind in c("seg", "jump", "jseg")) {
      for (npsi in c(1L, 3L)) {
        for (linear in c(TRUE, FALSE)) {
          npar <- (if (linear) 1L else 0L) +
            (if (kind == "jseg") 3L else 2L) * npsi
          cf <- runif(npar, 0.5, 2) * sample(c(-1, 1), npar, TRUE)
          a <- modelterms7:::.seg_block_r(kind, xv, cf, npsi, linear,
                                          cs[seq_len(npsi)], rr[1], rr[2],
                                          lim)
          b <- modelterms7:::.seg_block(kind, xv, cf, npsi, linear,
                                        cs[seq_len(npsi)], rr[1], rr[2],
                                        lim)
          info <- paste(kind, npsi, linear, n)
          expect_equal(b$X, a$X, tolerance = 1e-15, info = info)
          expect_equal(b$value, a$value, tolerance = 1e-15, info = info)
          expect_identical(b$psi, a$psi, info = info)
        }
      }
    }
  }

  # A level of `by` can be empty on a subset, which is how check_term
  # reaches this: taking the address of the first element of a matrix
  # with no rows is out of bounds.
  e <- modelterms7:::.seg_block("jseg", numeric(0), c(1, 1, 2, -4), 1L,
                                TRUE, 0.05, 0, 10, c(0, 10))
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
  # the penalty is declared over the coefficients it covers, and is not
  # attached to the whole block: term_penalty() answers for a penalty over
  # all of it, and there is none
  expect_null(term_penalty(built))
  ent <- term_penalties(built)
  expect_length(ent, 1L)
  expect_identical(ent[[1L]]$name, "delta")
  expect_identical(ent[[1L]]$penalty@params, "lambda")
  expect_false(term_smooth(built))

  # the two slope changes out of the five coefficients, by position
  expect_identical(ent[[1L]]$index, 2:3)
  expect_identical(term_coef_names(built)[ent[[1L]]$index],
                   c("seg.delta1", "seg.delta2"))
  expect_identical(ent[[1L]]$penalty@n_coef, 2L)

  # named as coordinates the map is the identity, which is what keeps the
  # proximal operator available: a selection map is the generalized lasso
  expect_null(ent[[1L]]$penalty@map)
  expect_true(penalties7::has_prox(ent[[1L]]$penalty))

  # ridge is the smooth alternative and reaches the same coefficients
  br <- term_build(seg(x, npsi = 2, penalty = "ridge"), dd)
  expect_true(term_smooth(br))
  expect_identical(term_penalties(br)[[1L]]$index, 2:3)
  # with no penalty there is none
  expect_length(term_penalties(term_build(seg(x), dd)), 0L)
  expect_null(term_penalty(term_build(seg(x), dd)))

  # a specification has nothing to index yet and reports no penalty, as an
  # unbuilt ridge() does, rather than raising
  expect_length(term_penalties(seg(x, penalty = "lasso")), 0L)
  expect_true(term_smooth(seg(x, penalty = "lasso")))
})

test_that("a joint term penalizes its two kinds of change separately", {
  set.seed(81)
  dd <- data.frame(x = sort(runif(120, 0, 10)))
  # a slope change and a jump are not comparable and cannot share a
  # hyperparameter, so jseg declares two penalties
  built <- term_build(jseg(x, penalty = "lasso"), dd)
  ent <- term_penalties(built)
  expect_length(ent, 2L)
  expect_identical(vapply(ent, function(e) e$name, character(1)),
                   c("delta", "kappa"))
  cn <- term_coef_names(built)
  expect_identical(cn[ent[[1L]]$index], "jseg.delta1")
  expect_identical(cn[ent[[2L]]$index], "jseg.kappa1")
  # the break-point pair g is never penalized: psi is read off it
  expect_false(any(grepl("\\.g", cn[unlist(lapply(ent, function(e) e$index))])))

  # a pure jump has only the one kind
  bj <- term_build(jump(x, npsi = 2, penalty = "lasso"), dd)
  ej <- term_penalties(bj)
  expect_length(ej, 1L)
  expect_identical(ej[[1L]]$name, "kappa")
  expect_identical(term_coef_names(bj)[ej[[1L]]$index],
                   c("jump.kappa1", "jump.kappa2"))
})

test_that("a penalty covers the levels of by under one hyperparameter", {
  set.seed(82)
  dd <- data.frame(x = sort(runif(160, 0, 10)),
                   g = factor(rep(c("a", "b"), length.out = 160)))
  built <- term_build(seg(x, by = g, penalty = "lasso"), dd)
  ent <- term_penalties(built)
  expect_length(ent, 1L)
  expect_identical(term_coef_names(built)[ent[[1L]]$index],
                   c("seg.a.delta1", "seg.b.delta1"))
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

test_that("seg_start finds a start every construction converges from", {
  # The grid initialization of Fasola et al. is what turns the basin
  # into a non-issue: measured over eight samples on a joint jump and
  # change of slope, a single start recovers the break-point in between
  # none and half of them depending on where it is put, and the grid in
  # all of them.
  mk <- function(seed, kind) {
    set.seed(seed)
    d <- data.frame(x = sort(runif(400, 0, 10)))
    d$y <- 0.3 * d$x +
      (if (kind != "jump") 1.5 * pmax(d$x - 5, 0) else 0) +
      (if (kind != "seg") 2 * (d$x > 5) else 0) + rnorm(400, sd = 0.3)
    d
  }
  for (kind in c("seg", "jump", "jseg")) {
    hits <- vapply(1:4, function(sd) {
      d <- mk(sd, kind)
      spec <- switch(kind, seg = seg(x), jump = jump(x), jseg = jseg(x))
      st <- seg_start(spec, d, d$y)
      b <- seg_iterate(term_build(st, d), d$y, iters = 200)
      abs(seg_psi(term_build(st, d), b) - 5) < 0.3
    }, logical(1))
    expect_true(all(hits), info = kind)
  }
})

test_that("seg_start scores a grid and handles several points and by", {
  set.seed(21)
  dd <- data.frame(x = sort(runif(600, 0, 10)))
  dd$y <- 0.3 * dd$x + 2 * (dd$x > 3.5) + 2 * (dd$x > 6.5) +
    rnorm(600, sd = 0.3)

  st <- seg_start(jump(x, npsi = 2), dd, dd$y)
  expect_length(st@spec$psi, 2L)
  expect_true(all(abs(sort(st@spec$psi) - c(3.5, 6.5)) < 0.6))
  # the grid lies inside the interval a break-point is held in
  lim <- as.numeric(stats::quantile(dd$x, c(0.05, 0.95), names = FALSE))
  expect_true(all(st@spec$psi >= lim[1] & st@spec$psi <= lim[2]))
  # and the iteration started there reaches the two break-points
  b <- seg_iterate(term_build(st, dd), dd$y, iters = 200)
  expect_equal(sort(seg_psi(term_build(st, dd), b)), c(3.5, 6.5),
               tolerance = 0.15)

  # one set of positions is carried when by splits the rows
  set.seed(22)
  d2 <- data.frame(x = sort(runif(400, 0, 10)),
                   g = factor(sample(c("a", "b"), 400, TRUE)))
  d2$y <- 0.3 * d2$x + ifelse(d2$g == "a", 2 * (d2$x > 4), 2 * (d2$x > 7)) +
    rnorm(400, sd = 0.3)
  s2 <- seg_start(jump(x, by = g), d2, d2$y)
  expect_length(s2@spec$psi, 1L)
  b2 <- seg_iterate(term_build(s2, d2), d2$y, iters = 200)
  expect_equal(seg_psi(term_build(s2, d2), b2), c(4, 7), tolerance = 0.3)
})

test_that("seg_start rejects what it cannot use", {
  dd <- data.frame(x = sort(runif(50, 0, 10)), y = rnorm(50))
  expect_error(seg_start(linpar(~x), dd, dd$y), "segmented term")
  expect_error(seg_start(term_build(seg(x), dd), dd, dd$y), "unbuilt")
  expect_error(seg_start(seg(x, npsi = 3), dd, dd$y, k = 3), "at least 4")
  expect_error(seg_start(seg(x), dd, dd$y[1:10]), "one value per row")
})

test_that("the interpreter routes the three and print reports the points", {
  set.seed(10)
  dd <- data.frame(y = rnorm(60), x = sort(runif(60, 0, 10)))
  out <- interpret_formula(y ~ seg(x) + jump(x, label = "j"), dd)
  expect_named(out$terms, c("linpar", "seg(x)", 'jump(x, label = "j")'))
  expect_output(print(out$terms[["seg(x)"]]), "specification")
  expect_output(print(term_build(out$terms[["seg(x)"]], dd)), "break-point")
})

test_that("term_value on other rows is the contribution, not the block", {
  set.seed(14)
  dd <- data.frame(x = sort(runif(200, 0, 10)))
  dd$y <- 1 + 0.5 * dd$x + 2 * pmax(dd$x - 6, 0) + rnorm(200, sd = 0.3)
  built <- term_build(seg(x, psi = 6), dd)
  cf <- c(0.5, 2, 6)
  nd <- data.frame(x = c(1, 5, 6.5, 9))

  # the segmented function itself, which is continuous at the break-point
  expect_equal(term_value(built, coef = cf, newdata = nd),
               cf[1L] * nd$x + cf[2L] * pmax(nd$x - cf[3L], 0),
               tolerance = 1e-12)
  # the block times the coefficients is the LINEARIZATION and is not that:
  # its third column is the Jacobian in the break-point, so it carries a
  # step of -delta*psi where the construction is continuous
  lin <- as.numeric(term_predict(built, nd) %*% cf)
  expect_gt(max(abs(lin - term_value(built, coef = cf, newdata = nd))), 1)

  # a discontinuous term has no such gap: there the block times the
  # coefficients IS the contribution, the break-point being read off them
  bj <- term_build(jump(x, psi = 6), dd)
  cj <- bj@blueprint$coef
  expect_equal(as.numeric(term_predict(bj, nd) %*% cj),
               term_value(bj, coef = cj, newdata = nd), tolerance = 1e-10)
})
