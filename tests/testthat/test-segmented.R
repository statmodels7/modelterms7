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
  expect_identical(term_coef_names(built), c("seg.beta", "seg.gamma1",
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
  b <- c(kappa, -kappa * psi)
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

  built <- term_build(jump(x, psi = 3), dd)
  expect_identical(term_coef_names(built),
                   c("jump.delta1", "jump.g1"))
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
                   c("jseg.beta", "jseg.gamma1", "jseg.delta1", "jseg.g1"))
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
  # does not. The compiled route is taken wherever no coefficient carries
  # a development, which is the ordinary case and the one it was written
  # for; a development goes through the R form.
  set.seed(11)
  for (n in c(37L, 1000L)) {
    dd <- data.frame(x = sort(runif(n, 0, 10)))
    for (kind in c("seg", "jump", "jseg")) {
      for (npsi in c(1L, 3L)) {
        for (linear in c(TRUE, FALSE)) {
          if (kind == "jump" && linear) next
          ctor <- switch(kind, seg = seg, jump = jump, jseg = jseg)
          spec <- if (kind == "jump") ctor(x, npsi = npsi) else
            ctor(x, npsi = npsi, linear = linear)
          built <- term_build(spec, dd)
          bp <- built@blueprint
          npar <- ncol(built@X)
          cf <- runif(npar, 0.5, 2) * sample(c(-1, 1), npar, TRUE)
          a <- modelterms7:::.seg_block(bp, dd$x, cf, bp$cscale)
          b <- modelterms7:::.seg_block_cpp(bp, dd$x, cf, bp$cscale)
          info <- paste(kind, npsi, linear, n)
          expect_equal(unname(as.matrix(a$X)), unname(b$X),
                       tolerance = 1e-15, info = info)
          expect_equal(a$value, b$value, tolerance = 1e-15, info = info)
          expect_equal(a$psi, b$psi, tolerance = 1e-15, info = info)
        }
      }
    }
  }

  # A subset can be empty, which is how check_term reaches this: taking the
  # address of the first element of a matrix with no rows is out of bounds.
  built <- term_build(jseg(x), data.frame(x = sort(runif(50, 0, 10))))
  e <- modelterms7:::.seg_block_cpp(built@blueprint, numeric(0),
                                    c(1, 1, 2, -4), built@blueprint$cscale)
  expect_identical(dim(e$X), c(0L, 4L))
  expect_identical(e$value, numeric(0))
})

test_that("several break-points are carried and reported in order", {
  set.seed(6)
  dd <- data.frame(x = sort(runif(200, 0, 10)))
  built <- term_build(seg(x, npsi = 3), dd)
  expect_identical(term_coef_names(built),
                   c("seg.beta", "seg.gamma1", "seg.gamma2", "seg.gamma3",
                     "seg.psi1", "seg.psi2", "seg.psi3"))
  # the default starting positions are the interior quantiles
  expect_equal(seg_psi(built),
               as.numeric(stats::quantile(dd$x, c(0.25, 0.5, 0.75),
                                          names = FALSE)),
               tolerance = 1e-10)

  bj <- term_build(jump(x, npsi = 2), dd)
  expect_identical(term_coef_names(bj),
                   c("jump.delta1", "jump.delta2", "jump.g1", "jump.g2"))
  expect_identical(length(seg_psi(bj)), 2L)
})

test_that("by = ~0 + g gives an independent set of everything per level", {
  # the shorthand develops EVERY coefficient of the term on the group
  # indicators, which is the per-level model written as a development, so
  # the columns are parameter-major: one per level within each coefficient
  set.seed(7)
  n <- 400
  dd <- data.frame(x = sort(runif(n, 0, 10)),
                   g = factor(sample(c("a", "b"), n, TRUE)))
  built <- term_build(seg(x, by = ~0 + g), dd)
  expect_identical(term_npar(built), 6L)
  expect_identical(term_coef_names(built),
                   c("seg.beta.ga", "seg.beta.gb", "seg.gamma1.ga",
                     "seg.gamma1.gb", "seg.psi1.ga", "seg.psi1.gb"))

  # a level's columns vanish off its own rows
  X <- as.matrix(term_matrix(built))
  expect_true(all(X[dd$g != "a", c(1L, 3L, 5L)] == 0))
  expect_true(all(X[dd$g != "b", c(2L, 4L, 6L)] == 0))

  res <- check_term(seg(x, by = ~0 + g), dd, verbose = FALSE)
  expect_true(all(res$status == "OK"),
              info = paste(res$check[res$status != "OK"], collapse = ", "))

  # the discontinuous constructions reach the same model, the product the
  # break-point is read off collapsing group by group over indicators
  bj <- term_build(jump(x, by = ~0 + g), dd)
  expect_identical(term_coef_names(bj),
                   c("jump.delta1.ga", "jump.delta1.gb", "jump.g1.ga",
                     "jump.g1.gb"))

  # a bare variable is rejected with the formula it stands for
  expect_error(seg(x, by = g), "by = ~0 \\+ g")
})

test_that("a penalty on the changes is the development on a penalized intercept", {
  set.seed(8)
  dd <- data.frame(x = sort(runif(120, 0, 10)))
  # `gamma ~ 0 + lasso(~1)` says the changes themselves carry a lasso. The
  # `0 +` removes the subformula's own unpenalized intercept, which for an
  # intercept-only development would be the same column twice.
  built <- term_build(seg(x, npsi = 2, gamma ~ 0 + lasso(~1)), dd)
  expect_null(term_penalty(built))
  ent <- term_penalties(built)
  expect_length(ent, 1L)
  # one entry over BOTH changes: a subformula shared by every coefficient
  # of a kind is one penalized block under one hyperparameter, which is
  # what selecting how many break-points there are asks for
  expect_identical(ent[[1L]]$penalty@params, "lambda")
  expect_identical(ent[[1L]]$penalty@n_coef, 2L)
  expect_identical(term_coef_names(built)[ent[[1L]]$index],
                   c("seg.gamma1.lasso.(Intercept)",
                     "seg.gamma2.lasso.(Intercept)"))
  expect_false(term_smooth(built))

  # named as coordinates the map is the identity, which is what keeps the
  # proximal operator available: a selection map is the generalized lasso
  expect_null(ent[[1L]]$penalty@map)
  expect_true(penalties7::has_prox(ent[[1L]]$penalty))

  # the break-points are not touched by it
  expect_false(any(grepl("psi", term_coef_names(built)[ent[[1L]]$index])))

  # ridge is the smooth alternative over the same coefficients
  br <- term_build(seg(x, npsi = 2, gamma ~ 0 + ridge(~1)), dd)
  expect_true(term_smooth(br))
  expect_identical(term_penalties(br)[[1L]]$penalty@n_coef, 2L)

  # with no development there is no penalty
  expect_length(term_penalties(term_build(seg(x), dd)), 0L)
  expect_null(term_penalty(term_build(seg(x), dd)))
  # and a specification has nothing to index yet
  expect_length(term_penalties(seg(x, npsi = 2, gamma ~ 0 + lasso(~1))), 0L)

  # an unshared subformula names one coefficient and is not pooled
  one <- term_build(seg(x, npsi = 2, gamma1 ~ 0 + lasso(~1)), dd)
  e1 <- term_penalties(one)
  expect_length(e1, 1L)
  expect_identical(e1[[1L]]$penalty@n_coef, 1L)
  expect_true(startsWith(e1[[1L]]$name, "gamma1::"))
})

test_that("a joint term penalizes its two kinds of change separately", {
  set.seed(81)
  dd <- data.frame(x = sort(runif(120, 0, 10)))
  # a slope change and a change of level are not comparable and cannot
  # share a hyperparameter, so two subformulas give two penalties
  built <- term_build(jseg(x, npsi = 2, gamma ~ 0 + lasso(~1),
                           delta ~ 0 + lasso(~1)), dd)
  ent <- term_penalties(built)
  expect_length(ent, 2L)
  expect_identical(vapply(ent, function(e) e$name, character(1)),
                   c("gamma::lasso(~1)", "delta::lasso(~1)"))
  cn <- term_coef_names(built)
  expect_true(all(grepl("^jseg\\.gamma", cn[ent[[1L]]$index])))
  expect_true(all(grepl("^jseg\\.delta", cn[ent[[2L]]$index])))
  # the auxiliary pair g is never penalized: the break-point is read off it
  expect_false(any(grepl("\\.g[0-9]",
                         cn[unlist(lapply(ent, function(e) e$index))])))
})

test_that("prediction reapplies the break-points and the terms validate", {
  set.seed(9)
  dd <- data.frame(x = sort(runif(150, 0, 10)),
                   g = factor(rep(c("a", "b"), length.out = 150)))
  for (spec in list(seg(x), jump(x), jseg(x), seg(x, npsi = 2),
                    seg(x, linear = FALSE))) {
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
  expect_error(seg(x, linear = NA), "TRUE or FALSE")
  expect_error(jseg(x, label = ""), "non-empty")
  # jump has no linear effect by construction, so it has no such argument
  expect_error(jump(x, linear = FALSE), "not an argument")
  # and none of the three has a penalty argument any more: a penalty on the
  # changes is the development on a penalized intercept
  expect_error(seg(x, penalty = "lasso"), "not an argument")
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
    # a step model has no linear effect, by construction and by
    # definition, so the covariate enters its truth only through the steps
    d$y <- (if (kind != "jump") 0.3 * d$x else 0) +
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
  dd$y <- 2 * (dd$x > 3.5) + 2 * (dd$x > 6.5) + rnorm(600, sd = 0.3)

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
  d2$y <- ifelse(d2$g == "a", 2 * (d2$x > 4), 2 * (d2$x > 7)) +
    rnorm(400, sd = 0.3)
  s2 <- seg_start(jump(x, by = ~0 + g), d2, d2$y)
  expect_length(s2@spec$psi, 1L)
  t2 <- term_build(s2, d2)
  b2 <- seg_iterate(t2, d2$y, iters = 200)
  # a developed break-point carries one position per observation; the
  # positions of a group are its own
  psi2 <- seg_psi(t2, b2)
  got <- vapply(levels(d2$g), function(l) mean(psi2[d2$g == l, 1L]),
                numeric(1))
  expect_equal(unname(got), c(4, 7), tolerance = 0.3)
})

test_that("seg_start rejects what it cannot use", {
  dd <- data.frame(x = sort(runif(50, 0, 10)), y = rnorm(50))
  expect_error(seg_start(linpar(~x), dd, dd$y), "break-point term")
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

test_that("a development of ~1 reproduces the scalar construction", {
  set.seed(21)
  n <- 200
  dd <- data.frame(x = sort(runif(n, 0, 10)))
  dd$y <- 1 + 0.4 * dd$x + 2 * pmax(dd$x - 6, 0) + 1.5 * (dd$x > 6) +
    rnorm(n, sd = 0.3)

  for (ctor in list(seg, jump)) {
    plain <- term_build(ctor(x, psi = 5), dd)
    devd <- term_build(ctor(x, psi ~ 1, psi = 5), dd)
    # same block and same contribution at the starting coefficients, the
    # development being the scalar break-point under another name
    expect_equal(unname(as.matrix(term_matrix(devd))),
                 unname(as.matrix(term_matrix(plain))), tolerance = 1e-10)
    expect_equal(term_value(devd), term_value(plain), tolerance = 1e-10)
  }

  # and the iteration walks the two to the same break-point
  for (ctor in list(seg, jump)) {
    b1 <- seg_iterate(term_build(ctor(x, psi = 4), dd), dd$y, iters = 60)
    t2 <- term_build(ctor(x, psi ~ 1, psi = 4), dd)
    b2 <- seg_iterate(t2, dd$y, iters = 60)
    p1 <- seg_psi(term_build(ctor(x, psi = 4), dd), b1)
    p2 <- seg_psi(t2, b2)
    expect_equal(mean(p2), as.numeric(p1), tolerance = 1e-6)
  }

  # jseg rejects a development: its quadratic reading of the break-point
  # does not split over the columns, and the componentwise reading that
  # remains diverges whenever the jump size passes near zero (measured)
  # jseg takes a development of ~1, which is the scalar case under another
  # name, and rejects one on a design that is not a partition: its
  # quadratic reading of the break-point splits observation by observation
  # only where each observation belongs to one group
  expect_error(term_build(jseg(x, psi ~ id),
                          data.frame(x = sort(runif(60, 0, 10)),
                                     id = factor(rep(c("a", "b"), 30)))),
               "SAME development")
})

test_that("a break-point developed by group finds per-group positions", {
  set.seed(22)
  n <- 240
  truth <- c(a = 4, b = 5.5, c = 6.5)
  dd <- data.frame(x = runif(n, 0, 10),
                   id = factor(rep(c("a", "b", "c"), length.out = n)))
  dd$y <- 0.5 * dd$x + 2.5 * pmax(dd$x - truth[dd$id], 0) +
    rnorm(n, sd = 0.3)

  spec <- seg(x, psi ~ id)
  built <- term_build(seg_start(spec, dd, dd$y), dd)
  # one gamma vector per break-point over the shared design: intercept
  # plus the factor's contrasts, with the slopes shared across groups
  expect_identical(term_coef_names(built),
                   c("seg.beta", "seg.gamma1", "seg.psi1.(Intercept)",
                     "seg.psi1.idb", "seg.psi1.idc"))
  b <- seg_iterate(built, dd$y, iters = 80)
  psi <- seg_psi(built, b)
  expect_true(is.matrix(psi))
  got <- vapply(levels(dd$id), function(l) mean(psi[dd$id == l, 1L]),
                numeric(1))
  expect_equal(unname(got), unname(truth), tolerance = 0.25)
})

test_that("a jump developed by group finds per-group positions", {
  set.seed(23)
  n <- 300
  truth <- c(a = 4, b = 6.5)
  dd <- data.frame(x = runif(n, 0, 10),
                   id = factor(rep(c("a", "b"), length.out = n)))
  dd$y <- 3 * (dd$x > truth[dd$id]) + rnorm(n, sd = 0.4)

  spec <- jump(x, psi ~ id)
  built <- term_build(seg_start(spec, dd, dd$y), dd)
  expect_identical(term_coef_names(built),
                   c("jump.delta1", "jump.g1.(Intercept)", "jump.g1.idb"))
  b <- seg_iterate(built, dd$y, iters = 80)
  psi <- seg_psi(built, b)
  got <- vapply(levels(dd$id), function(l) mean(psi[dd$id == l, 1L]),
                numeric(1))
  expect_equal(unname(got), unname(truth), tolerance = 0.3)
  expect_equal(b[1], 3, tolerance = 0.4)
})

test_that("a penalized development is the random-changepoint model, on seg", {
  set.seed(24)
  n <- 200
  dd <- data.frame(x = runif(n, 0, 10),
                   id = factor(rep(sprintf("g%d", 1:5), length.out = n)))
  dd$y <- 0.5 * dd$x + 2 * pmax(dd$x - 5, 0) + rnorm(n, sd = 0.3)

  built <- term_build(seg(x, psi ~ random(~1 | id)), dd)
  ent <- term_penalties(built)
  expect_length(ent, 1L)
  expect_true(startsWith(ent[[1L]]$name, "psi1::"))
  # the entry covers the indicator columns of the development, in the
  # term's own numbering
  expect_true(all(startsWith(term_coef_names(built)[ent[[1L]]$index],
                             "seg.psi1.")))

  # a discontinuous construction rejects it: the penalty would act on
  # c = -kappa * gamma rather than on the development
  expect_error(term_build(jump(x, psi ~ random(~1 | id)), dd),
               "scaled by the change of level")
})

test_that("the development's syntax is validated where it is written", {
  expect_error(seg(x, psi ~ id, by = ~g), "does not combine")
  expect_error(seg(x, tau ~ id), "not a coefficient of a seg term")
  expect_error(seg(x, psi ~ id, psi ~ z), "developed twice")
  expect_error(seg(x, psi1 ~ id, psi ~ z), "developed twice")
  expect_error(seg(x, ~id), "two-sided formulas")
  # a stem names every coefficient of a kind, a number names one of them,
  # and a term carries only the kinds its construction has
  expect_error(jump(x, gamma ~ id), "not a coefficient of a jump term")
  expect_error(seg(x, delta1 ~ id), "not a coefficient of a seg term")
  expect_error(seg(x, npsi = 2, psi3 ~ id), "not a coefficient")
})

test_that("a developed term predicts by reapplying the sub-design", {
  set.seed(25)
  n <- 120
  dd <- data.frame(x = runif(n, 0, 10),
                   id = factor(rep(c("a", "b"), length.out = n)))
  dd$y <- 0.5 * dd$x + 2 * pmax(dd$x - c(a = 4, b = 6)[dd$id], 0) +
    rnorm(n, sd = 0.3)
  built <- term_build(seg(x, psi ~ id), dd)
  b <- seg_iterate(built, dd$y, iters = 40)
  cur <- term_refresh(built, b)
  keep <- which(dd$id == "a")
  sub <- dd[keep, , drop = FALSE]
  expect_equal(term_value(cur, coef = b, newdata = sub),
               term_value(cur, coef = b)[keep], tolerance = 1e-10)
  expect_equal(unname(as.matrix(term_predict(cur, sub))),
               unname(as.matrix(term_matrix(cur)))[keep, , drop = FALSE],
               tolerance = 1e-10)
})

test_that("a break-point term says where its own coefficients begin", {
  set.seed(31)
  dd <- data.frame(x = sort(runif(120, 0, 10)))

  # zero is degenerate rather than neutral here, so the term answers with
  # the start term_build() computed: unit changes and the break-points at
  # the positions asked for
  b <- term_build(jump(x, npsi = 2, psi = c(3, 7)), dd)
  expect_identical(term_coef_start(b), b@blueprint$coef)
  expect_equal(term_coef_start(b), c(1, 1, -3, -7), tolerance = 1e-12)
  expect_equal(seg_psi(term_refresh(b, term_coef_start(b))), c(3, 7),
               tolerance = 1e-12)

  bs <- term_build(seg(x, psi = 4), dd)
  expect_equal(term_coef_start(bs), c(0, 1, 4), tolerance = 1e-12)

  # an ordinary block wants zeros, which is what the base method gives
  expect_identical(term_coef_start(term_build(linpar(~x), dd)), c(0, 0))
  expect_identical(term_coef_start(term_build(ridge(~x), dd)), 0)

  # and a vector of zeros really is degenerate for a jump: every
  # break-point collapses onto the same clamped position
  z <- term_refresh(b, numeric(4))
  expect_equal(diff(seg_psi(z)), 0, tolerance = 1e-12)
  expect_lt(qr(as.matrix(term_matrix(z)))$rank, 4L)
})

test_that("a fitted break-point term reports psi, not the pair it is read from", {
  set.seed(32)
  n <- 400
  dd <- data.frame(x = sort(runif(n, 0, 10)))
  dd$y <- 1 + 0.4 * dd$x + 2 * pmax(dd$x - 6, 0) + 1.5 * (dd$x > 6) +
    rnorm(n, sd = 0.3)

  for (kind in c("seg", "jump", "jseg")) {
    spec <- switch(kind, seg = seg(x), jump = jump(x), jseg = jseg(x))
    built <- term_build(seg_start(spec, dd, dd$y), dd)
    cf <- seg_iterate(built, dd$y, iters = 200)
    cur <- term_refresh(built, cf)
    rd <- term_readable(cur, cf)

    # the positions it reports ARE the ones seg_psi() gives
    expect_equal(rd$value[grep("^psi", rd$name)],
                 as.numeric(seg_psi(cur, cf)), tolerance = 1e-10,
                 info = kind)
    # every quantity of the model is named, and the auxiliary pair is not
    expect_true(all(grepl("^(beta|gamma|delta|psi)[0-9]*$", rd$name)),
                info = kind)
    expect_false(any(grepl("^g[0-9]", rd$name)), info = kind)

    # the jacobian against numDeriv, which shares no arithmetic with the
    # delta method written into the method. The joint construction reads
    # its position from a quadratic whose previous iterate is a state of
    # the run and not a function of the coefficients, so what is
    # differentiated there is the fixed point's reading; the changes are
    # checked whatever the construction.
    f <- function(v) {
      q <- term_readable(term_refresh(cur, v), v)
      q$value
    }
    J <- numDeriv::jacobian(f, cf)
    rows <- if (kind == "jseg") !grepl("^psi", rd$name) else
      rep(TRUE, length(rd$name))
    expect_equal(unname(rd$jacobian[rows, , drop = FALSE]),
                 J[rows, , drop = FALSE], tolerance = 1e-7, info = kind)
  }

  # a developed coefficient has no single position to report, and the
  # method says so by answering NULL rather than by inventing one
  d2 <- dd
  d2$g <- factor(rep(c("a", "b"), length.out = n))
  bd <- term_build(seg(x, by = ~0 + g), d2)
  expect_null(term_readable(bd, bd@blueprint$coef))
})


test_that("a continuous break-point term reports how its block moves", {
  # term_block_contract() is what a marginal criterion's gradient needs from a
  # term whose block is not a fixed design. The base method's zeros are right
  # for a fixed one and wrong here, and left in place they cost the outer
  # gradient of a penalized seg 6.6e-04 relative against a central difference
  # of the criterion, where an nl sits at 4.3e-10.
  #
  # The reference is a brute-force dX/dbeta and NOT the criterion: at two
  # break-points the criterion's own central difference disagrees with itself
  # by 7.37 relative, so it cannot judge anything there.
  set.seed(31)
  m <- 6L
  n <- 402L
  id <- factor(rep(seq_len(m), each = n / m))
  x <- stats::runif(n, 0, 10)
  d <- data.frame(x = x, id = id,
                  y = 1 + 0.3 * x + 1.5 * pmax(x - 5, 0) +
                      stats::rnorm(n, sd = 0.4))

  brute <- function(tm, cf, A, cols, dirs, h = 1e-6) {
    vapply(dirs, function(c1) {
      cp <- cf; cp[c1] <- cp[c1] + h
      cm <- cf; cm[c1] <- cm[c1] - h
      Xp <- as.matrix(term_matrix(term_refresh(tm, cp)))[, cols, drop = FALSE]
      Xm <- as.matrix(term_matrix(term_refresh(tm, cm)))[, cols, drop = FALSE]
      sum(A[, cols, drop = FALSE] * (Xp - Xm)) / (2 * h)
    }, numeric(1))
  }

  for (call in list(quote(seg(x)), quote(seg(x, npsi = 2)),
                    quote(seg(x, gamma1 ~ 0 + id)),
                    quote(seg(x, psi ~ 0 + id)))) {
    tm <- term_build(eval(call), d)
    cf <- term_coef_start(tm)
    tm <- term_refresh(tm, cf)
    nm <- term_coef_names(tm)
    set.seed(5)
    A <- matrix(stats::rnorm(n * length(nm)), n, length(nm))
    gcols <- grep("[.]gamma", nm)
    pcols <- grep("[.]psi", nm)

    # the change columns are (x - psi)_+, differentiable in every direction
    A1 <- A; A1[, -gcols] <- 0
    ex <- term_block_contract(tm, cf, A1)
    bf <- brute(tm, cf, A1, gcols, seq_along(nm))
    expect_lt(max(abs(ex - bf)), 1e-6 * max(1, max(abs(bf))))

    # the break-point columns differentiated in the CHANGES, where they are
    # linear. In the break-point their derivative is a delta and the formula
    # takes its almost-everywhere value, which no difference can confirm.
    A2 <- A; A2[, -pcols] <- 0
    ex2 <- term_block_contract(tm, cf, A2)
    bf2 <- brute(tm, cf, A2, pcols, gcols)
    expect_lt(max(abs(ex2[gcols] - bf2)), 1e-6 * max(1, max(abs(bf2))))
  }
})


test_that("a discontinuous break-point term keeps the base method's zeros", {
  # Their position is read off a product of the unknowns and the weight they
  # carry has an unbounded derivative in the break-point, so the block is a
  # working linearization rather than a Jacobian. Zeros are what the base
  # class gives and are the honest answer; this pins that they are what comes
  # back, so a later partial implementation cannot arrive unnoticed.
  set.seed(31)
  n <- 402L
  x <- stats::runif(n, 0, 10)
  d <- data.frame(x = x,
                  y = 1 + 0.3 * x + 2 * (x > 5) + stats::rnorm(n, sd = 0.4))
  for (call in list(quote(jump(x)), quote(jseg(x)))) {
    tm <- term_build(eval(call), d)
    cf <- term_coef_start(tm)
    tm <- term_refresh(tm, cf)
    nm <- term_coef_names(tm)
    A <- matrix(1, n, length(nm))
    expect_identical(term_block_contract(tm, cf, A), numeric(length(nm)))
  }
})


test_that("the block of a discontinuous construction is not a Jacobian", {
  # what the fitting layer routes on: the continuous construction's block is
  # the exact derivative of its contribution, the frozen-weight ones' is a
  # working linearization, and an ordinary term answers TRUE trivially
  expect_true(term_jacobian_block(seg(x)))
  expect_false(term_jacobian_block(jump(x)))
  expect_false(term_jacobian_block(jseg(x)))
  expect_true(term_jacobian_block(linpar(~x)))
})

test_that("crossed break-point lineages are relabeled at a refresh", {
  set.seed(5)
  d <- data.frame(x = sort(stats::runif(300, 0, 10)))
  # a jump whose coefficients imply psi = (7, 3): delta = (1, 1) and
  # g = (-7, -3), i.e. the triples out of order
  tm <- term_build(jump(x, npsi = 2, psi = c(2, 8)), d)
  cf <- term_coef_start(tm)
  cf[tm@blueprint$index[["psi1"]]] <- -7
  cf[tm@blueprint$index[["psi2"]]] <- -3
  tm2 <- term_refresh(tm, cf)
  psi <- seg_psi(tm2)
  expect_equal(psi, sort(psi))
  expect_equal(psi, c(3, 7), tolerance = 1e-8)
  # relabeling moves no value: the contribution is the same sum over k
  expect_equal(tm2@blueprint$value,
               1 * (d$x > 3) + 1 * (d$x > 7), tolerance = 1e-12)
  # and the coefficients the term stored are the relabeled ones, which is
  # what a caller continues from
  cf2 <- tm2@blueprint$coef
  expect_equal(unname(cf2[tm@blueprint$index[["psi1"]]]), -3)
  expect_equal(unname(cf2[tm@blueprint$index[["psi2"]]]), -7)
})

test_that("ordered lineages are left exactly alone", {
  set.seed(6)
  d <- data.frame(x = sort(stats::runif(200, 0, 10)))
  tm <- term_build(jseg(x, npsi = 2, psi = c(3, 7)), d)
  cf <- term_coef_start(tm)
  tm2 <- term_refresh(tm, cf)
  psi <- seg_psi(tm2)
  expect_equal(psi, sort(psi))
  expect_equal(unname(tm2@blueprint$coef), unname(cf))
})

test_that("n_boot is declared on the term and validated", {
  expect_identical(seg(x)@spec$n_boot, 10L)
  expect_identical(jump(x, n_boot = 0)@spec$n_boot, 0L)
  expect_identical(jseg(x, n_boot = 3)@spec$n_boot, 3L)
  expect_error(seg(x, n_boot = -1), "n_boot")
  expect_error(jump(x, n_boot = 1.5), "n_boot")
  expect_error(jseg(x, n_boot = c(1, 2)), "n_boot")
})


test_that("seg_reheat resets the schedule and seg_relocate places psi", {
  set.seed(9)
  d <- data.frame(x = sort(stats::runif(200, 0, 10)))
  b <- term_build(jump(x, psi = 4), d)
  cf <- term_coef_start(b)
  for (i in 1:6) b <- term_refresh(b, cf)
  b2 <- seg_reheat(b)
  expect_equal(b2@blueprint$cscale, rep(0.05, 1))
  expect_identical(b2@blueprint$nref, 0L)
  expect_true(is.na(seg_step(b2)))
  # relocation: positions given come back exactly, the changes kept
  b3 <- seg_relocate(b, 6)
  expect_equal(unname(seg_psi(b3)), 6)
  expect_equal(b3@blueprint$value, as.numeric(d$x > 6) *
                 b@blueprint$coef[b@blueprint$index[["delta1"]]],
               tolerance = 1e-12)
  # confined and sorted
  bj <- term_build(jseg(x, npsi = 2), d)
  bj2 <- seg_relocate(bj, c(9.9, 0.01))
  psi <- unname(seg_psi(bj2))
  expect_equal(psi, sort(psi))
  expect_gte(psi[1], bj@blueprint$lim[1])
  expect_lte(psi[2], bj@blueprint$lim[2])
  expect_error(seg_relocate(bj, 1), "positions")
})


test_that("the profile is read at a point and swept, weighted or not", {
  set.seed(12)
  d <- data.frame(x = sort(stats::runif(400, 0, 10)))
  d$y <- 2 * (d$x > 3) - 1.5 * (d$x > 7) + stats::rnorm(400, sd = 0.3)
  b <- term_build(jump(x, npsi = 2, psi = c(1, 2)), d)
  # the truth's profile beats a wrong placement's
  expect_gt(seg_profile_rss(b, d$y),
            seg_profile_rss(seg_relocate(b, c(3, 7)), d$y))
  # the sweep finds it from the wrong placement
  psi <- unname(seg_psi(seg_polish(b, d$y)))
  expect_equal(psi, c(3, 7), tolerance = 0.15)
  # unit weights are the unweighted profile exactly
  expect_equal(seg_profile_rss(b, d$y, weights = rep(1, 400)),
               seg_profile_rss(b, d$y), tolerance = 1e-12)
  # and a weighted sweep runs and returns positions in the interval
  w <- tabulate(sample.int(400, 400, replace = TRUE), nbins = 400)
  pw <- unname(seg_psi(seg_polish(b, d$y, weights = w)))
  expect_true(all(pw >= b@blueprint$lim[1] & pw <= b@blueprint$lim[2]))
  expect_error(seg_profile_rss(b, d$y, weights = rep(-1, 400)), "weights")
})
