# The second derivative of a design block in its coefficients.
#
# It enters the HESSIAN of a marginal criterion and nothing else, so none of
# these checks can be reached through a fit: they are put to the generic
# directly, against a reference that shares no arithmetic with the closed form.

set.seed(4)
nb2 <- 60
db2 <- data.frame(x = seq(0.2, 3, length.out = nb2),
                  z = runif(nb2, -1, 1), w = runif(nb2, 0.5, 2),
                  g = factor(rep(c("a", "b", "c"), length.out = nb2)),
                  id = factor(rep(seq_len(6), length.out = nb2)))
db2$y <- 2 * exp(-1.3 * db2$x) + stats::rnorm(nb2, sd = 0.1)

# ONE central difference of term_block_deriv along u. A single stencil on an
# analytic quantity is the reference the toolkit sanctions; differencing
# term_matrix twice would be the nesting it forbids everywhere.
bd2_ref <- function(term, cf, v, u, h) {
  (term_block_deriv(term, coef = cf + h * u, v = v) -
   term_block_deriv(term, coef = cf - h * u, v = v)) / (2 * h)
}

bd2_cases <- function() {
  list(
    "scalar parameters" =
      term_build(nl(~ a * exp(-r * x), start = list(a = 2, r = 1.3)), db2),
    "log link on r" =
      term_build(nl(~ a * exp(-r * x), start = list(a = 2, r = 1.3),
                    links = list(r = linkfunctions7::log_link())), db2),
    "logit on a, log on r" =
      term_build(nl(~ a * exp(-r * x), start = list(a = 0.4, r = 1.3),
                    links = list(a = linkfunctions7::logit_link(),
                                 r = linkfunctions7::log_link())), db2),
    "a subformula" =
      term_build(nl(~ a * exp(-r * x), a ~ 0 + g, start = list(r = 1.3)), db2),
    "a SPARSE subformula" =
      term_build(nl(~ a * exp(-r * x), a ~ 0 + random(~ 1 | id),
                    start = list(r = 1.3)), db2),
    "three parameters" =
      term_build(nl(~ a * exp(-r * x) + s * z,
                    start = list(a = 2, r = 1.3, s = 0.5)), db2))
}

test_that("the second derivative of the block is symmetric in its directions", {
  # Mixed partial derivatives are equal, so the two directions may be
  # exchanged. The equality is to a rounding of the accumulated value and NOT
  # to the bit: the (p2, p3) grid is summed in a fixed order, so exchanging
  # the directions re-orders that sum. Measured, the gap is 1e-16 relative;
  # an implementation pairing a direction with the wrong parameter is out by
  # order one and fails this by fifteen orders.
  set.seed(7)
  for (nm in names(bd2_cases())) {
    tm <- bd2_cases()[[nm]]
    k <- term_npar(tm)
    v <- rnorm(k)
    u <- rnorm(k)
    A <- term_block_deriv2(tm, v = v, u = u)
    B <- term_block_deriv2(tm, v = u, u = v)
    expect_lt(max(abs(A - B)) / max(abs(A)), 1e-13)
  }
})

test_that("a fixed design gets exactly zero, and a curved nl does not", {
  # Zero is exact here rather than approximate: a design that does not move
  # with its coefficients has a second derivative that is identically zero.
  fixed <- list(
    linpar = term_build(linpar(~ x + z), db2),
    smooth = term_build(s(x, k = 8), db2),
    random = term_build(random(~ 1 | id), db2),
    ridge = term_build(ridge(~ x + z), db2),
    tensor = term_build(te(x, z, k = 4), db2))
  for (nm in names(fixed)) {
    tm <- fixed[[nm]]
    k <- term_npar(tm)
    D <- term_block_deriv2(tm, v = rep(1, k), u = rep(1, k))
    expect_identical(dim(D), dim(as.matrix(term_matrix(tm))))
    expect_true(all(D == 0))
  }
  # the negative control: without it the checks above are satisfied by a
  # method returning zeros everywhere
  tm <- bd2_cases()[["scalar parameters"]]
  expect_gt(max(abs(term_block_deriv2(tm, v = c(1, 1), u = c(1, 1)))), 1e-3)
})

test_that("the second derivative matches one stencil on the first", {
  # A gap that falls as h^2 is the reference's truncation and not an error in
  # the closed form; the ratio between two steps is what says so. Measured, it
  # is 100 to a tenth on every analytic case.
  set.seed(11)
  for (nm in names(bd2_cases())) {
    tm <- bd2_cases()[[nm]]
    cf <- tm@blueprint$coef
    k <- length(cf)
    v <- rnorm(k)
    u <- rnorm(k)
    A <- term_block_deriv2(tm, v = v, u = u)
    gaps <- vapply(c(1e-3, 1e-4), function(h) {
      R <- bd2_ref(tm, cf, v, u, h)
      max(abs(A - R)) / max(1e-12, max(abs(R)))
    }, numeric(1))
    expect_lt(gaps[2L], 1e-6)
    expect_gt(gaps[1L] / gaps[2L], 10)
  }
})

test_that("an opaque function is right to the accuracy its stencils have", {
  # nl() reads an opaque f by differencing, so BOTH sides of the check above
  # carry that noise and the reference amplifies it by 1/h -- measured, the
  # gap GROWS by ten per decade there. What judges the closed form instead is
  # the same f written symbolically, which is exact and shares no arithmetic
  # with the stencils. The ladder is one stencil per order: the block agrees
  # to 3.8e-13, the first derivative to 5.4e-10, the second to 3.3e-07.
  sym <- term_build(nl(~ a * exp(-r * x), start = list(a = 2, r = 1.3)), db2)
  opq <- term_build(nl(function(x, theta) theta$a * exp(-theta$r * x),
                       params = c("a", "r"), x = x,
                       start = list(a = 2, r = 1.3)), db2)
  expect_identical(sym@blueprint$mode, "symbolic")
  expect_identical(opq@blueprint$mode, "numeric")
  set.seed(3)
  v <- rnorm(2)
  u <- rnorm(2)
  A <- term_block_deriv2(sym, v = v, u = u)
  B <- term_block_deriv2(opq, v = v, u = u)
  expect_lt(max(abs(A - B)) / max(abs(A)), 1e-5)
})

test_that("a bilinear f under identity links gives exactly zero", {
  # f = a x + b z + a b w is jointly QUADRATIC in (a, b), so its third
  # derivative in the parameters vanishes and every addend of the closed form
  # that survives carries an h''. Under the identity link h'' and h''' are
  # zero, so the answer is zero exactly and not merely small. It is the
  # control that identifies what this method adds over term_block_deriv.
  bil <- list(
    plain = term_build(nl(~ a * x + b * z + a * b * w,
                          start = list(a = 1.4, b = 0.8)), db2),
    developed = term_build(nl(~ a * x + b * z + a * b * w, a ~ 0 + g,
                              start = list(b = 0.8)), db2),
    sparse = term_build(nl(~ a * x + b * z + a * b * w,
                           a ~ 0 + random(~ 1 | id),
                           start = list(b = 0.8)), db2))
  set.seed(13)
  for (nm in names(bil)) {
    tm <- bil[[nm]]
    k <- term_npar(tm)
    expect_true(all(term_block_deriv2(tm, v = rnorm(k), u = rnorm(k)) == 0))
  }
})

test_that("the qualification on the links is necessary", {
  # With a link other than the identity the addends carrying h'' and h'''
  # survive even though f itself is quadratic, so the bilinear cell stops
  # being a zero control. Without this test the one above could be satisfied
  # by a method that dropped those addends.
  tm <- term_build(nl(~ a * x + b * z + a * b * w,
                      start = list(a = 1.4, b = 0.8),
                      links = list(b = linkfunctions7::log_link())), db2)
  set.seed(13)
  k <- term_npar(tm)
  v <- rnorm(k)
  u <- rnorm(k)
  D <- term_block_deriv2(tm, v = v, u = u)
  expect_gt(max(abs(D)), 1e-6)
  R <- bd2_ref(tm, tm@blueprint$coef, v, u, 1e-4)
  expect_lt(max(abs(D - R)) / max(abs(R)), 1e-6)
})

test_that("a sparse development is scaled in its own storage", {
  # the return is dense by contract, but the intermediate must not densify
  # the sub-design: at order three this runs once per PAIR of hyperparameters
  tm <- term_build(nl(~ a * exp(-r * x), a ~ 0 + random(~ 1 | id),
                      start = list(r = 1.3)), db2)
  Z <- tm@blueprint$Z[["a"]]
  expect_true(inherits(Z, "Matrix"))
  acc <- runif(nrow(db2))
  expect_true(inherits(Z * acc, "Matrix"))
  expect_equal(as.matrix(Z * acc), as.matrix(Z) * acc)
})

test_that("the directions are checked against the term's own width", {
  tm <- term_build(nl(~ a * exp(-r * x), start = list(a = 2, r = 1.3)), db2)
  expect_error(term_block_deriv2(tm, v = 1, u = c(1, 1)), "'v' must have")
  expect_error(term_block_deriv2(tm, v = c(1, 1), u = 1), "'u' must have")
})


# ---------------------------------------------------------------------------
# The smoothed break-point terms, the second implementer.  With smoothed = an
# abs_smoother the block IS the true Jacobian, so the same closed forms apply
# one order further up; the sharp constructions answer zeros.

set.seed(31)
nsm <- 240
dsm <- data.frame(x = sort(runif(nsm, 0, 10)),
                  id = factor(rep(seq_len(6), length.out = nsm)),
                  g = factor(rep(c("a", "b"), length.out = nsm)))
dsm$y <- 1 + 0.3 * dsm$x + 1.5 * pmax(dsm$x - 5, 0) +
  2 * (dsm$x > 5) + stats::rnorm(nsm, sd = 0.4)

sm_cases <- function() {
  sp <- penalties7::smooth_probit()
  sh <- penalties7::smooth_hyperbolic()
  list(
    "seg probit" = term_build(seg(x, smoothed = sp), dsm),
    "jump probit" = term_build(jump(x, smoothed = sp), dsm),
    "jseg probit" = term_build(jseg(x, smoothed = sp), dsm),
    "seg hyperbolic" = term_build(seg(x, smoothed = sh), dsm),
    "jseg two break-points" = term_build(jseg(x, npsi = 2, smoothed = sp),
                                         dsm),
    "jseg, psi ~ 0 + id" = term_build(jseg(x, psi ~ 0 + id, smoothed = sp),
                                      dsm),
    "jseg, psi random" = term_build(jseg(x, psi ~ random(~ 1 | id),
                                         smoothed = sp), dsm),
    "seg, gamma1 ~ 0 + g" = term_build(seg(x, gamma1 ~ 0 + g, smoothed = sp),
                                       dsm))
}

test_that("a smoothed break-point term is symmetric in its directions", {
  set.seed(5)
  for (nm in names(sm_cases())) {
    tm <- sm_cases()[[nm]]
    k <- length(term_coef_names(tm))
    v <- rnorm(k)
    u <- rnorm(k)
    A <- term_block_deriv2(tm, v = v, u = u)
    B <- term_block_deriv2(tm, v = u, u = v)
    expect_lt(max(abs(A - B)) / max(abs(A)), 1e-13)
  }
})

test_that("a sharp break-point term gets zeros and a smoothed one does not", {
  # The two zeros have different reasons. For the continuous construction the
  # second derivative really is zero away from the break-points: the truncated
  # line's derivative in the position is an indicator, whose own derivative is
  # a point mass, and the position column is linear in the change. For a jump
  # or a jseg the block is a working linearization with a frozen weight rather
  # than a Jacobian, so the question is not the one its columns answer.
  for (tm in list(term_build(seg(x), dsm), term_build(jump(x), dsm),
                  term_build(jseg(x), dsm))) {
    k <- length(term_coef_names(tm))
    expect_true(all(term_block_deriv2(tm, v = rep(1, k),
                                      u = rep(1, k)) == 0))
  }
  for (nm in names(sm_cases())) {
    tm <- sm_cases()[[nm]]
    k <- length(term_coef_names(tm))
    expect_gt(max(abs(term_block_deriv2(tm, v = rep(1, k), u = rep(1, k)))),
              1e-6)
  }
})

test_that("the smoothed second derivative matches one stencil on the first", {
  set.seed(17)
  for (nm in names(sm_cases())) {
    tm <- sm_cases()[[nm]]
    cf <- tm@blueprint$coef
    k <- length(cf)
    v <- rnorm(k)
    u <- rnorm(k)
    A <- term_block_deriv2(tm, v = v, u = u)
    gaps <- vapply(c(1e-3, 1e-4), function(h) {
      R <- bd2_ref(tm, cf, v, u, h)
      max(abs(A - R)) / max(1e-12, max(abs(R)))
    }, numeric(1))
    expect_lt(gaps[2L], 1e-4)
    expect_gt(gaps[1L] / gaps[2L], 10)
  }
})

test_that("a confined break-point contributes exactly zero, unlike order one", {
  # Every addend of the second derivative carries a direction in the
  # break-point, so where the position sits against a confinement limit the
  # whole contribution vanishes exactly. The FIRST derivative does not: the
  # position column is -gamma S(u), whose derivative in the CHANGE carries no
  # break-point direction at all. That asymmetry is what makes this a control
  # on the second derivative rather than on the gate the two share.
  sp <- penalties7::smooth_probit()
  for (tm in list(term_build(seg(x, smoothed = sp), dsm),
                  term_build(jump(x, smoothed = sp), dsm),
                  term_build(jseg(x, smoothed = sp), dsm))) {
    cf <- tm@blueprint$coef
    nms <- term_coef_names(tm)
    cf[grep("psi[0-9]*$", nms)] <- tm@blueprint$lim[2L] + 5
    k <- length(cf)
    set.seed(2)
    expect_true(all(term_block_deriv2(tm, coef = cf, v = rnorm(k),
                                      u = rnorm(k)) == 0))
    expect_gt(max(abs(term_block_deriv(tm, coef = cf, v = rnorm(k)))), 1e-6)
  }
})

test_that("under the quintic the answer is exactly zero outside the width", {
  # smooth_quintic() is exact outside [-h, h], so every derivative of the
  # smoother above the first vanishes there and so does the whole second
  # derivative of the block -- exactly, not nearly. It is the control that
  # says the method reads the smoother's own higher orders and not something
  # standing in for them.
  tm <- term_build(seg(x, smoothed = penalties7::smooth_quintic()), dsm)
  bp <- tm@blueprint
  k <- length(bp$coef)
  set.seed(3)
  D <- term_block_deriv2(tm, v = rnorm(k), u = rnorm(k))
  uu <- bp$xv - bp$coef[bp$index[["psi1"]]]
  far <- abs(uu) > bp$smooth$width
  expect_gt(sum(far), 0)
  expect_gt(sum(!far), 0)
  expect_true(all(D[far, ] == 0))
  expect_gt(max(abs(D[!far, , drop = FALSE])), 1e-6)
})

test_that("the smoothed directions are checked against the term's width", {
  tm <- term_build(seg(x, smoothed = penalties7::smooth_probit()), dsm)
  k <- length(term_coef_names(tm))
  expect_error(term_block_deriv2(tm, v = 1, u = rep(1, k)), "'v' must have")
  expect_error(term_block_deriv2(tm, v = rep(1, k), u = 1), "'u' must have")
})


# ---------------------------------------------------------------------------
# The sharp break-point terms: what licenses the zeros.
#
# For the CONTINUOUS construction the claim is that the second derivative is
# zero ALMOST everywhere, and "almost" is the whole content: away from the
# break-points it is zero exactly, at every step, and at an observation
# sitting on a break-point it is a point mass, which is not a number. For the
# discontinuous ones the claim is different -- their block is not a Jacobian
# at all -- and the two are checked separately.

set.seed(31)
nsh <- 400
dsh <- data.frame(x = sort(runif(nsh, 0, 10)))
dsh$y <- 1 + 0.3 * dsh$x + 1.5 * pmax(dsh$x - 5, 0) +
  2 * (dsh$x > 5) + stats::rnorm(nsh, sd = 0.4)

bd2_steps <- c(1e-3, 2.5e-4, 6.25e-5)     # h, h/4, h/16

test_that("a sharp seg's second derivative is exactly zero away from psi", {
  # Not small: zero. The truncated line's derivative in the position is an
  # indicator, whose own derivative is a point mass, and the position column
  # is linear in the change -- so on every observation the perturbation does
  # not carry across the break-point the difference is identically zero,
  # whatever the step.
  tm <- term_build(seg(x), dsh)
  bp <- tm@blueprint
  ip <- bp$index[["psi1"]]
  k <- length(bp$coef)
  set.seed(7)
  v <- rnorm(k)
  u <- rnorm(k)
  j <- which.min(abs(bp$xv - bp$coef[ip]))
  for (place in c("between", "on")) {
    cf <- bp$coef
    cf[ip] <- if (place == "on") bp$xv[j] else (bp$xv[j] + bp$xv[j + 1L]) / 2
    for (h in bd2_steps) {
      D <- (term_block_deriv(tm, coef = cf + h * u, v = v) -
            term_block_deriv(tm, coef = cf - h * u, v = v)) / (2 * h)
      far <- abs(bp$xv - cf[ip]) > abs(h * u[ip])
      expect_true(all(D[far, ] == 0))
    }
  }
})

test_that("the one exception is a point mass and behaves like one", {
  # The positive control. With an observation exactly on the break-point the
  # difference on THAT row grows as 1/h -- measured 598, 2394, 9574 at h, h/4,
  # h/16, four times per step -- which is what a point mass looks like on a
  # difference quotient and is exactly what "almost everywhere" sets aside.
  # Without this the test above is satisfied by a sample that happens to have
  # no observation near the break-point, which asks nothing.
  tm <- term_build(seg(x), dsh)
  bp <- tm@blueprint
  ip <- bp$index[["psi1"]]
  k <- length(bp$coef)
  set.seed(7)
  v <- rnorm(k)
  u <- rnorm(k)
  cf <- bp$coef
  cf[ip] <- bp$xv[which.min(abs(bp$xv - bp$coef[ip]))]
  m <- vapply(bd2_steps, function(h) {
    D <- (term_block_deriv(tm, coef = cf + h * u, v = v) -
          term_block_deriv(tm, coef = cf - h * u, v = v)) / (2 * h)
    max(abs(D))
  }, numeric(1))
  expect_gt(m[2L] / m[1L], 3)
  expect_gt(m[3L] / m[2L], 3)
})

test_that("order one IS a derivative for seg and is not for jump or jseg", {
  # This is why the zeros differ in kind. For the continuous construction the
  # block's own difference converges onto term_block_deriv(), so its first
  # derivative exists and the SECOND is what vanishes. For a jump or a jseg
  # term_block_deriv() answers zeros while the block moves by five orders of
  # magnitude, because the block is a working linearization with a frozen
  # weight rather than a Jacobian: the zeros there are a refusal and not a
  # description.
  blk <- function(tm, cf, cs) {
    b <- tm@blueprint
    as.matrix(.seg_assemble(b, b$xv, cf, cscale = cs)$X)
  }
  tm <- term_build(seg(x), dsh)
  bp <- tm@blueprint
  set.seed(7)
  u <- rnorm(length(bp$coef))
  A <- term_block_deriv(tm, coef = bp$coef, v = u)
  for (h in bd2_steps) {
    D <- (blk(tm, bp$coef + h * u, bp$cscale) -
          blk(tm, bp$coef - h * u, bp$cscale)) / (2 * h)
    expect_lt(mean(abs(D - A)), 1e-10)
  }
  for (kind in c("jump", "jseg")) {
    tj <- term_build(if (kind == "jump") jump(x) else jseg(x), dsh)
    bpj <- tj@blueprint
    kj <- length(bpj$coef)
    set.seed(7)
    uj <- rnorm(kj)
    expect_true(all(term_block_deriv(tj, v = uj) == 0))
    # at the factor the annealing schedule descends to, the weight is steep
    D <- (blk(tj, bpj$coef + 1e-3 * uj, 5e-4) -
          blk(tj, bpj$coef - 1e-3 * uj, 5e-4)) / 2e-3
    expect_gt(max(abs(D)), 1e3)
  }
})
