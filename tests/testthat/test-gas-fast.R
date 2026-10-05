# The fast route of the score-driven filter: the C registries of
# distributions7 and linkfunctions7 against the R callbacks they replace.
# The scalar entry points mirror the vector kernels expression by expression
# and the composition mirrors distrib_kernel()'s -- but a compiler is free
# to CONTRACT a multiply-add into an FMA where R's interpreter never does,
# and clang on the arm64 macOS runner did exactly that (last-bit differences
# in the jacobian, 2026-08-19), so the cross-route twin carries a tolerance
# of 1e-13: tight enough that a wrong composition, off by the factor it
# drops, still fails it. What stays identical() is everything computed by
# ONE route both times -- the inert-context fallbacks, and the threads over
# groups, which may not move a bit.

fast_panel <- function(seed, groups, per) {
  set.seed(seed)
  id <- factor(rep(seq_len(groups), each = per))
  n <- length(id)
  d <- data.frame(id = id, t = rep(seq_len(per), groups),
                  y = abs(rnorm(n, 3, 1)) + 0.1)
  d
}

cb_from_kernel <- function(d7, param, theta, y) {
  k <- distributions7::distrib_kernel(d7, param)
  th_at <- function(i) lapply(theta, function(v) v[[min(i, length(v))]])
  list(score = function(e, i) as.numeric(k$score(y[i], th_at(i), e)),
       curvature = function(e, i) as.numeric(k$curvature(y[i], th_at(i), e)))
}

# the registry's name and the constants after the parameters, as
# statmodels7's structural_callbacks() builds them
fast_spec <- function(d7, param, theta, y) {
  params <- d7@params
  lk <- d7@link_params[[param]]
  rt <- distributions7::distrib_scalar_route(d7)
  lr <- linkfunctions7::link_scalar_route(lk)
  list(family = rt$name, link = lr$name, link_par = lr$par,
       k = match(param, params), bounds = as.numeric(lk@link_bounds),
       y = as.numeric(y), theta = c(theta, rt$constants))
}

run_both <- function(term, dd, psi, d7, param, theta, threads = 1L) {
  tb <- term_build(term, dd)
  eta0 <- rep(0.1, nrow(dd))
  cb <- cb_from_kernel(d7, param, theta, dd$y)
  ref <- term_filter(tb, eta0, dd$y, cb$score, cb$curvature, psi)
  fst <- term_filter(tb, eta0, dd$y, cb$score, cb$curvature, psi,
                     fast = fast_spec(d7, param, theta, dd$y),
                     threads = threads)
  list(ref = ref, fst = fst)
}

test_that("the fast route reproduces the callbacks, scalar filter", {
  dd <- fast_panel(1, groups = 10, per = 25)
  psi <- list(omega = 0.2, alpha1 = 0.15, pacf1 = 0.4)
  for (case in list(
    list(d7 = distributions7::gaussian1_distrib(), param = "mu",
         theta = list(mu = 0, sigma = 1.3)),
    list(d7 = distributions7::gaussian1_distrib(), param = "sigma",
         theta = list(mu = 0.5, sigma = 1)),
    list(d7 = distributions7::gamma1_distrib(), param = "mu",
         theta = list(mu = 3, phi = 0.4)))) {
    both <- run_both(gas(p = 1, q = 1, by = id, time = t), dd, psi,
                     case$d7, case$param, case$theta)
    expect_equal(both$ref, both$fst, tolerance = 1e-13)
  }
})

test_that("the fast route reads a family's constants after its parameters", {
  set.seed(4)
  dd <- fast_panel(4, groups = 10, per = 25)
  dd$y <- stats::rbinom(nrow(dd), size = 10, prob = 0.4)
  psi <- list(omega = 0.1, alpha1 = 0.05, pacf1 = 0.4)
  d7 <- distributions7::binomial_distrib(size = 10)
  both <- run_both(gas(p = 1, q = 1, by = id, time = t), dd, psi, d7, "mu",
                   list(mu = 0.4))
  expect_equal(both$ref, both$fst, tolerance = 1e-13)
  # the size reaches the entry: read as one trial, the filter moves
  cb <- cb_from_kernel(d7, "mu", list(mu = 0.4), dd$y)
  wrong <- term_filter(term_build(gas(p = 1, q = 1, by = id, time = t), dd),
                       rep(0.1, nrow(dd)), dd$y, cb$score, cb$curvature, psi,
                       fast = utils::modifyList(
                         fast_spec(d7, "mu", list(mu = 0.4), dd$y),
                         list(theta = list(mu = 0.4, size = 1))))
  expect_gt(max(abs(wrong$eta - both$fst$eta)), 1e-3)
})

test_that("the fast route reproduces the callbacks on the submodel route too", {
  dd <- fast_panel(2, groups = 10, per = 25)
  tb <- term_build(gas(p = 1, q = 1, omega ~ random(~1 | id),
                       by = id, time = t), dd)
  nm <- term_params(tb)
  psi <- stats::setNames(as.list(c(0.2, rep(0.05, length(nm) - 3L),
                                   0.15, 0.4)),
                         nm)
  d7 <- distributions7::gaussian1_distrib()
  theta <- list(mu = 0, sigma = 1.2)
  cb <- cb_from_kernel(d7, "mu", theta, dd$y)
  eta0 <- rep(0.1, nrow(dd))
  ref <- term_filter(tb, eta0, dd$y, cb$score, cb$curvature, psi)
  fst <- term_filter(tb, eta0, dd$y, cb$score, cb$curvature, psi,
                     fast = fast_spec(d7, "mu", theta, dd$y))
  expect_equal(ref, fst, tolerance = 1e-13)
})

test_that("the result does not depend on threads, bit for bit", {
  dd <- fast_panel(3, groups = 12, per = 20)   # above the group threshold
  psi <- list(omega = 0.2, alpha1 = 0.15, pacf1 = 0.4)
  d7 <- distributions7::gaussian1_distrib()
  theta <- list(mu = 0, sigma = 1.1)
  b1 <- run_both(gas(p = 1, q = 1, by = id, time = t), dd, psi, d7, "mu",
                 theta, threads = 1L)
  b2 <- run_both(gas(p = 1, q = 1, by = id, time = t), dd, psi, d7, "mu",
                 theta, threads = 2L)
  expect_identical(b1$fst, b2$fst)
  expect_equal(b1$ref, b2$fst, tolerance = 1e-13)
})

test_that("the adjoint rides the fast forward pass and calls nothing back", {
  dd <- fast_panel(5, groups = 10, per = 25)
  d7 <- distributions7::gaussian1_distrib()
  theta <- list(mu = 0, sigma = 1.2)
  cb <- cb_from_kernel(d7, "mu", theta, dd$y)
  eta0 <- rep(0.1, nrow(dd))
  set.seed(6)
  gw <- rnorm(nrow(dd))
  calls <- new.env()
  calls$k <- 0L
  sc_cnt <- function(e, i) { calls$k <- calls$k + 1L; cb$score(e, i) }
  cu_cnt <- function(e, i) { calls$k <- calls$k + 1L; cb$curvature(e, i) }

  # the scalar route
  tb <- term_build(gas(p = 1, q = 1, by = id, time = t), dd)
  psi <- list(omega = 0.2, alpha1 = 0.15, pacf1 = 0.4)
  ref <- term_adjoint(tb, eta0, dd$y, sc_cnt, cu_cnt, psi, g = gw)
  expect_gt(calls$k, 0L)
  calls$k <- 0L
  fst <- term_adjoint(tb, eta0, dd$y, sc_cnt, cu_cnt, psi, g = gw,
                      fast = fast_spec(d7, "mu", theta, dd$y), threads = 2L)
  expect_identical(calls$k, 0L)
  expect_equal(ref, fst, tolerance = 1e-13)

  # the submodel route
  tbs <- term_build(gas(p = 1, q = 1, omega ~ random(~1 | id),
                        by = id, time = t), dd)
  nm <- term_params(tbs)
  psis <- stats::setNames(as.list(c(0.2, rep(0.05, length(nm) - 3L),
                                    0.15, 0.4)), nm)
  calls$k <- 0L
  refs <- term_adjoint(tbs, eta0, dd$y, sc_cnt, cu_cnt, psis, g = gw)
  expect_gt(calls$k, 0L)
  calls$k <- 0L
  fsts <- term_adjoint(tbs, eta0, dd$y, sc_cnt, cu_cnt, psis, g = gw,
                       fast = fast_spec(d7, "mu", theta, dd$y), threads = 2L)
  expect_identical(calls$k, 0L)
  expect_equal(refs, fsts, tolerance = 1e-13)
})

test_that("an uncovered family or link leaves the context inert", {
  dd <- fast_panel(4, groups = 6, per = 20)
  psi <- list(omega = 0.2, alpha1 = 0.15, pacf1 = 0.4)
  d7 <- distributions7::gaussian1_distrib()
  theta <- list(mu = 0, sigma = 1.3)
  tb <- term_build(gas(p = 1, q = 1, by = id, time = t), dd)
  eta0 <- rep(0.1, nrow(dd))
  cb <- cb_from_kernel(d7, "mu", theta, dd$y)
  ref <- term_filter(tb, eta0, dd$y, cb$score, cb$curvature, psi)
  sp <- fast_spec(d7, "mu", theta, dd$y)
  sp$family <- "NoSuchDistrib"
  expect_identical(term_filter(tb, eta0, dd$y, cb$score, cb$curvature, psi,
                               fast = sp), ref)
  sp <- fast_spec(d7, "mu", theta, dd$y)
  sp$link <- "no-such-link"
  expect_identical(term_filter(tb, eta0, dd$y, cb$score, cb$curvature, psi,
                               fast = sp), ref)
})

# gas(scaling = d): the callbacks compose u = s I^-d and its derivative from
# distrib_kernel()'s information and dinformation exactly as statmodels7's
# structural_callbacks() does, and the fast route composes the same quantities
# from the family's d7_info_dinfo
cb_scaled <- function(d7, param, theta, y, d) {
  k <- distributions7::distrib_kernel(d7, param)
  th_at <- function(i) lapply(theta, function(v) v[[min(i, length(v))]])
  list(score = function(e, i) {
         th <- th_at(i)
         as.numeric(k$score(y[i], th, e)) *
           as.numeric(k$information(y[i], th, e))^(-d)
       },
       curvature = function(e, i) {
         th <- th_at(i)
         s <- as.numeric(k$score(y[i], th, e))
         info <- as.numeric(k$information(y[i], th, e))
         as.numeric(k$curvature(y[i], th, e)) * info^(-d) -
           d * s * info^(-d - 1) * as.numeric(k$dinformation(y[i], th, e))
       })
}

test_that("the fast route reproduces the scaled callbacks", {
  dd <- fast_panel(3, groups = 10, per = 25)
  n <- nrow(dd)
  psi <- list(omega = 0.2, alpha1 = 0.15, pacf1 = 0.4)
  set.seed(4)
  for (case in list(
    # the information of a gaussian mean is 1/sigma^2, varying by observation
    list(d7 = distributions7::gaussian1_distrib(), param = "mu", d = 0.5,
         theta = list(mu = 0, sigma = stats::runif(n, 0.8, 1.6))),
    list(d7 = distributions7::gaussian1_distrib(), param = "sigma", d = 1,
         theta = list(mu = 0.5, sigma = 1)),
    list(d7 = distributions7::gamma1_distrib(), param = "mu", d = 1,
         theta = list(mu = 3, phi = 0.4)),
    list(d7 = distributions7::gamma1_distrib(), param = "phi", d = -0.5,
         theta = list(mu = stats::runif(n, 2, 4), phi = 0.4)))) {
    tb <- term_build(gas(p = 1, q = 1, by = id, time = t, scaling = case$d), dd)
    eta0 <- rep(0.1, n)
    if (case$param == "phi") eta0 <- rep(log(0.4), n)
    cb <- cb_scaled(case$d7, case$param, case$theta, dd$y, case$d)
    ref <- term_filter(tb, eta0, dd$y, cb$score, cb$curvature, psi)
    fs <- fast_spec(case$d7, case$param, case$theta, dd$y)
    fs$scaling <- case$d
    fst <- term_filter(tb, eta0, dd$y, cb$score, cb$curvature, psi, fast = fs)
    expect_equal(ref, fst, tolerance = 1e-13,
                 label = paste(case$param, case$d))
    # and the threads over groups move no bit
    fst4 <- term_filter(tb, eta0, dd$y, cb$score, cb$curvature, psi, fast = fs,
                        threads = 4L)
    expect_identical(fst, fst4)
  }
})

test_that("the fast route reproduces the scaled callbacks on the submodel route", {
  dd <- fast_panel(5, groups = 10, per = 25)
  tb <- term_build(gas(p = 1, q = 1, omega ~ random(~1 | id),
                       by = id, time = t, scaling = 1), dd)
  nm <- term_params(tb)
  psi <- stats::setNames(as.list(c(0.2, rep(0.05, length(nm) - 3L),
                                   0.15, 0.4)), nm)
  d7 <- distributions7::gaussian1_distrib()
  set.seed(6)
  theta <- list(mu = 0, sigma = stats::runif(nrow(dd), 0.8, 1.6))
  cb <- cb_scaled(d7, "mu", theta, dd$y, 1)
  eta0 <- rep(0.1, nrow(dd))
  ref <- term_filter(tb, eta0, dd$y, cb$score, cb$curvature, psi)
  fs <- fast_spec(d7, "mu", theta, dd$y)
  fs$scaling <- 1
  fst <- term_filter(tb, eta0, dd$y, cb$score, cb$curvature, psi, fast = fs)
  expect_equal(ref, fst, tolerance = 1e-13)
})
