# The fast route of the score-driven filter: the C registries of
# distributions7 and linkfunctions7 against the R callbacks they replace.
# The two are held IDENTICAL, not close: the scalar entry points mirror the
# vector kernels expression by expression and the composition mirrors
# distrib_kernel()'s, so a last-bit difference is a defect. An uncovered
# family leaves the context inert, and the threads over groups may not move
# a bit either.

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

fast_spec <- function(d7, param, theta, y) {
  params <- d7@params
  lk <- d7@link_params[[param]]
  list(family = attr(S7::S7_class(d7), "name"), link = lk@link_name,
       k = match(param, params), bounds = as.numeric(lk@link_bounds),
       y = as.numeric(y), theta = theta)
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

test_that("the fast route is bit-identical to the callbacks, scalar filter", {
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
    expect_identical(both$ref, both$fst)
  }
})

test_that("the fast route is bit-identical on the submodel route too", {
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
  expect_identical(ref, fst)
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
  expect_identical(b1$ref, b2$fst)
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
  expect_identical(ref, fst)

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
  expect_identical(refs, fsts)
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
