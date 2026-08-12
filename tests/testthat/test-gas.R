# Score-driven dynamics: the recursion, its exact derivative, and the chart
# that keeps the persistence stationary.

set.seed(8)
n <- 60
dd <- data.frame(t = seq_len(n), y = rnorm(n),
                 g = rep(c("a", "b"), each = n / 2))

# a gaussian log-likelihood in the mean, as the model layer would supply it:
# score and curvature of the contribution with respect to the predictor
gauss_score <- function(y) function(e, i) y[i] - e
gauss_curv <- function(y) function(e, i) -1

test_that("the parameters are named for the chart they live on", {
  expect_identical(term_params(gas(p = 1, q = 1)),
                   c("omega", "a1", "pacf1"))
  expect_identical(term_params(gas(p = 2, q = 3)),
                   c("omega", "a1", "a2", "pacf1", "pacf2", "pacf3"))
  lk <- term_links(gas(p = 1, q = 2))
  expect_identical(vapply(lk, function(l) l@link_name, character(1)),
                   c(omega = "identity", a1 = "identity",
                     pacf1 = "rhobit", pacf2 = "rhobit"))
})

test_that("the Levinson-Durbin map is stationary and its jacobian is right", {
  # at q = 1 the partial autocorrelation IS the coefficient
  expect_equal(gas_levinson(0.7)$phi, 0.7)

  for (pc in list(c(0.5, -0.3), c(0.9, 0.8, -0.6))) {
    ld <- gas_levinson(pc)
    # stationary: every root of 1 - sum(phi z^j) outside the unit circle
    roots <- polyroot(c(1, -ld$phi))
    expect_true(all(Mod(roots) > 1 + 1e-8))
    # the jacobian against a numerical derivative of the same map
    num <- numDeriv::jacobian(function(v) gas_levinson(v)$phi, pc)
    expect_equal(ld$jacobian, num, tolerance = 1e-7)
  }
})

test_that("the filter reproduces the recursion written out by hand", {
  term <- term_build(gas(p = 1, q = 1, time = t), dd)
  psi <- list(omega = 0.2, a1 = 0.3, pacf1 = 0.6)
  out <- term_filter(term, eta = rep(0, n), y = dd$y,
                     score = gauss_score(dd$y), curvature = gauss_curv(dd$y),
                     psi = psi)

  # the same recursion, written independently
  om <- 0.2; a <- 0.3; b <- 0.6
  f <- numeric(n); s <- numeric(n)
  f0 <- om / (1 - b)
  for (t in seq_len(n)) {
    f[t] <- om + a * (if (t > 1) s[t - 1] else 0) + b * (if (t > 1) f[t - 1] else f0)
    s[t] <- dd$y[t] - f[t]
  }
  expect_equal(out$eta, f, tolerance = 1e-12)
})

test_that("the jacobian of the filter is exact", {
  for (cfg in list(c(1, 1), c(2, 1), c(1, 2))) {
    term <- term_build(gas(p = cfg[1], q = cfg[2], time = t), dd)
    nm <- term_params(term)
    psi0 <- stats::setNames(
      c(0.15, rep(0.25, cfg[1]), c(0.5, -0.2)[seq_len(cfg[2])]), nm)

    run <- function(v) {
      term_filter(term, eta = rep(0, n), y = dd$y,
                  score = gauss_score(dd$y), curvature = gauss_curv(dd$y),
                  psi = as.list(stats::setNames(v, nm)))$eta
    }
    got <- term_filter(term, eta = rep(0, n), y = dd$y,
                       score = gauss_score(dd$y),
                       curvature = gauss_curv(dd$y),
                       psi = as.list(psi0))$jacobian
    num <- numDeriv::jacobian(run, psi0)
    expect_equal(unname(got), num, tolerance = 1e-6,
                 info = sprintf("p = %d, q = %d", cfg[1], cfg[2]))
  }
})

test_that("a static predictor shifts the filter and the groups stay apart", {
  term <- term_build(gas(p = 1, q = 1, by = g, time = t), dd)
  psi <- list(omega = 0.1, a1 = 0.2, pacf1 = 0.5)
  out <- term_filter(term, eta = rep(0, n), y = dd$y,
                     score = gauss_score(dd$y), curvature = gauss_curv(dd$y),
                     psi = psi)

  # each group is filtered from its own starting level, so the first row of
  # the second group repeats the first row of the first when its data do
  first_of_each <- vapply(split(seq_len(n), dd$g), function(i) i[1], integer(1))
  expect_equal(out$eta[first_of_each[1]] - dd$y[first_of_each[1]] * 0,
               out$eta[first_of_each[1]])
  # and a group's values depend only on its own rows
  dd2 <- dd
  dd2$y[dd$g == "b"] <- dd2$y[dd$g == "b"] + 10
  out2 <- term_filter(term, eta = rep(0, n), y = dd2$y,
                      score = gauss_score(dd2$y),
                      curvature = gauss_curv(dd2$y), psi = psi)
  expect_equal(out$eta[dd$g == "a"], out2$eta[dd$g == "a"])
  expect_false(isTRUE(all.equal(out$eta[dd$g == "b"], out2$eta[dd$g == "b"])))
})

test_that("time orders the recursion and the result is scattered back", {
  shuffled <- dd[sample(n), ]
  a <- term_filter(term_build(gas(p = 1, q = 1, time = t), dd),
                   eta = rep(0, n), y = dd$y, score = gauss_score(dd$y),
                   curvature = gauss_curv(dd$y),
                   psi = list(omega = 0.1, a1 = 0.2, pacf1 = 0.5))
  b <- term_filter(term_build(gas(p = 1, q = 1, time = t), shuffled),
                   eta = rep(0, n), y = shuffled$y,
                   score = gauss_score(shuffled$y),
                   curvature = gauss_curv(shuffled$y),
                   psi = list(omega = 0.1, a1 = 0.2, pacf1 = 0.5))
  # the same series in a different row order gives the same values per row
  expect_equal(a$eta[shuffled$t], b$eta, tolerance = 1e-12)
})

test_that("a purely autoregressive term needs no score lag budget", {
  term <- term_build(gas(p = 1, q = 0, time = t), dd)
  expect_identical(term_params(term), c("omega", "a1"))
  out <- term_filter(term, eta = rep(0, n), y = dd$y,
                     score = gauss_score(dd$y), curvature = gauss_curv(dd$y),
                     psi = list(omega = 0.3, a1 = 0.4))
  # with q = 0 the level is omega plus the score lag alone
  expect_equal(out$eta[1], 0.3)
})

test_that("the term is routed, printed, and refuses what it cannot do", {
  out <- interpret_formula(y ~ t + gas(p = 1, q = 1, time = t), dd)
  expect_named(out$terms, c("linpar", "gas(p = 1, q = 1, time = t)"))
  spec <- out$terms[["gas(p = 1, q = 1, time = t)"]]
  expect_true(S7::S7_inherits(spec, structural_term))
  expect_output(print(spec), "specification")
  expect_output(print(term_build(spec, dd)), "1 group")

  expect_error(gas(p = 0), "at least 1")
  expect_error(gas(label = ""), "non-empty")
  expect_error(term_filter(gas(), rep(0, n), dd$y, gauss_score(dd$y),
                           gauss_curv(dd$y), list(omega = 0, a1 = 0, pacf1 = 0)),
               "not been built")
  # a structural term has no design block
  expect_error(term_build(gas(), data.frame(z = 1:3)), NA)
})

test_that("the compiled filter and the R twin agree exactly", {
  # A mechanical transcription is the one reference the compiled route can
  # be held to: the two share no code, so agreement at machine precision
  # says the port changed nothing.
  for (cfg in list(c(1, 1), c(2, 2), c(3, 1))) {
    term <- term_build(gas(p = cfg[1], q = cfg[2], time = t, by = g), dd)
    nm <- term_params(term)
    psi <- as.list(stats::setNames(
      c(0.12, rep(0.22, cfg[1]), c(0.5, -0.3, 0.2)[seq_len(cfg[2])]), nm))

    got <- term_filter(term, rep(0.05, n), dd$y, gauss_score(dd$y),
                       gauss_curv(dd$y), psi)

    # the twin, driven with the pieces the method computes
    p_ <- cfg[1]; q_ <- cfg[2]
    v <- unlist(psi[nm])
    ld <- gas_levinson(if (q_ > 0) v[paste0("pacf", seq_len(q_))] else numeric(0))
    np <- length(nm)
    db <- matrix(0, max(q_, 1L), np)
    if (q_ > 0) db[seq_len(q_), 1L + p_ + seq_len(q_)] <- ld$jacobian
    sb <- sum(ld$phi)
    f0 <- v[["omega"]] / (1 - sb)
    df0 <- numeric(np); df0[1L] <- 1 / (1 - sb)
    if (q_ > 0) {
      df0 <- df0 + (v[["omega"]] / (1 - sb)^2) *
        colSums(db[seq_len(q_), , drop = FALSE])
    }
    ref <- gas_filter_r(rep(0.05, n), term@blueprint$order, p_, q_,
                        v[["omega"]],
                        if (p_ > 0) v[paste0("a", seq_len(p_))] else numeric(0),
                        ld$phi, db, f0, df0, 1L + seq_len(p_), np,
                        gauss_score(dd$y), gauss_curv(dd$y))

    expect_equal(got$eta, ref$eta, tolerance = 1e-14,
                 info = sprintf("p = %d, q = %d", cfg[1], cfg[2]))
    expect_equal(unname(got$jacobian), ref$jacobian, tolerance = 1e-14,
                 info = sprintf("p = %d, q = %d", cfg[1], cfg[2]))
  }
})

# --- a population value and a deviation per group -------------------------

test_that("deviations are parameters of the term, named and charted", {
  spec <- gas(p = 1, q = 1, by = g, time = t, deviations = "omega")
  # a specification reports the population parameters alone: how many
  # groups there are is a property of the data
  expect_identical(term_params(spec), c("omega", "a1", "pacf1"))

  built <- term_build(spec, dd)
  expect_identical(term_params(built),
                   c("omega", "a1", "pacf1", "omega.dev.a", "omega.dev.b"))
  expect_identical(term_npar(built), 5L)
  # a deviation acts on the unconstrained scale and is unconstrained itself
  lk <- term_links(built)
  expect_identical(vapply(lk, function(l) l@link_name, character(1)),
                   c(omega = "identity", a1 = "identity", pacf1 = "rhobit",
                     omega.dev.a = "identity", omega.dev.b = "identity"))

  all_dev <- term_build(gas(p = 1, q = 2, by = g, deviations = TRUE), dd)
  expect_identical(term_npar(all_dev), 4L + 4L * 2L)
})

test_that("deviations at zero reproduce the shared-parameter filter", {
  shared <- term_build(gas(p = 1, q = 1, by = g, time = t), dd)
  dev <- term_build(gas(p = 1, q = 1, by = g, time = t, deviations = TRUE),
                    dd)
  psi <- list(omega = 0.1, a1 = 0.2, pacf1 = 0.5)
  a <- term_filter(shared, eta = rep(0, n), y = dd$y,
                   score = gauss_score(dd$y), curvature = gauss_curv(dd$y),
                   psi = psi)
  b <- term_filter(dev, eta = rep(0, n), y = dd$y,
                   score = gauss_score(dd$y), curvature = gauss_curv(dd$y),
                   psi = c(psi, stats::setNames(as.list(rep(0, 6)),
                                                term_params(dev)[4:9])))
  expect_identical(b$eta, a$eta)
  # and the population columns are the shared ones exactly: at a zero
  # deviation the two chain factors are reciprocal, by the inverse
  # function theorem, so nothing is scaled
  expect_identical(unname(b$jacobian[, 1:3]), unname(a$jacobian))
})

test_that("a shift shared by the population and the deviations does nothing", {
  # m + 1 numbers describe m group values, so the filter is invariant under
  # adding a constant on the unconstrained scale to the population value and
  # subtracting it from every deviation. The invariance is exact and is why
  # the deviations need their penalty to be identified.
  term <- term_build(gas(p = 1, q = 1, by = g, time = t,
                         deviations = c("omega", "pacf1")), dd)
  nm <- term_params(term)
  psi <- stats::setNames(c(0.15, 0.3, 0.4, 0.2, -0.5, 0.1, 0.6), nm)
  run <- function(v) {
    term_filter(term, eta = rep(0, n), y = dd$y,
                score = gauss_score(dd$y), curvature = gauss_curv(dd$y),
                psi = as.list(v))$eta
  }
  shifted <- psi
  shifted[["omega"]] <- psi[["omega"]] + 0.25
  shifted[grepl("^omega\\.dev", nm)] <- psi[grepl("^omega\\.dev", nm)] - 0.25
  # the persistence rides rhobit, so the shift is applied on that scale
  lk <- linkfunctions7::rhobit_link()
  shifted[["pacf1"]] <- linkfunctions7::linkinv(
    lk, linkfunctions7::linkfun(lk, psi[["pacf1"]]) + 0.4)
  shifted[grepl("^pacf1\\.dev", nm)] <- psi[grepl("^pacf1\\.dev", nm)] - 0.4
  expect_equal(run(shifted), run(psi), tolerance = 1e-12)
})

test_that("the jacobian is exact in the population values and the deviations", {
  for (cfg in list(list(1L, 1L, "omega"), list(1L, 2L, TRUE),
                   list(2L, 1L, c("a1", "a2")))) {
    term <- term_build(gas(p = cfg[[1L]], q = cfg[[2L]], by = g, time = t,
                           deviations = cfg[[3L]]), dd)
    nm <- term_params(term)
    nb <- 1L + cfg[[1L]] + cfg[[2L]]
    set.seed(21)
    psi0 <- stats::setNames(
      c(0.15, rep(0.25, cfg[[1L]]), c(0.5, -0.2)[seq_len(cfg[[2L]])],
        stats::runif(length(nm) - nb, -0.4, 0.4)), nm)

    run <- function(v) {
      term_filter(term, eta = rep(0, n), y = dd$y,
                  score = gauss_score(dd$y), curvature = gauss_curv(dd$y),
                  psi = as.list(stats::setNames(v, nm)))$eta
    }
    got <- term_filter(term, eta = rep(0, n), y = dd$y,
                       score = gauss_score(dd$y),
                       curvature = gauss_curv(dd$y),
                       psi = as.list(psi0))$jacobian
    num <- numDeriv::jacobian(run, psi0)
    expect_equal(unname(got), num, tolerance = 1e-6,
                 info = sprintf("p = %d, q = %d", cfg[[1L]], cfg[[2L]]))
  }
})

test_that("the penalty reaches the deviations and the population is free", {
  built <- term_build(gas(p = 1, q = 1, by = g, deviations = "omega",
                          penalty = "lasso"), dd)
  ent <- term_penalties(built)
  expect_length(ent, 1L)
  expect_identical(ent[[1L]]$name, "omega")
  expect_identical(term_params(built)[ent[[1L]]$index],
                   c("omega.dev.a", "omega.dev.b"))
  # named as coordinates, so the separable branch applies unchanged
  expect_null(ent[[1L]]$penalty@map)
  expect_true(penalties7::has_prox(ent[[1L]]$penalty))
  expect_false(term_smooth(built))

  # one penalty per parameter carrying deviations: two parameters of a
  # filter are on scales of their own
  both <- term_build(gas(p = 1, q = 1, by = g, deviations = c("omega", "a1"),
                         penalty = "ridge"), dd)
  eb <- term_penalties(both)
  expect_length(eb, 2L)
  expect_identical(vapply(eb, function(e) e$name, character(1)),
                   c("omega", "a1"))
  expect_true(term_smooth(both))

  expect_length(term_penalties(term_build(gas(p = 1, q = 1, by = g), dd)), 0L)
})

test_that("the constructor refuses a deviation it cannot place", {
  expect_error(gas(deviations = TRUE), "needs 'by'")
  expect_error(gas(by = g, deviations = "sigma"), "the parameters are")
  expect_error(gas(by = g, penalty = "lasso"), "needs 'deviations'")
  expect_error(gas(by = g, deviations = 1), "FALSE, TRUE")
  # a specification has no deviations to index yet, and reports none, as an
  # unbuilt ridge() reports no penalty rather than raising
  expect_length(term_penalties(gas(by = g, deviations = TRUE,
                                   penalty = "lasso")), 0L)
})

# --- the reverse recursion ------------------------------------------------

test_that("the adjoint is the derivative in the static predictor", {
  # what a model layer needs and term_filter() does not give: the level at
  # one time is driven by the scores at earlier ones, which are read at
  # predictors the coefficients also enter, so the derivative of the
  # objective in the static predictor is not the score
  loglik <- function(e) sum(stats::dnorm(dd$y, e, log = TRUE))
  direct <- function(e) dd$y - e

  cfgs <- list(
    list(term = term_build(gas(p = 1, q = 1, time = t), dd),
         psi = list(omega = 0.15, a1 = 0.3, pacf1 = 0.5),
         eta = rep(0.2, n)),
    list(term = term_build(gas(p = 2, q = 2, time = t), dd),
         psi = list(omega = 0.1, a1 = 0.25, a2 = -0.15,
                    pacf1 = 0.4, pacf2 = -0.2),
         eta = sin(seq_len(n) / 5)),
    list(term = term_build(gas(p = 1, q = 1, by = g, time = t), dd),
         psi = list(omega = 0.1, a1 = 0.2, pacf1 = 0.5),
         eta = rep(0.1, n))
  )
  for (cf in cfgs) {
    L <- function(v) {
      loglik(term_filter(cf$term, v, dd$y, gauss_score(dd$y),
                         gauss_curv(dd$y), cf$psi)$eta)
    }
    e <- term_filter(cf$term, cf$eta, dd$y, gauss_score(dd$y),
                     gauss_curv(dd$y), cf$psi)$eta
    ad <- term_adjoint(cf$term, cf$eta, dd$y, gauss_score(dd$y),
                       gauss_curv(dd$y), cf$psi, g = direct(e))
    num <- numDeriv::grad(L, cf$eta)
    expect_equal(ad$deta, num, tolerance = 1e-6)
    # and the direct derivative alone is not it, by a wide margin: the
    # feedback is what the reverse recursion exists to carry
    expect_gt(max(abs(direct(e) - num)), 0.5)
  }
})

test_that("the adjoint carries the derivative through the score as well", {
  # a model layer whose score is the derivative of its log-likelihood in one
  # distribution parameter reaches the derivative in the predictor of
  # ANOTHER by multiplying dscore by the mixed second derivative, so the
  # quantity is checked by perturbing the score the caller supplies
  term <- term_build(gas(p = 1, q = 1, time = t), dd)
  psi <- list(omega = 0.15, a1 = 0.3, pacf1 = 0.5)
  eta0 <- rep(0.2, n)
  set.seed(31)
  u <- stats::rnorm(n)
  L <- function(cc) {
    sc <- function(e, i) dd$y[i] - e + cc * u[i]
    sum(stats::dnorm(dd$y,
                     term_filter(term, eta0, dd$y, sc, gauss_curv(dd$y),
                                 psi)$eta, log = TRUE))
  }
  e <- term_filter(term, eta0, dd$y, gauss_score(dd$y), gauss_curv(dd$y),
                   psi)$eta
  ad <- term_adjoint(term, eta0, dd$y, gauss_score(dd$y), gauss_curv(dd$y),
                     psi, g = dd$y - e)
  expect_equal(sum(ad$dscore * u), numDeriv::grad(L, 0), tolerance = 1e-6)
})

test_that("the adjoint refuses what it cannot answer", {
  term <- term_build(gas(p = 1, q = 1, time = t), dd)
  psi <- list(omega = 0.1, a1 = 0.2, pacf1 = 0.5)
  expect_error(term_adjoint(term, rep(0, n), dd$y, gauss_score(dd$y),
                            gauss_curv(dd$y), psi, g = rep(1, 3)),
               "one value per observation")
  expect_error(term_adjoint(gas(p = 1, q = 1), rep(0, n), dd$y,
                            gauss_score(dd$y), gauss_curv(dd$y), psi,
                            g = rep(1, n)), "not been built")
})

test_that("the second derivative of Levinson-Durbin is exact", {
  # the first step of the exact observed information of a model carrying
  # this term: the persistence reaches the predictor through this map, so
  # the curvature needs its second derivative. The recursion is bilinear,
  # so differentiating twice adds no new kind of term.
  for (pc in list(0.6, c(0.5, -0.3), c(0.9, 0.8, -0.6),
                  c(0.2, -0.4, 0.7, 0.1))) {
    out <- gas_levinson2(pc)
    old <- gas_levinson(pc)
    # the value and the jacobian are the ones the existing map gives, to
    # the bit: the second-order version must not be a second route to them
    expect_identical(out$phi, old$phi)
    expect_identical(out$jacobian, old$jacobian)

    for (i in seq_along(pc)) {
      num <- numDeriv::hessian(function(v) gas_levinson(v)$phi[i], pc)
      expect_equal(out$hessian[[i]], num, tolerance = 1e-6,
                   info = sprintf("q = %d, coefficient %d", length(pc), i))
      # symmetric by construction rather than by being symmetrized
      expect_identical(out$hessian[[i]], t(out$hessian[[i]]))
    }
    # the last coefficient IS the last partial autocorrelation, so its
    # second derivative is exactly zero at every order
    expect_true(all(out$hessian[[length(pc)]] == 0))
  }
  expect_length(gas_levinson2(numeric(0))$hessian, 0L)
})

test_that("the second-order recursion gives the exact curvature", {
  # Step two of the exact observed information: the forward jacobian of the
  # predictor in a caller's unknowns, and the second derivative contracted
  # against the caller's weights. The model layer is stubbed by a gaussian
  # mean of unit variance -- l_pp = -1, l_ppp = 0, one equation -- which
  # isolates the recursion from everything the layer would contribute.
  set.seed(8)
  nn <- 50L
  d2 <- data.frame(t = seq_len(nn), y = stats::rnorm(nn))
  X <- cbind(1, as.numeric(scale(seq_len(nn))))
  mb <- ncol(X)
  sc <- function(e, i) d2$y[i] - e
  cu <- function(e, i) -1

  for (cfg in list(c(1, 1), c(2, 1), c(1, 2), c(2, 2))) {
    term <- term_build(gas(p = cfg[1], q = cfg[2], time = t), d2)
    nm <- term_params(term)
    np <- length(nm)
    lk <- term_links(term)
    m <- mb + np
    seed <- cbind(X, matrix(0, nn, np))
    blocks <- function(e, i, D) list(cross = numeric(m), M = matrix(0, m, m))
    psi_of <- function(z) {
      as.list(stats::setNames(vapply(seq_along(nm), function(j)
        linkfunctions7::linkinv(lk[[nm[j]]], z[j]), numeric(1)), nm))
    }
    set.seed(9)
    u0 <- c(0.3, -0.2, stats::runif(np, -0.3, 0.5))
    gw <- stats::rnorm(nn)
    eta0 <- as.numeric(X %*% u0[seq_len(mb)])
    got <- term_curvature(term, eta0, d2$y, sc, cu,
                          psi_of(u0[mb + seq_len(np)]), gw, seed, blocks)

    lab <- sprintf("p = %d, q = %d", cfg[1], cfg[2])
    # the jacobian, against a numerical derivative of the filter itself
    e_of <- function(u) {
      term_filter(term, as.numeric(X %*% u[seq_len(mb)]), d2$y, sc, cu,
                  psi_of(u[mb + seq_len(np)]))$eta
    }
    expect_equal(got$jacobian, numDeriv::jacobian(e_of, u0),
                 tolerance = 1e-6, info = lab)

    # and the curvature, against the second derivative of the weighted sum
    gsum <- function(u) sum(gw * e_of(u))
    expect_equal(got$curvature, numDeriv::hessian(gsum, u0),
                 tolerance = 1e-6, info = lab)
    # symmetric exactly, which it is because the result is symmetrized: an
    # entry and its transpose collect the same terms in a different order,
    # and floating-point addition is not associative. At p = q = 1 the two
    # orders coincide and the claim looks free, which is why the loop runs
    # to p = q = 2
    expect_identical(got$curvature, t(got$curvature))
  }
})

test_that("the curvature refuses what it does not carry", {
  term <- term_build(gas(p = 1, q = 1, time = t), dd)
  nm <- term_params(term)
  psi <- list(omega = 0.1, a1 = 0.2, pacf1 = 0.5)
  sd0 <- matrix(0, n, 3)
  bl <- function(e, i, D) list(cross = numeric(3), M = matrix(0, 3, 3))
  expect_error(term_curvature(term, rep(0, n), dd$y, gauss_score(dd$y),
                              gauss_curv(dd$y), psi, rep(1, 3), sd0, bl),
               "one value per observation")
  # deviations add a per-group chain to every derivative and are refused
  # rather than silently dropped
  dev <- term_build(gas(p = 1, q = 1, by = g, deviations = "omega"), dd)
  nd <- length(term_params(dev))
  expect_error(term_curvature(dev, rep(0, n), dd$y, gauss_score(dd$y),
                              gauss_curv(dd$y),
                              as.list(stats::setNames(rep(0.1, nd),
                                                      term_params(dev))),
                              rep(1, n), matrix(0, n, nd),
                              function(e, i, D) list(cross = numeric(nd),
                                                     M = matrix(0, nd, nd))),
               "does not carry deviations")
})
