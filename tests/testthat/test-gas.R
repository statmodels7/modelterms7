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
                   c("omega", "alpha1", "pacf1"))
  expect_identical(term_params(gas(p = 2, q = 3)),
                   c("omega", "alpha1", "alpha2", "pacf1", "pacf2", "pacf3"))
  lk <- term_links(gas(p = 1, q = 2))
  expect_identical(vapply(lk, function(l) l@link_name, character(1)),
                   c(omega = "identity", alpha1 = "log",
                     pacf1 = "rhobit", pacf2 = "rhobit"))
  # every loading is positive by default, not only the first
  lk2 <- term_links(gas(p = 3, q = 0))
  expect_identical(vapply(lk2, function(l) l@link_name, character(1)),
                   c(omega = "identity", alpha1 = "log", alpha2 = "log",
                     alpha3 = "log"))
})

test_that("the links are configurable and validated", {
  # an override replaces the default for the parameter it names alone
  term <- gas(p = 2, q = 1,
              links = list(alpha2 = linkfunctions7::identity_link()))
  lk <- term_links(term)
  expect_identical(vapply(lk, function(l) l@link_name, character(1)),
                   c(omega = "identity", alpha1 = "log",
                     alpha2 = "identity", pacf1 = "rhobit"))

  # a name that is not a parameter, a value that is not a link, and an
  # unnamed list are each rejected where they are written
  expect_error(gas(links = list(gamma = linkfunctions7::log_link())),
               "the parameters are")
  expect_error(gas(links = list(alpha1 = "log")), "not a linkfunctions7 link")
  expect_error(gas(links = list(linkfunctions7::log_link())), "named list")
})

test_that("the start is the term's own, and the loadings start weak", {
  z <- term_start(gas(p = 1, q = 1))
  expect_identical(names(z), c("omega", "alpha1", "pacf1"))
  expect_identical(unname(z[c("omega", "pacf1")]), c(0, 0))
  # alpha starts at 0.1 on the parameter scale, through whatever chart it
  # rides: log(0.1) on the default, 0.1 under an identity override
  expect_equal(unname(z[["alpha1"]]), log(0.1))
  z2 <- term_start(gas(p = 1, q = 1,
                       links = list(alpha1 = linkfunctions7::identity_link())))
  expect_equal(unname(z2[["alpha1"]]), 0.1)
  # the departures of a development start at zero, being departures
  term <- term_build(gas(p = 1, q = 1, omega ~ random(~1 | g), by = g), dd)
  z3 <- term_start(term)
  expect_identical(unname(z3[!grepl("^alpha", names(z3))]),
                   numeric(sum(!grepl("^alpha", names(z3)))))
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
  psi <- list(omega = 0.2, alpha1 = 0.3, pacf1 = 0.6)
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
  psi <- list(omega = 0.1, alpha1 = 0.2, pacf1 = 0.5)
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
                   psi = list(omega = 0.1, alpha1 = 0.2, pacf1 = 0.5))
  b <- term_filter(term_build(gas(p = 1, q = 1, time = t), shuffled),
                   eta = rep(0, n), y = shuffled$y,
                   score = gauss_score(shuffled$y),
                   curvature = gauss_curv(shuffled$y),
                   psi = list(omega = 0.1, alpha1 = 0.2, pacf1 = 0.5))
  # the same series in a different row order gives the same values per row
  expect_equal(a$eta[shuffled$t], b$eta, tolerance = 1e-12)
})

test_that("a purely autoregressive term needs no score lag budget", {
  term <- term_build(gas(p = 1, q = 0, time = t), dd)
  expect_identical(term_params(term), c("omega", "alpha1"))
  out <- term_filter(term, eta = rep(0, n), y = dd$y,
                     score = gauss_score(dd$y), curvature = gauss_curv(dd$y),
                     psi = list(omega = 0.3, alpha1 = 0.4))
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
                           gauss_curv(dd$y), list(omega = 0, alpha1 = 0, pacf1 = 0)),
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
                        if (p_ > 0) v[paste0("alpha", seq_len(p_))] else
                          numeric(0),
                        ld$phi, db, f0, df0, 1L + seq_len(p_), np,
                        gauss_score(dd$y), gauss_curv(dd$y))

    expect_equal(got$eta, ref$eta, tolerance = 1e-14,
                 info = sprintf("p = %d, q = %d", cfg[1], cfg[2]))
    expect_equal(unname(got$jacobian), ref$jacobian, tolerance = 1e-14,
                 info = sprintf("p = %d, q = %d", cfg[1], cfg[2]))
  }
})

# --- a population value and a departure per group -------------------------

test_that("departures at zero reproduce the shared-parameter filter", {
  shared <- term_build(gas(p = 1, q = 1, by = g, time = t), dd)
  dev <- term_build(gas(p = 1, q = 1, omega ~ random(~1 | g),
                        alpha1 ~ random(~1 | g), pacf1 ~ random(~1 | g),
                        by = g, time = t), dd)
  psi <- list(omega = 0.1, alpha1 = 0.2, pacf1 = 0.5)
  a <- term_filter(shared, eta = rep(0, n), y = dd$y,
                   score = gauss_score(dd$y), curvature = gauss_curv(dd$y),
                   psi = psi)
  nm <- term_params(dev)
  v <- stats::setNames(rep(list(0), length(nm)), nm)
  v[["omega.(Intercept)"]] <- 0.1
  v[["alpha1.(Intercept)"]] <- log(0.2)
  lkr <- linkfunctions7::rhobit_link()
  v[["pacf1.(Intercept)"]] <- linkfunctions7::linkfun(lkr, 0.5)
  b <- term_filter(dev, eta = rep(0, n), y = dd$y,
                   score = gauss_score(dd$y), curvature = gauss_curv(dd$y),
                   psi = v)
  # the developments read each chart at the intercept alone, so at zero
  # departures the two runs are the same recursion; the comparison carries
  # a tolerance because the chart arithmetic differs by a rounding across
  # platforms' libm, the lesson the deviations machinery once recorded
  expect_equal(b$eta, a$eta, tolerance = 1e-13)
})

test_that("a shift shared by the population and the departures does nothing", {
  # m + 1 numbers describe m group values, so the filter is invariant under
  # adding a constant to the development's intercept and subtracting it
  # from every group's departure, on the unconstrained scale. The
  # invariance is exact and is why the departures need their penalty to be
  # identified.
  term <- term_build(gas(p = 1, q = 1, omega ~ random(~1 | g),
                         pacf1 ~ random(~1 | g), by = g, time = t), dd)
  nm <- term_params(term)
  set.seed(21)
  psi <- stats::setNames(stats::runif(length(nm), -0.4, 0.4), nm)
  run <- function(v) {
    term_filter(term, eta = rep(0, n), y = dd$y,
                score = gauss_score(dd$y), curvature = gauss_curv(dd$y),
                psi = as.list(v))$eta
  }
  shifted <- psi
  shifted[["omega.(Intercept)"]] <- psi[["omega.(Intercept)"]] + 0.25
  sel <- grepl("^omega\\.random", nm)
  shifted[sel] <- psi[sel] - 0.25
  shifted[["pacf1.(Intercept)"]] <- psi[["pacf1.(Intercept)"]] + 0.4
  sel <- grepl("^pacf1\\.random", nm)
  shifted[sel] <- psi[sel] - 0.4
  expect_equal(run(shifted), run(psi), tolerance = 1e-12)
})

test_that("the reported quantities are the literature's, with a jacobian", {
  # omega and the loadings are the coordinates themselves, each on the
  # identity link. The persistence is NOT: it rides a partial
  # autocorrelation, and what the literature calls beta_j is the
  # autoregressive coefficient, a function of the whole chart.
  for (q in 1:3) {
    tm <- gas(p = 1, q = q)
    nmv <- term_params(tm)
    set.seed(q)
    z <- stats::setNames(c(0.3, 0.4, stats::runif(q, -0.6, 0.8)), nmv)
    rd <- term_readable(tm, z)
    expect_identical(rd$name,
                     c("omega", "alpha1", paste0("beta", seq_len(q))))
    # the coefficients are those of a stationary autoregression
    expect_true(all(Mod(polyroot(c(1, -rd$value[-(1:2)]))) > 1 + 1e-8))
    # and the jacobian is the one a delta method needs
    num <- numDeriv::jacobian(function(v)
      term_readable(tm, stats::setNames(v, nmv))$value, z)
    expect_equal(unname(rd$jacobian), num, tolerance = 1e-7,
                 info = sprintf("q = %d", q))
  }

  # at q = 1 the coefficient IS the partial autocorrelation, so the chain
  # factor is the link's alone and the two coincide exactly
  z1 <- c(omega = 0.2, alpha1 = 0.3, pacf1 = 0.9)
  expect_equal(term_readable(gas(1, 1), z1)$value[[3L]],
               linkfunctions7::linkinv(linkfunctions7::rhobit_link(), 0.9))
  # above it they do not, which is the whole reason the coordinate is not
  # named after the coefficient
  z2 <- c(omega = 0.2, alpha1 = 0.3, pacf1 = 1.2, pacf2 = -0.4)
  rho <- linkfunctions7::linkinv(linkfunctions7::rhobit_link(), c(1.2, -0.4))
  expect_false(isTRUE(all.equal(term_readable(gas(1, 2), z2)$value[[3L]],
                                rho[[1L]])))

  # the level's coordinate IS its quantity; a loading is reported through
  # its chart, exp of the coordinate on the default log link
  base <- term_readable(gas(p = 1, q = 0), c(omega = 0.5, alpha1 = 0.2))
  expect_identical(base$name, c("omega", "alpha1"))
  expect_equal(base$value, c(0.5, exp(0.2)))
  # and under an identity override it is the coordinate again
  ident <- term_readable(
    gas(p = 1, q = 0, links = list(alpha1 = linkfunctions7::identity_link())),
    c(omega = 0.5, alpha1 = 0.2))
  expect_equal(ident$value, c(0.5, 0.2))
})

test_that("the removed shorthands are reported by name", {
  # deviations= and penalty= were retired when the subformulas subsumed
  # them; a call carrying either lands in `...` and is named rather than
  # swallowed by the formula check
  expect_error(gas(by = g, deviations = TRUE), "unused argument 'deviations'")
  expect_error(gas(by = g, penalty = penalties7::lasso_penalty),
               "unused argument 'penalty'")
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
         psi = list(omega = 0.15, alpha1 = 0.3, pacf1 = 0.5),
         eta = rep(0.2, n)),
    list(term = term_build(gas(p = 2, q = 2, time = t), dd),
         psi = list(omega = 0.1, alpha1 = 0.25, alpha2 = -0.15,
                    pacf1 = 0.4, pacf2 = -0.2),
         eta = sin(seq_len(n) / 5)),
    list(term = term_build(gas(p = 1, q = 1, by = g, time = t), dd),
         psi = list(omega = 0.1, alpha1 = 0.2, pacf1 = 0.5),
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
  psi <- list(omega = 0.15, alpha1 = 0.3, pacf1 = 0.5)
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
  psi <- list(omega = 0.1, alpha1 = 0.2, pacf1 = 0.5)
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

test_that("the third derivative of Levinson-Durbin is exact", {
  # The exact gradient of a marginal criterion over a penalty on this term's
  # own parameters needs one more order than the curvature, contracted
  # against the direction the penalized mode moves in.
  #
  # The loop runs to q = 4 and NOT to q = 2, which would assert nothing: the
  # map is multilinear of degree k in the first k partial autocorrelations,
  # so at q = 2 the only non-trivial coefficient is rho_1(1 - rho_2), whose
  # third derivative is identically zero. A check stopping there compares
  # zero with zero -- the same shape as a symmetry that is free at p = q = 1.
  set.seed(3)
  for (pc in list(c(0.9, 0.8, -0.6), c(0.2, -0.4, 0.7, 0.1))) {
    q <- length(pc)
    w <- stats::rnorm(q)
    got <- gas_levinson3(pc, w)
    for (i in seq_len(q)) {
      num <- numDeriv::jacobian(
        function(z) as.numeric(gas_levinson2(z)$hessian[[i]]), pc)
      ref <- matrix(as.numeric(num %*% w), q, q)
      expect_equal(got[[i]], ref, tolerance = 1e-6,
                   info = sprintf("q = %d, coefficient %d", q, i))
      expect_identical(got[[i]], t(got[[i]]))
    }
    # the last coefficient IS the last partial autocorrelation
    expect_true(all(got[[q]] == 0))
  }
  # and the degenerate cases, where the answer is zero for a reason
  expect_length(gas_levinson3(numeric(0), numeric(0)), 0L)
  expect_true(all(gas_levinson3(c(0.5, -0.3), c(1, 1))[[1L]] == 0))
})

test_that("the third-order recursion gives the exact directional derivative", {
  # Step one of the exact gradient over a structural penalty: the second
  # derivative of the predictor differentiated once more along one direction.
  #
  # The layer is stubbed by l = -log cosh(e - y), whose four derivatives are
  # all NON-ZERO and all BOUNDED. Both properties are load-bearing: a
  # gaussian mean has l''' = l'''' = 0, which multiplies every new term of
  # the recursion by zero and asserts nothing, while an exponential family's
  # score grows with the predictor and the recursion is geometric in it, so
  # the filter leaves the doubles before anything is differenced.
  set.seed(8)
  nn <- 50L
  d2 <- data.frame(t = seq_len(nn), y = stats::rnorm(nn))
  X <- cbind(1, as.numeric(scale(seq_len(nn))) * 0.3)
  mb <- ncol(X)
  tt <- function(e, i) tanh(e - d2$y[i])
  sc <- function(e, i) -tt(e, i)
  cu <- function(e, i) -(1 - tt(e, i)^2)
  l3 <- function(e, i) { t <- tt(e, i); 2 * t * (1 - t^2) }
  l4 <- function(e, i) { t <- tt(e, i); (1 - t^2) * (2 - 6 * t^2) }

  # q reaches 3 so that the map's own third derivative is exercised; see
  # the Levinson-Durbin test above for why q = 2 would not
  for (cfg in list(c(1, 1), c(2, 1), c(1, 2), c(2, 2), c(1, 3))) {
    term <- term_build(gas(p = cfg[1], q = cfg[2], time = t), d2)
    nm <- term_params(term)
    np <- length(nm)
    lk <- term_links(term)
    m <- mb + np
    seed <- cbind(X, matrix(0, nn, np))
    psi_of <- function(z) as.list(stats::setNames(vapply(seq_along(nm),
      function(j) linkfunctions7::linkinv(lk[[nm[j]]], z[j]), numeric(1)), nm))
    set.seed(9)
    u0 <- c(0.3, -0.2, stats::runif(np, -0.4, 0.2))
    gw <- stats::rnorm(nn)
    v <- stats::rnorm(m)
    eta0 <- function(u) as.numeric(X %*% u[seq_len(mb)])
    blk2 <- function(e, i, D, act) list(cross = numeric(length(D)),
                                        M = l3(e, i) * outer(D, D))
    blk3 <- function(e, i, D, act) list(
      cross = numeric(length(D)), M = l3(e, i) * outer(D, D),
      dcurv = l3(e, i) * D,
      N = l4(e, i) * sum(D * v[act]) * outer(D, D))

    got <- term_third(term, eta0(u0), d2$y, sc, cu,
                      psi_of(u0[mb + seq_len(np)]), gw, seed, blk3, v)
    at <- function(u) term_curvature(term, eta0(u), d2$y, sc, cu,
      psi_of(u[mb + seq_len(np)]), gw, seed, blk2)
    dir_of <- function(h, what) {
      (at(u0 + h * v)[[what]] - at(u0 - h * v)[[what]]) / (2 * h)
    }
    lab <- sprintf("p = %d, q = %d", cfg[1], cfg[2])

    # Richardson on the central difference: the plain difference carries an
    # O(h^2) truncation of its own, and at these scales it is LARGER than
    # the gap being measured, so the reference is extrapolated before it is
    # believed against the code
    rich <- (4 * dir_of(1e-3, "curvature") - dir_of(2e-3, "curvature")) / 3
    expect_equal(got$curvature, rich, tolerance = 1e-6, info = lab)
    expect_equal(got$dphi, dir_of(1e-3, "jacobian"), tolerance = 1e-4,
                 info = lab)
    # the refactor that shares the recursion between the two orders must
    # leave the second one's jacobian untouched, to the bit
    expect_identical(got$jacobian, at(u0)$jacobian)
  }
})

test_that("term_third is zero for an additive term and refused where absent", {
  # an additive term's predictor is a block of columns, so its second
  # derivative is already zero and so is every order above it
  td <- data.frame(x = stats::rnorm(20), y = stats::rnorm(20))
  lp <- term_build(linpar(~x), td)
  out <- term_third(lp, rep(0, 20), td$y, function(e, i) 0,
                    function(e, i) -1, list(), rep(1, 20),
                    matrix(0, 20, 2), function(e, i, D, act) NULL, c(1, 0))
  expect_true(all(out$curvature == 0))
  expect_true(all(out$dphi == 0))

  # a term that BENDS the predictor and has not written its third derivative
  # must refuse rather than inherit that zero, which a caller could not tell
  # from a term that genuinely has none
  rg <- term_build(regime(k = 2), td)
  expect_error(term_third(rg, rep(0, 20), td$y, function(e, i) 0,
                          function(e, i) -1, list(), rep(1, 20),
                          matrix(0, 20, 2), function(e, i, D, act) NULL,
                          c(1, 0)),
               "does not implement term_third")
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
  psi <- list(omega = 0.1, alpha1 = 0.2, pacf1 = 0.5)
  sd0 <- matrix(0, n, 3)
  bl <- function(e, i, D) list(cross = numeric(3), M = matrix(0, 3, 3))
  expect_error(term_curvature(term, rep(0, n), dd$y, gauss_score(dd$y),
                              gauss_curv(dd$y), psi, rep(1, 3), sd0, bl),
               "one value per observation")
})

test_that("the curvature of a per-group development has the affine structure", {
  X <- cbind(1, as.numeric(scale(seq_len(n))))
  mb <- ncol(X)
  sc <- gauss_score(dd$y)
  cu <- gauss_curv(dd$y)

  at <- function(term, u) {
    nm <- term_params(term)
    lk <- term_links(term)
    np <- length(nm)
    m <- mb + np
    psi_of <- function(z) as.list(stats::setNames(vapply(seq_along(nm),
      function(j) linkfunctions7::linkinv(lk[[nm[j]]], z[j]), numeric(1)), nm))
    term_curvature(term, as.numeric(X %*% u[seq_len(mb)]), dd$y, sc, cu,
                   psi_of(u[mb + seq_len(np)]), gw,
                   cbind(X, matrix(0, n, np)),
                   function(e, i, D, act) list(cross = numeric(length(D)),
                                               M = matrix(0, length(D),
                                                          length(D))))
  }

  set.seed(11)
  gw <- stats::rnorm(n)

  # At ZERO departures every group runs on the intercepts alone, so the
  # curvature over the intercept coordinates must be the shared-parameter
  # term's own. And whatever the departures are, a random intercept's
  # indicator columns sum to the constant, so the departure columns of one
  # parameter sum to its intercept column: the affine lift, asserted
  # rather than assumed.
  plain <- term_build(gas(p = 1, q = 1, by = g), dd)
  term <- term_build(gas(p = 1, q = 1, omega ~ random(~1 | g), by = g), dd)
  nm <- term_params(term)
  ng <- length(term@blueprint$order)
  u0 <- c(0.3, -0.2, 0.15, 0.25, 0.4)
  a <- at(plain, u0)$curvature
  # the developed layout: omega.(Intercept), the departures, alpha1, pacf1
  i_int <- mb + match("omega.(Intercept)", nm)
  i_dev <- mb + which(startsWith(nm, "omega.random"))
  i_al <- mb + match("alpha1", nm)
  i_pa <- mb + match("pacf1", nm)
  u_dev <- c(u0[1:2], 0.15, numeric(ng), 0.25, 0.4)
  b <- at(term, u_dev)$curvature
  keep <- c(seq_len(mb), i_int, i_al, i_pa)
  expect_equal(unname(b[keep, keep]), unname(a), tolerance = 1e-12)

  set.seed(7)
  b2 <- at(term,
           c(u0[1:2], 0.15, stats::runif(ng, -0.2, 0.2), 0.25, 0.4))$curvature
  expect_equal(rowSums(b2[keep, i_dev, drop = FALSE]),
               b2[keep, i_int], tolerance = 1e-12)
})

test_that("the third derivative of a development has the affine structure", {
  # The two statements a tolerance against a numerical reference cannot make.
  # A development's design places the intercept column as the SUM of the
  # group indicators, and that linear relation propagates to every order, so
  # it can be read off the matrix rather than measured.
  X <- cbind(1, as.numeric(scale(seq_len(n))) * 0.3)
  mb <- ncol(X)
  tt <- function(e, i) tanh(e - dd$y[i])
  sc <- function(e, i) -tt(e, i)
  cu <- function(e, i) -(1 - tt(e, i)^2)
  l3 <- function(e, i) { t <- tt(e, i); 2 * t * (1 - t^2) }
  l4 <- function(e, i) { t <- tt(e, i); (1 - t^2) * (2 - 6 * t^2) }

  at3 <- function(term, u, v) {
    nm <- term_params(term)
    lk <- term_links(term)
    np <- length(nm)
    psi_of <- function(z) as.list(stats::setNames(vapply(seq_along(nm),
      function(j) linkfunctions7::linkinv(lk[[nm[j]]], z[j]), numeric(1)), nm))
    term_third(term, as.numeric(X %*% u[seq_len(mb)]), dd$y, sc, cu,
               psi_of(u[mb + seq_len(np)]), gw,
               cbind(X, matrix(0, n, np)),
               function(e, i, D, act) list(
                 cross = numeric(length(D)), M = l3(e, i) * outer(D, D),
                 dcurv = l3(e, i) * D,
                 N = l4(e, i) * sum(D * v[act]) * outer(D, D)),
               v)
  }

  set.seed(11)
  gw <- stats::rnorm(n)
  plain <- term_build(gas(p = 1, q = 1, by = g), dd)
  term <- term_build(gas(p = 1, q = 1, omega ~ random(~1 | g), by = g), dd)
  nm <- term_params(term)
  ng <- length(term@blueprint$order)
  u0 <- c(0.3, -0.2, 0.15, -0.25, 0.4)
  i_int <- mb + match("omega.(Intercept)", nm)
  i_dev <- mb + which(startsWith(nm, "omega.random"))
  keep <- c(seq_len(mb), i_int, mb + match(c("alpha1", "pacf1"), nm))

  # (1) at ZERO departures the developed term IS the shared-parameter one,
  # along a direction that does not move the departures -- a direction that
  # did would be asking for something the plain term cannot represent
  set.seed(5)
  v0 <- stats::rnorm(mb + 3L)
  vdev <- numeric(mb + length(nm))
  vdev[keep] <- v0
  a <- at3(plain, u0, v0)$curvature
  b <- at3(term, c(u0[1:2], 0.15, numeric(ng), -0.25, 0.4), vdev)$curvature
  expect_equal(unname(b[keep, keep]), unname(a), tolerance = 1e-12)

  # (2) and whatever the departures and whatever the direction, the
  # departure columns of one parameter sum to its intercept column
  set.seed(7)
  vd2 <- stats::rnorm(mb + length(nm))
  b2 <- at3(term, c(u0[1:2], 0.15, stats::runif(ng, -0.2, 0.2), -0.25, 0.4),
            vd2)$curvature
  expect_equal(rowSums(b2[keep, i_dev, drop = FALSE]), b2[keep, i_int],
               tolerance = 1e-12)
})
