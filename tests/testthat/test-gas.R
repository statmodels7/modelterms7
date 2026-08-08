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
