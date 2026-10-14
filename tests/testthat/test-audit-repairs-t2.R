# Repairs from the documentation audit of 0.95.1 (second tranche). Every block
# fails on 0.95.1.

set.seed(1)
bd <- data.frame(x = seq(0, 10, length.out = 200))
bd$y <- 1 + 0.5 * pmax(bd$x - 5, 0) + 2 * (bd$x > 5) + rnorm(200, 0, 0.3)

test_that("a marginal seg() integrates a narrow prior on the break-point", {
  md <- data.frame(id = rep(1:3, each = 8), x = rep(1:8, 3))
  set.seed(2)
  md$y <- rnorm(24, 0.5 * pmax(md$x - 4.5, 0), 0.4)
  s <- 0.4
  logdens <- function(e, i) dnorm(md$y[i], e, s, log = TRUE)
  score <- function(e, i) (md$y[i] - e) / s^2
  tm <- term_build(seg(x, psi ~ random(~ 1 | id), marginal = TRUE), md)
  # the group's marginal log-likelihood, integrated piece by piece between
  # the observations
  ref <- function(m, tau) {
    xs <- md$x[1:8]
    ys <- md$y[1:8]
    lc <- function(g) vapply(g, function(gg) {
      sum(dnorm(ys, 0.2 * xs + 0.5 * pmax(xs - gg, 0), s, log = TRUE))
    }, 0)
    c0 <- lc(m)
    br <- sort(unique(c(m - 40 * tau, xs[abs(xs - m) < 40 * tau],
                        m + 40 * tau)))
    tot <- 0
    for (j in seq_len(length(br) - 1L)) {
      tot <- tot + stats::integrate(function(g) {
        exp(lc(g) - c0) * dnorm(g, m, tau)
      }, br[j], br[j + 1L], rel.tol = 1e-11)$value
    }
    log(tot) + c0
  }
  psi <- list(beta = 0.2, psi1.mean = 2.2, psi1.sd = 0.01, gamma1 = 0.5)
  ll <- term_loglik(tm, rep(0, 24), md$y, logdens, score, psi)$loglik
  # the capped panels gave an error of 6.5e-4
  expect_lt(abs(sum(ll[1:8]) - ref(2.2, 0.01)), 1e-8)
  psi$psi1.sd <- 1e-5
  expect_error(term_loglik(tm, rep(0, 24), md$y, logdens, score, psi),
               "too small against the spacing")
})

test_that("a regime path is drawn one chain per group", {
  n <- 2000
  dd <- data.frame(id = seq_len(n), t = 1, y = 0)
  tm <- term_build(regime(2, by = id, time = t), dd)
  psi <- list(level1 = 0, gap2 = 5, alr1.1 = 6, alr2.1 = -6)
  set.seed(3)
  s <- term_simulate(tm, psi, rep(0, n), draw = FALSE)$latent
  # groups of one row are independent draws from the stationary law; one
  # chain over the rows agreed in 99.8 per cent of consecutive pairs
  expect_lt(abs(mean(s[-1] == s[-n]) - 0.5), 0.05)
})

test_that("a gas() continuation reads the starting level before the series", {
  dat <- data.frame(t = c(1, 1:6), id = c("A", rep("B", 6)),
                    y = c(0.5, 0.2, -0.1, 0.4, 0.3, -0.2, 0.1))
  tb <- term_build(gas(p = 1, q = 2, by = id, time = t), dat)
  psi <- list(omega = 0.2, kappa1 = 0.3, pacf1 = 0.4, pacf2 = -0.2)
  sc <- function(e, i) dat$y[i] - e
  fo <- term_filter(tb, rep(0, 7), dat$y, sc, function(e, i) -1, psi,
                    fast = NULL, threads = 1L)
  f <- as.numeric(fo$eta)
  newd <- data.frame(t = c(2, 3, 7), id = c("A", "A", "B"))
  cont <- as.numeric(term_continue(tb, psi, f, dat$y - f, newd))
  # the filter on the data with the new rows appended at a zero score
  aug <- rbind(dat, data.frame(t = newd$t, id = newd$id, y = NA))
  ta <- term_build(gas(p = 1, q = 2, by = id, time = t), aug)
  fa <- term_filter(ta, rep(0, 10), aug$y,
                    function(e, i) ifelse(i > 7, 0, aug$y[i] - e),
                    function(e, i) -1, psi, fast = NULL, threads = 1L)
  # the group with one observed row gave 0.4 where the filter gives 0.3444
  expect_equal(cont, as.numeric(fa$eta)[8:10], tolerance = 1e-14)
})

test_that("term_readable() of gas() applies each partial autocorrelation's link", {
  dd <- data.frame(y = rnorm(30))
  tg <- term_build(gas(p = 1, q = 2,
                       links = list(pacf2 = linkfunctions7::identity_link())),
                   dd)
  zeta <- list(omega = 0.1, kappa1 = 0.2, pacf1 = 0.5, pacf2 = 0.2)
  rd <- term_readable(tg, zeta)
  rho <- c(linkfunctions7::linkinv(linkfunctions7::rhobit_link(), 0.5), 0.2)
  expect_equal(rd$value[rd$name %in% c("xi1", "xi2")], gas_levinson(rho)$phi,
               tolerance = 1e-14)
})

test_that("a likelihood-shape term rejects weights of the wrong length", {
  dd <- data.frame(t = 1:12, y = c(rnorm(6), rnorm(6, 3)))
  tm <- term_build(regime(2, time = t), dd)
  psi <- list(level1 = 0.3, gap2 = 2.5, alr1.1 = 1.2, alr2.1 = -0.7)
  seedM <- matrix(0, 12, 5)
  seedM[, 1] <- 1
  run <- function(w, ld = function(e, i) dnorm(dd$y[i], e, log = TRUE),
                  hs = function(e, i) array(-1, c(length(i), 1L, 1L))) {
    term_hessian(tm, eta = rep(0, 12), y = dd$y, logdens = ld,
                 grad = function(e, i) matrix(dd$y[i] - e, ncol = 1L),
                 hess = hs, psi = psi, seed = list(seedM), cols = 2:5,
                 level = 1L, weights = w)
  }
  expect_error(run(1:5), "'weights' must be one number or 12")
  expect_equal(run(2)$gradient, run(rep(2, 12))$gradient)
  md <- data.frame(id = rep(1:3, each = 8), x = rep(1:8, 3))
  md$y <- rnorm(24, 2 * (md$x >= 4.5), 0.4)
  mt <- term_build(jump(x, psi ~ random(~ 1 | id), marginal = TRUE), md)
  S <- matrix(0, 24, 4)
  S[, 1] <- 1
  expect_error(
    term_hessian(mt, rep(0, 24), md$y,
                 function(e, i) dnorm(md$y[i], e, 0.4, log = TRUE),
                 function(e, i) matrix((md$y[i] - e) / 0.16, ncol = 1),
                 function(e, i) array(-1 / 0.16, c(length(i), 1, 1)),
                 list(psi1.mean = 4.5, psi1.sd = 0.5, delta1 = 2), list(S),
                 cols = 2:4, level = 1, weights = rep(3, 6)),
    "'weights' must be one number or 24")
})

test_that("regime() checks what its callbacks return", {
  dd <- data.frame(t = 1:12, y = c(rnorm(6), rnorm(6, 3)))
  tm <- term_build(regime(2, time = t), dd)
  psi <- list(level1 = 0.3, gap2 = 2.5, alr1.1 = 1.2, alr2.1 = -0.7)
  bad <- function(e, i) dnorm(dd$y[1], e[1], log = TRUE)
  expect_error(term_posterior(tm, rep(0, 12), dd$y, logdens = bad, psi = psi),
               "one value per observation")
  seedM <- matrix(0, 12, 5)
  seedM[, 1] <- 1
  hs <- function(h) {
    term_hessian(tm, eta = rep(0, 12), y = dd$y,
                 logdens = function(e, i) dnorm(dd$y[i], e, log = TRUE),
                 grad = function(e, i) matrix(dd$y[i] - e, ncol = 1L),
                 hess = h, psi = psi, seed = list(seedM), cols = 2:5,
                 level = 1L)
  }
  expect_error(hs(function(e, i) -1), "'hess' must return an 12 by 1 by 1")
})

test_that("a break-point development that repeats the intercept is rejected", {
  expect_error(term_build(jump(x, npsi = 2, delta ~ lasso(~ 1)), bd),
               "carries the intercept twice")
  expect_error(term_build(seg(x, gamma ~ lasso(~ 1)), bd),
               "carries the intercept twice")
  expect_s3_class(
    term_build(jump(x, npsi = 2, delta ~ 0 + lasso(~ 1)), bd),
    "modelterms7::SegTerm")
})

test_that("the profile helpers of a break-point term reject non-finite input", {
  b <- term_build(jump(x, psi = 4), bd)
  expect_error(seg_relocate(b, Inf), "finite positions")
  yna <- bd$y
  yna[3] <- NA
  expect_error(seg_start(jump(x), bd, yna), "non-finite")
  expect_error(seg_polish(b, yna), "non-finite")
  expect_error(seg_profile_rss(b, yna), "non-finite")
})

test_that("the prior-location message of a marginal term names no other family", {
  md <- data.frame(id = rep(1:10, each = 8), x = rep(1:8, 10), y = 0)
  err <- tryCatch(
    jump(x, psi ~ random(~ 1 | id,
                         distrib = distributions7::gaussian1_distrib()),
         marginal = TRUE),
    error = conditionMessage)
  expect_match(err, "fixed at zero")
  expect_false(grepl("student_t1", err))
})

test_that("regime() and gas() check their integer arguments", {
  expect_error(suppressWarnings(regime(Inf)), "'k' must be a single integer")
  expect_error(suppressWarnings(regime(2, n_start = Inf)),
               "'n_start' must be a single whole number")
  expect_error(gas(Inf), "'p' must be a single integer")
  expect_error(gas(kappa1 ~ x), "a subformula is given by name after them")
})

test_that("the gas() continuation message is one line", {
  dg <- data.frame(y = rnorm(20))
  tg <- term_build(gas(p = 1, q = 1), dg)
  err <- tryCatch(
    term_continue(tg, list(omega = 0, kappa1 = 0.1, pacf1 = 0.3), rep(0, 20),
                  rep(0, 20), data.frame(y = 1)),
    error = conditionMessage)
  expect_match(err, "needs its `time`")
  expect_false(grepl("\n", err))
})

test_that("a held break-point term draws no position", {
  bh <- seg_hold(term_build(jump(x, psi = 5), bd))
  dr <- term_coef_draw(bh, c(2, -10))
  expect_identical(dr$index, integer(0))
  expect_identical(dr$coef, c(2, -10))
})

test_that("the chart derivatives of gas() apply each partial autocorrelation's link", {
  dd <- data.frame(y = rnorm(30))
  tg <- term_build(gas(p = 1, q = 2,
                       links = list(pacf2 = linkfunctions7::identity_link())),
                   dd)
  links <- term_links(tg)
  zeta <- c(omega = 0.1, kappa1 = 0.2, pacf1 = 0.4, pacf2 = 0.2)
  bz <- function(z) {
    gas_levinson(c(linkfunctions7::linkinv(links$pacf1, z[[3]]), z[[4]]))$phi
  }
  h <- 1e-5
  e <- function(k) replace(numeric(4), k, h)
  ch <- .gas_chart_derivs(zeta, 1L, 2L, links)
  for (j in 1:2) {
    num <- vapply(1:4, function(k) {
      (bz(zeta + e(k))[j] - bz(zeta - e(k))[j]) / (2 * h)
    }, 0)
    # the first link applied to both gave 0.9610 for d b2 / d pacf2, not 1
    expect_equal(as.numeric(ch$d_b[[j]]), num, tolerance = 1e-8)
  }
  v <- c(0.7, -0.4, 1.3, 0.5)
  w <- c(-0.3, 0.9, 0.6, -1.1)
  hb <- function(z) .gas_chart_derivs(z, 1L, 2L, links)$h_b
  t3 <- .gas_chart_derivs3(zeta, 1L, 2L, links, v)$t_b
  tb <- function(z) .gas_chart_derivs3(z, 1L, 2L, links, v)$t_b
  q4 <- .gas_chart_derivs4(zeta, 1L, 2L, links, v, w)$q_b
  for (j in 1:2) {
    expect_equal(t3[[j]], (hb(zeta + h * v)[[j]] - hb(zeta - h * v)[[j]]) /
                   (2 * h), tolerance = 1e-7)
    expect_equal(q4[[j]], (tb(zeta + h * w)[[j]] - tb(zeta - h * w)[[j]]) /
                   (2 * h), tolerance = 1e-7)
  }
})

test_that("term_simulate() of a marginal term checks the length of eta", {
  md <- data.frame(id = rep(1:3, each = 8), x = rep(1:8, 3), y = rnorm(24))
  mt <- term_build(jump(x, psi ~ random(~ 1 | id), marginal = TRUE), md)
  psi <- list(psi1.mean = 4.5, psi1.sd = 0.5, delta1 = 2)
  expect_error(term_simulate(mt, psi, rep(0, 5), draw = FALSE),
               "one value per row of the data the term was built on")
})
