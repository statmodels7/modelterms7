# The pieces of the fifth order of a score-driven filter that do not depend
# on the model: the Levinson-Durbin map's fifth derivative, the chart's
# derivative of any order up to five and the product rule over a set of
# directions. Each is checked against one central difference of the order
# below (Richardson-extrapolated), against the existing orders where they
# overlap, and with a negative control.

test_that("the map's fifth derivative is the derivative of its fourth", {
  set.seed(31)
  for (q in c(5L, 6L)) {
    pc <- stats::runif(q, -0.6, 0.6)
    v <- stats::rnorm(q)
    w <- stats::rnorm(q)
    x <- stats::rnorm(q)
    got <- gas_levinson5(pc, v, w, x)
    d4 <- function(h) gas_levinson4(pc + h * x, v, w)
    dif <- function(h) Map(function(a, b) (a - b) / (2 * h), d4(h), d4(-h))
    rich <- Map(function(a, b) (4 * a - b) / 3, dif(1e-3), dif(2e-3))
    for (i in seq_len(q)) {
      expect_equal(got[[i]], rich[[i]], tolerance = 1e-7,
                   info = sprintf("q = %d, coefficient %d", q, i))
      expect_identical(got[[i]], t(got[[i]]))
    }
    # symmetric in the three directions
    expect_equal(gas_levinson5(pc, x, v, w), got, tolerance = 1e-13)
    expect_true(all(got[[q]] == 0))
  }
  # identically zero below degree five
  expect_true(all(unlist(gas_levinson5(c(0.5, -0.3, 0.2, 0.1),
                                       rep(1, 4), rep(1, 4), rep(1, 4))) == 0))
  # the negative control: a perturbed reference is told apart
  pc <- c(0.3, -0.2, 0.25, 0.1, -0.15)
  g <- gas_levinson5(pc, rep(1, 5), rep(0.5, 5), rep(-1, 5))
  expect_gt(max(abs(unlist(g))), 1e-3)
})

chart_links <- function(nm) {
  stats::setNames(lapply(nm, function(j) {
    if (startsWith(j, "pacf")) linkfunctions7::rhobit_link()
    else if (startsWith(j, "kappa")) linkfunctions7::log_link()
    else linkfunctions7::identity_link()
  }), nm)
}

test_that("the chart of any order reproduces the orders two to four", {
  set.seed(32)
  for (q in c(1L, 3L, 4L)) {
    p <- 2L
    nm <- .gas_base_params(p, q)
    links <- chart_links(nm)
    zeta <- stats::setNames(stats::rnorm(length(nm), 0, 0.4), nm)
    vz <- stats::rnorm(length(nm))
    wz <- stats::rnorm(length(nm))
    c2 <- .gas_chart_derivs(zeta, p, q, links)
    c3 <- .gas_chart_derivs3(zeta, p, q, links, vz)
    c4 <- .gas_chart_derivs4(zeta, p, q, links, vz, wz)
    g2 <- .gas_chart_contract(zeta, p, q, links, list())
    g3 <- .gas_chart_contract(zeta, p, q, links, list(vz))
    g4 <- .gas_chart_contract(zeta, p, q, links, list(vz, wz))
    lab <- sprintf("q = %d", q)
    expect_equal(g2$omega, c2$h_omega, tolerance = 1e-13, info = lab)
    expect_equal(g2$b, c2$h_b, tolerance = 1e-13, info = lab)
    expect_equal(g3$omega, c3$t_omega, tolerance = 1e-13, info = lab)
    expect_equal(g3$a, c3$t_a, tolerance = 1e-13, info = lab)
    expect_equal(g3$b, c3$t_b, tolerance = 1e-13, info = lab)
    expect_equal(g4$a, c4$q_a, tolerance = 1e-13, info = lab)
    expect_equal(g4$b, c4$q_b, tolerance = 1e-13, info = lab)
  }
})

test_that("the chart's fifth derivative is the derivative of its fourth", {
  set.seed(33)
  for (q in c(2L, 5L)) {
    p <- 2L
    nm <- .gas_base_params(p, q)
    links <- chart_links(nm)
    zeta <- stats::setNames(stats::rnorm(length(nm), 0, 0.4), nm)
    vz <- stats::rnorm(length(nm))
    wz <- stats::rnorm(length(nm))
    xz <- stats::rnorm(length(nm))
    got <- .gas_chart_contract(zeta, p, q, links, list(vz, wz, xz))
    c4 <- function(h) .gas_chart_derivs4(zeta + h * xz, p, q, links, vz, wz)
    dif <- function(h, f) (f(c4(h)) - f(c4(-h))) / (2 * h)
    rich <- function(f) (4 * dif(1e-3, f) - dif(2e-3, f)) / 3
    lab <- sprintf("q = %d", q)
    expect_equal(got$omega, rich(function(r) r$q_omega), tolerance = 1e-7,
                 info = lab)
    for (i in seq_len(p)) {
      expect_equal(got$a[[i]], rich(function(r) r$q_a[[i]]), tolerance = 1e-7,
                   info = lab)
    }
    for (j in seq_len(q)) {
      ref <- rich(function(r) r$q_b[[j]])
      expect_equal(got$b[[j]], ref, tolerance = 1e-6, info = lab)
      # the negative control: one entry moved by one per cent is detected
      bad <- got$b[[j]]
      k <- which.max(abs(bad))
      bad[k] <- bad[k] * 1.01
      if (abs(ref[k]) > 1e-3) {
        expect_false(isTRUE(all.equal(bad, ref, tolerance = 1e-6)))
      }
    }
    expect_equal(.gas_chart_contract(zeta, p, q, links, list(xz, vz, wz)),
                 got, tolerance = 1e-12)
  }
})

test_that("the product rule over directions reproduces the fourth order's", {
  set.seed(34)
  mk <- 4L
  vks <- list(stats::rnorm(mk), stats::rnorm(mk))
  sym <- function() { z <- matrix(stats::rnorm(mk * mk), mk); (z + t(z)) / 2 }
  m <- replicate(8, sym(), simplify = FALSE)
  v0 <- stats::rnorm(2)
  d1 <- list(stats::rnorm(mk), stats::rnorm(mk))
  A4 <- .gas_jet(v0[1], d1[[1]], m[[1]], m[[2]], m[[3]], m[[4]], vks)
  B4 <- .gas_jet(v0[2], d1[[2]], m[[5]], m[[6]], m[[7]], m[[8]], vks)
  An <- .gas_jetn(v0[1], d1[[1]], list(s = m[[1]], s1 = m[[2]], s2 = m[[3]],
                                       s12 = m[[4]]), vks)
  Bn <- .gas_jetn(v0[2], d1[[2]], list(s = m[[5]], s1 = m[[6]], s2 = m[[7]],
                                       s12 = m[[8]]), vks)
  expect_equal(.gas_prodn(An, Bn), .gas_prod4(A4, B4), tolerance = 1e-13)
})

test_that("the product rule of the fifth order is exact on exponentials", {
  # A = exp(a.u) and B = exp(b.u): every derivative along S with two free
  # slots is the value times prod(a.s) a a', and the product is exp((a+b).u)
  set.seed(35)
  mk <- 3L
  a <- stats::rnorm(mk, 0, 0.5)
  b <- stats::rnorm(mk, 0, 0.5)
  dirs <- list(stats::rnorm(mk), stats::rnorm(mk), stats::rnorm(mk))
  jet_exp <- function(c0) {
    val <- 1.3
    mats <- list()
    for (S in .gas_subsets(3L)) {
      f <- val
      for (s in S) f <- f * sum(c0 * dirs[[s]])
      mats[[paste0("s", paste(S, collapse = ""))]] <- f * outer(c0, c0)
    }
    .gas_jetn(val, val * c0, mats, dirs)
  }
  got <- .gas_prodn(jet_exp(a), jet_exp(b))
  ab <- a + b
  ref <- 1.3^2 * prod(vapply(dirs, function(d) sum(ab * d), numeric(1))) *
    outer(ab, ab)
  expect_equal(got, ref, tolerance = 1e-13)
  expect_equal(.gas_prodn(jet_exp(b), jet_exp(a)), got, tolerance = 1e-13)
  # the negative control: dropping one direction's term changes the result
  A <- jet_exp(a)
  A$mats$s123 <- A$mats$s123 * 1.01
  expect_false(isTRUE(all.equal(.gas_prodn(A, jet_exp(b)), ref,
                                tolerance = 1e-8)))
})
