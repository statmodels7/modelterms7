# The compiled route of the joint filter against the R route. The R route
# is fed by callbacks that compose the family's link-scale derivatives from
# the distributions7 generics, as a fitting layer does; the compiled route
# reads the same quantities through the joint registry. The two are compared
# on the filtered predictors, the Jacobian, the curvature and the reverse
# pass, for the covered families, the three scalings, cross loadings and
# groups.

# the driving vector u = I_GG^-d g_G and its derivative C in the filtered
# predictors, at the predictors `e` of the filtered parameters, the others
# held at their static values
jc_callbacks <- function(d, y, Eall, G, scaling) {
  params <- d@params
  P <- length(params)
  key <- function(a, b) {
    nm <- distributions7::hess_names(params)
    k <- paste(params[a], params[b], sep = "_")
    if (k %in% nm) k else paste(params[b], params[a], sep = "_")
  }
  at <- function(e, i) {
    full <- Eall[i, ]
    full[G] <- e
    th <- stats::setNames(lapply(seq_len(P), function(a) {
      lk <- d@link_params[[params[a]]]
      linkfunctions7::link_bounds_clamp(linkfunctions7::linkinv(lk, full[a]),
                                        lk@link_bounds)
    }), params)
    g <- distributions7::distrib_gradient(d, y[i], th, scale = "link")
    h <- distributions7::distrib_hessian(d, y[i], th, scale = "link")
    gG <- vapply(G, function(a) g[[params[a]]], 0)
    H <- outer(G, G, Vectorize(function(a, b) h[[key(a, b)]]))
    if (scaling == 0) return(list(u = gG, C = H))
    E <- distributions7::distrib_expected_hessian(d, y[i], th, scale = "link")
    dE <- distributions7::distrib_dexpected_hessian(d, y[i], th,
                                                    scale = "link")
    I <- -outer(G, G, Vectorize(function(a, b) E[[key(a, b)]]))
    dI <- lapply(G, function(c) -outer(G, G, Vectorize(function(a, b)
      dE[[distributions7::dexpected_key(params, a, b, c)]])))
    if (scaling == 1) {
      S <- solve(I)
      dS <- lapply(dI, function(x) -S %*% x %*% S)
    } else {
      ev <- eigen(I, symmetric = TRUE)
      Q <- ev$vectors
      r <- sqrt(ev$values)
      S <- Q %*% (t(Q) / r)
      dS <- lapply(dI, function(x)
        -S %*% (Q %*% ((t(Q) %*% x %*% Q) / outer(r, r, `+`)) %*% t(Q)) %*% S)
    }
    C <- S %*% H
    for (k in seq_along(G)) C[, k] <- C[, k] + as.numeric(dS[[k]] %*% gG)
    list(u = as.numeric(S %*% gG), C = C)
  }
  list(score = function(e, i) at(e, i)$u,
       curvature = function(e, i) at(e, i)$C)
}

jc_fast <- function(d, y, Eall, G, scaling) {
  params <- d@params
  lks <- d@link_params[params]
  route <- distributions7::distrib_scalar_route(d)
  list(family = route$name,
       constants = lapply(route$constants, as.numeric),
       links = vapply(lks, function(lk)
         linkfunctions7::link_scalar_route(lk)$name, ""),
       link_par = lapply(lks, function(lk)
         linkfunctions7::link_scalar_route(lk)$par),
       bounds = lapply(lks, function(lk) as.numeric(lk@link_bounds)),
       y = as.numeric(y), eta = Eall, filtered = as.integer(G),
       scaling = scaling)
}

jc_case <- function(fam, n, seed) {
  set.seed(seed)
  if (fam == "gaussian1") {
    d <- distributions7::gaussian1_distrib()
    Eall <- cbind(1 + 0.2 * sin(seq_len(n)), log(1.3) + 0.1 * cos(seq_len(n)))
    y <- stats::rnorm(n, 1, 1.3)
    G <- 1:2
  } else {
    d <- distributions7::student_t1_distrib()
    Eall <- cbind(0.1 + 0.05 * sin(seq_len(n)), 0.1 * cos(seq_len(n)),
                  log(6) + 0.2 * sin(seq_len(n) / 3))
    y <- 0.1 + stats::rt(n, df = 6)
    G <- 2:3
  }
  list(d = d, Eall = Eall, y = y, G = G)
}

# the terms of the parameters `pnames`; cross NULL (no cross loadings),
# "sym" (every level loads every other score at the lags 1..p) or "upper"
# (the first level loads the second score at the lags 1..p + 2, beyond its
# own p, and the others load nothing)
jc_terms <- function(K, p, q, cross, by, data, scaling = 0,
                     pnames = c("a", "b", "c")[seq_len(K)]) {
  mk <- function(k) {
    pc <- if (identical(cross, "sym")) {
      stats::setNames(rep(p, K - 1L), setdiff(pnames, k))
    } else if (identical(cross, "upper") && k == pnames[1L]) {
      stats::setNames(p + 2, pnames[2L])
    }
    if (by) gas(p = p, q = q, by = g, time = t, p_cross = pc,
                scaling = scaling)
    else gas(p = p, q = q, time = t, p_cross = pc, scaling = scaling)
  }
  stats::setNames(lapply(pnames, function(k) term_build(mk(k), data)), pnames)
}

jc_psi <- function(jt, seed) {
  set.seed(seed)
  lapply(stats::setNames(jt@params, jt@params), function(k) {
    nm <- term_params(jt@terms[[k]])
    v <- stats::setNames(numeric(length(nm)), nm)
    for (j in nm) {
      v[[j]] <- if (j == "omega") stats::runif(1, -0.05, 0.05) else
        if (grepl("_", j)) stats::runif(1, -0.03, 0.03) else
        if (startsWith(j, "kappa")) stats::runif(1, 0.02, 0.06) else
        stats::runif(1, 0.3, 0.7)
    }
    as.list(v)
  })
}

test_that("the compiled joint filter is the R route", {
  n <- 60
  dd <- data.frame(t = seq_len(n), g = rep(c("a", "b"), each = n / 2))
  for (fam in c("gaussian1", "student_t1")) {
    cs <- jc_case(fam, n, 2)
    for (scaling in c(0, 0.5, 1)) {
      for (cfg in list(list(1, 1, NULL, FALSE), list(2, 1, "sym", TRUE),
                       list(1, 2, "sym", FALSE), list(1, 1, "upper", TRUE))) {
        tms <- jc_terms(2, cfg[[1]], cfg[[2]], cfg[[3]], cfg[[4]], dd,
                        scaling = scaling, pnames = cs$d@params[cs$G])
        jt <- gas_joint(tms)
        psi <- jc_psi(jt, 9)
        eta <- cs$Eall[, cs$G]
        cb <- jc_callbacks(cs$d, cs$y, cs$Eall, cs$G, scaling)
        fast <- jc_fast(cs$d, cs$y, cs$Eall, cs$G, scaling)
        lab <- sprintf("%s, d = %s, p = %d, q = %d, cross = %s, by = %s", fam,
                       format(scaling), cfg[[1]], cfg[[2]], format(cfg[[3]]),
                       cfg[[4]])
        r <- term_filter(jt, eta, cs$y, cb$score, cb$curvature, psi)
        cc <- term_filter(jt, eta, cs$y, cb$score, cb$curvature, psi,
                          fast = fast)
        # the two routes invert the information differently (solve() against
        # an eigen decomposition), and the last bits that the platform's
        # LAPACK leaves accumulate along the recursion: 6e-13 on Linux with
        # R 4.5, so the predictors are held to the Jacobian's tolerance
        expect_equal(cc$eta, r$eta, tolerance = 1e-11, info = lab)
        for (k in jt@params) {
          expect_equal(cc$jacobian[[k]], r$jacobian[[k]], tolerance = 1e-11,
                       info = paste(lab, k))
        }
        expect_equal(cc$curv, r$curv, tolerance = 1e-11, info = lab)
        set.seed(3)
        gw <- matrix(stats::rnorm(2 * n), n, 2)
        ar <- term_adjoint(jt, eta, cs$y, cb$score, cb$curvature, psi, gw)
        ac <- term_adjoint(jt, eta, cs$y, cb$score, cb$curvature, psi, gw,
                           stored_curv = cc$curv)
        expect_equal(ac$deta, ar$deta, tolerance = 1e-11, info = lab)
        expect_equal(ac$dscore, ar$dscore, tolerance = 1e-11, info = lab)
      }
    }
  }
})

# every family that the joint registry of distributions7 covers, with the
# orders it covers: 4 (all entries) or 2 (the gradient and the Hessian, for an
# unscaled filter only)
jc_families <- list(
  gaussian1 = 4, gamma1 = 4, negbin2 = 4, beta1 = 4, cauchy = 4,
  logistic = 4, gaussian2 = 4, gaussian3 = 4, lognormal1 = 4, invgauss1 = 4,
  invgauss2 = 4, gamma2 = 4, gpd = 4, vonmises1 = 4, vonmises2 = 4,
  weibull3 = 4, lognormal2 = 4, student_t1 = 4, student_t2 = 4,
  gengamma1 = 4, gengamma2 = 4, gumbel = 4, laplace = 4, laplace2 = 4,
  weibull1 = 4, beta2 = 4, enet = 4, negbin1 = 4, pig1 = 4, pig2 = 4,
  betabinom1 = 4, betabinom2 = 4, skewnormal1 = 4, skewt = 4,
  pseudohuber = 4, pseudohuber2 = 4, skewnormal2 = 4)

jc_family <- function(nm) {
  ctor <- get(paste0(nm, "_distrib"), envir = asNamespace("distributions7"))
  if (startsWith(nm, "betabinom")) ctor(size = 9) else ctor()
}

# one covered family, its first two parameters filtered, with cross loadings
# and groups, at the scalings its coverage `cov` admits; the static
# predictors vary around `base` (0.3 for every parameter by default)
jc_check_family <- function(nm, d, cov, base = NULL) {
  n <- 40
  dd <- data.frame(t = seq_len(n), g = rep(c("a", "b"), each = n / 2))
  P <- length(d@params)
  G <- 1:2
  if (is.null(base)) base <- rep(0.3, P)
  Eall <- vapply(seq_len(P), function(a) base[a] + 0.1 * sin(seq_len(n) + a),
                 numeric(n))
  th <- lapply(seq_len(P), function(a)
    linkfunctions7::linkinv(d@link_params[[d@params[a]]], Eall[, a]))
  names(th) <- d@params
  set.seed(5)
  y <- vapply(seq_len(n), function(i)
    distributions7::distrib_rng(d, 1, lapply(th, `[`, i)), 0)
  for (scaling in if (cov == 4) c(0, 0.5, 1) else c(0, 0.5)) {
    tms <- jc_terms(2, 1, 1, "sym", TRUE, dd, scaling = scaling,
                    pnames = d@params[G])
    jt <- gas_joint(tms)
    psi <- jc_psi(jt, 9)
    eta <- Eall[, G]
    fast <- jc_fast(d, y, Eall, G, scaling)
    lab <- sprintf("%s, d = %s", nm, format(scaling))
    cf <- .gas_joint_coefs(jt, psi)
    direct <- gas_joint_filter_cpp(eta, cf, length(term_params(jt)),
                                   jt@build_info$order, fast)
    if (cov == 2 && scaling != 0) {
      # the expected information is not in the registry for this family
      expect_null(direct, info = lab)
      next
    }
    expect_false(is.null(direct), info = lab)
    cb <- jc_callbacks(d, y, Eall, G, scaling)
    r <- term_filter(jt, eta, y, cb$score, cb$curvature, psi)
    cc <- term_filter(jt, eta, y, cb$score, cb$curvature, psi, fast = fast)
    expect_equal(cc$eta, r$eta, tolerance = 1e-11, info = lab)
    for (k in jt@params) {
      expect_equal(cc$jacobian[[k]], r$jacobian[[k]], tolerance = 1e-10,
                   info = paste(lab, k))
    }
    expect_equal(cc$curv, r$curv, tolerance = 1e-10, info = lab)
  }
}

test_that("the compiled joint filter is the R route for every covered family", {
  for (nm in names(jc_families)) {
    jc_check_family(nm, jc_family(nm), jc_families[[nm]])
  }
})

# wrapped families, whose joint entries distributions7 composes from the inner
# family's, with their coverage
jc_wrapped <- list(
  fixed_t1 = list(function() distributions7::fixed(
    distributions7::student_t1_distrib(), nu = 5), 4),
  fixed_skewt = list(function() distributions7::fixed(
    distributions7::skewt_distrib(), nu = 6), 4),
  zi_poisson = list(function() distributions7::zero_inflated(
    distributions7::poisson_distrib()), 4),
  zi_negbin2 = list(function() distributions7::zero_inflated(
    distributions7::negbin2_distrib()), 4),
  fixed_zi_negbin2 = list(function() distributions7::fixed(
    distributions7::zero_inflated(distributions7::negbin2_distrib()),
    theta = 2), 4),
  za_negbin2 = list(function() distributions7::zero_adjusted(
    distributions7::negbin2_distrib()), 4),
  za_gamma1 = list(function() distributions7::zero_adjusted(
    distributions7::gamma1_distrib()), 4),
  tf_gaussian1 = list(function() distributions7::transformation(
    distributions7::gaussian1_distrib(), distributions7::exp_transform()), 4),
  # a folded family has no information in its center at a center of zero,
  # so the center is held two scales from the fold
  fo_gaussian1 = list(function() distributions7::folded(
    distributions7::gaussian1_distrib()), 4, c(2, 0)),
  fo_skewnormal1 = list(function() distributions7::folded(
    distributions7::skewnormal1_distrib()), 4, c(2, 0, 0.3)),
  tr_negbin2 = list(function() distributions7::truncated(
    distributions7::negbin2_distrib(), lower = 1), 4),
  tr_gaussian1 = list(function() distributions7::truncated(
    distributions7::gaussian1_distrib(), lower = -1, upper = 3), 4)
)

test_that("the compiled joint filter is the R route for wrapped families", {
  for (nm in names(jc_wrapped)) {
    w <- jc_wrapped[[nm]]
    jc_check_family(nm, w[[1]](), w[[2]], if (length(w) > 2L) w[[3]])
  }
})

test_that("the compiled joint filter gives the same result on any thread count", {
  n <- 120
  dd <- data.frame(t = rep(seq_len(12), 10), g = rep(letters[1:10], each = 12))
  cs <- jc_case("student_t1", n, 5)
  tms <- jc_terms(2, 1, 1, NULL, TRUE, dd, scaling = 0.5)
  names(tms) <- c("sigma", "nu")
  jt <- gas_joint(tms)
  psi <- jc_psi(jt, 2)
  fast <- jc_fast(cs$d, cs$y, cs$Eall, cs$G, 0.5)
  one <- term_filter(jt, cs$Eall[, cs$G], cs$y, NULL, NULL, psi, fast = fast)
  four <- term_filter(jt, cs$Eall[, cs$G], cs$y, NULL, NULL, psi, fast = fast,
                      threads = 4L)
  expect_identical(four, one)
})

test_that("a family that the registry does not know leaves the context inert", {
  n <- 40
  dd <- data.frame(t = seq_len(n))
  d <- distributions7::gamma1_distrib()
  Eall <- cbind(rep(log(2), n), rep(log(0.5), n))
  set.seed(1)
  y <- stats::rgamma(n, 2, 1)
  tms <- jc_terms(2, 1, 1, NULL, FALSE, dd)
  names(tms) <- d@params
  jt <- gas_joint(tms)
  psi <- jc_psi(jt, 4)
  cb <- jc_callbacks(d, y, Eall, 1:2, 0)
  fast <- jc_fast(d, y, Eall, 1:2, 0)
  fast$family <- "NoSuchDistrib"
  expect_null(gas_joint_filter_cpp(Eall, .gas_joint_coefs(jt, psi),
                                   length(term_params(jt)), jt@build_info$order,
                                   fast))
  expect_identical(
    term_filter(jt, Eall, y, cb$score, cb$curvature, psi, fast = fast),
    term_filter(jt, Eall, y, cb$score, cb$curvature, psi))
})

test_that("a curvature of the wrong shape is rejected", {
  n <- 40
  dd <- data.frame(t = seq_len(n))
  tms <- jc_terms(2, 1, 1, NULL, FALSE, dd)
  names(tms) <- c("mu", "sigma")
  jt <- gas_joint(tms)
  psi <- jc_psi(jt, 4)
  expect_error(term_adjoint(jt, matrix(0, n, 2), NULL, NULL, NULL, psi,
                            matrix(0, n, 2), stored_curv = array(0, c(n, 2, 1))),
               "dimension")
})
