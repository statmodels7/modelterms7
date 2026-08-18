# Nonlinear parametric terms: the Jacobian as the design block, the two
# ways of supplying the function, links, and parameter submodels.

set.seed(13)
n <- 60
dd <- data.frame(x = seq(0, 3, length.out = n),
                 g = factor(rep(c("a", "b"), length.out = n)))
dd$y <- 2 * exp(-1.3 * dd$x) + rnorm(n, sd = 0.05)

test_that("the block is the jacobian, and it is exact", {
  spec <- nl(~ a * exp(-r * x), start = list(a = 2, r = 1.3))
  built <- term_build(spec, dd)
  expect_identical(term_coef_names(built), c("nl.a", "nl.r"))
  expect_identical(built@deriv_mode, "symbolic")

  b <- c(2, 1.3)
  # the contribution and its derivative, from the definition
  f <- function(v) v[1] * exp(-v[2] * dd$x)
  expect_equal(term_value(built), f(b), tolerance = 1e-12)
  expect_equal(unname(term_matrix(built)), numDeriv::jacobian(f, b),
               tolerance = 1e-7)
})

test_that("refreshing moves the block to the new parameters", {
  built <- term_build(nl(~ a * exp(-r * x), start = list(a = 2, r = 1.3)), dd)
  b2 <- c(0.5, 2.5)
  moved <- term_refresh(built, b2)

  f <- function(v) v[1] * exp(-v[2] * dd$x)
  expect_equal(term_value(moved), f(b2), tolerance = 1e-12)
  expect_equal(unname(term_matrix(moved)), numDeriv::jacobian(f, b2),
               tolerance = 1e-7)
  # the term it came from is untouched
  expect_equal(term_value(built), f(c(2, 1.3)), tolerance = 1e-12)
  expect_error(term_refresh(built, 1), "length 2")

  # and an ordinary term ignores the refresh, its block being data alone
  lin <- term_build(linpar(~x), dd)
  expect_identical(term_matrix(term_refresh(lin, c(0, 0))),
                   term_matrix(lin))
})

test_that("a Gauss-Newton step on the linearization finds the truth", {
  # the whole point of the Jacobian contract: iterate value + J, and the
  # nonlinear least squares problem is solved by linear steps
  built <- term_build(nl(~ a * exp(-r * x), start = list(a = 1, r = 0.5)), dd)
  b <- c(linkfunctions7::linkfun(linkfunctions7::identity_link(), 1), 0.5)
  for (it in 1:25) {
    cur <- term_refresh(built, b)
    J <- term_matrix(cur)
    resid <- dd$y - term_value(cur)
    b <- b + solve(crossprod(J) + 1e-10 * diag(2), crossprod(J, resid))
  }
  expect_equal(as.numeric(b), c(2, 1.3), tolerance = 0.05)
})

test_that("a link holds a parameter in its own set", {
  spec <- nl(~ a * exp(-r * x),
             links = list(r = linkfunctions7::log_link()),
             start = list(a = 2, r = 1.3))
  built <- term_build(spec, dd)
  # the coefficient is now log(r), so the block carries the chain rule
  b <- c(2, log(1.3))
  f <- function(v) v[1] * exp(-exp(v[2]) * dd$x)
  expect_equal(term_value(built), f(b), tolerance = 1e-10)
  expect_equal(unname(term_matrix(built)), numDeriv::jacobian(f, b),
               tolerance = 1e-6)
  # and no coefficient can produce a negative rate
  expect_gt(min(exp(seq(-50, 50, length.out = 11))), 0)
})

test_that("a parameter can be developed with covariates", {
  spec <- nl(~ a * exp(-r * x), subformulas = list(a = ~g),
             start = list(r = 1.3))
  built <- term_build(spec, dd)
  expect_identical(term_coef_names(built),
                   c("nl.a.(Intercept)", "nl.a.gb", "nl.r"))

  Z <- stats::model.matrix(~g, dd)
  b <- c(2, -0.4, 1.3)
  f <- function(v) as.numeric(Z %*% v[1:2]) * exp(-v[3] * dd$x)
  moved <- term_refresh(built, b)
  expect_equal(term_value(moved), f(b), tolerance = 1e-10)
  expect_equal(unname(term_matrix(moved)), numDeriv::jacobian(f, b),
               tolerance = 1e-6)

  res <- check_term(spec, dd, verbose = FALSE)
  expect_true(all(res$status == "OK"),
              info = paste(res$check[res$status != "OK"], collapse = ", "))
})

test_that("a function is accepted and differenced, and says so", {
  fx <- function(x, theta) theta$a * exp(-theta$r * x)
  spec <- nl(fx, params = c("a", "r"), x = x, start = list(a = 2, r = 1.3))
  built <- term_build(spec, dd)
  expect_identical(built@deriv_mode, "numeric")

  b <- c(2, 1.3)
  f <- function(v) v[1] * exp(-v[2] * dd$x)
  expect_equal(term_value(built), f(b), tolerance = 1e-12)
  # differenced, so held to the accuracy of a central difference
  expect_equal(unname(term_matrix(built)), numDeriv::jacobian(f, b),
               tolerance = 1e-6)

  res <- check_term(spec, dd, verbose = FALSE)
  expect_true(all(res$status == "OK"))
})

test_that("an expression deriv() cannot read falls back to differences", {
  # besselJ has no symbolic derivative in R's table, so the formula route
  # must notice and difference instead of failing
  spec <- nl(~ a * besselJ(x + 1, 0) + r * x, start = list(a = 1, r = 1))
  built <- term_build(spec, dd)
  expect_identical(built@deriv_mode, "formula_fd")
  b <- c(1.5, 0.4)
  f <- function(v) v[1] * besselJ(dd$x + 1, 0) + v[2] * dd$x
  moved <- term_refresh(built, b)
  expect_equal(term_value(moved), f(b), tolerance = 1e-12)
  expect_equal(unname(term_matrix(moved)), numDeriv::jacobian(f, b),
               tolerance = 1e-6)
})

test_that("what only a formula can do is refused to a function", {
  fx <- function(x, theta) theta$a * x
  expect_error(nl(fx), "'params' must name the parameters")
  expect_error(nl(fx, params = "a", subformulas = list(a = ~g)),
               "only a formula says")
  expect_error(nl(y ~ a * x), "one-sided")
  expect_error(nl("not a formula"), "formula or a function")
  expect_error(nl(~ a * x, links = list(nope = 1)), NA)
  expect_error(term_build(nl(~ a * x, links = list(nope = 1)), dd),
               "not a parameter")
  expect_error(term_build(nl(~ x + 1), dd), "no parameters to estimate")
})

test_that("the interpreter routes it and print says how it differentiates", {
  out <- interpret_formula(y ~ nl(~ a * exp(-r * x)), dd)
  expect_named(out$terms, c("linpar", "nl(~a * exp(-r * x))"))
  built <- term_build(out$terms[["nl(~a * exp(-r * x))"]], dd)
  expect_output(print(built), "symbolic derivatives")
  expect_output(print(nl(~ a * x)), "specification")
})

test_that("a penalty is asked for inside the subformula", {
  # the sub-term that carries it declares it, with its own hyperparameter;
  # an unpenalized term reports none
  expect_length(term_penalties(term_build(nl(~ a * exp(-r * x)), dd)), 0L)
  built <- term_build(nl(~ a * exp(-r * x), a ~ lasso(~g),
                         start = list(r = 1.3)), dd)
  ent <- term_penalties(built)
  expect_length(ent, 1L)
  expect_identical(ent[[1L]]$name, "a::lasso(~g)")
  expect_true(penalties7::has_prox(ent[[1L]]$penalty))
  expect_false(term_smooth(built))
  # the whole of the block is not penalized: the rate is left free
  expect_null(term_penalty(built))

  # the removed arguments are reported by name
  expect_error(nl(~ a * x, penalty = penalties7::lasso_penalty),
               "unused argument 'penalty'")
  expect_error(nl(~ a * x, penalize = "a"), "unused argument 'penalize'")
})

test_that("a two-sided formula in ... is the subformula it names", {
  a1 <- term_build(nl(~ a * exp(-r * x), a ~ g, start = list(r = 1.3)), dd)
  a2 <- term_build(nl(~ a * exp(-r * x), subformulas = list(a = ~g),
                      start = list(r = 1.3)), dd)
  expect_identical(term_coef_names(a1), term_coef_names(a2))
  expect_equal(term_matrix(a1), term_matrix(a2))

  # a parameter carries one subformula, whichever spelling supplies it
  expect_error(nl(~ a * exp(-r * x), a ~ g, subformulas = list(a = ~g)),
               "two subformulas")
  expect_error(nl(~ a * exp(-r * x), a ~ g, a ~ x), "two subformulas")
  # what lands in ... must be a two-sided formula naming a parameter
  expect_error(nl(~ a * exp(-r * x), ~g), "two-sided")
  expect_error(nl(~ a * exp(-r * x), "g"), "two-sided")
  expect_error(nl(~ a * exp(-r * x), (a + 1) ~ g), "two-sided")
})

test_that("a subformula of ~1 reproduces the plain parameter", {
  plain <- term_build(nl(~ a * exp(-r * x), start = list(a = 2, r = 1.3)), dd)
  one <- term_build(nl(~ a * exp(-r * x), a ~ 1,
                       start = list(a = 2, r = 1.3)), dd)
  expect_identical(term_coef_names(one), c("nl.a.(Intercept)", "nl.r"))
  expect_equal(unname(as.matrix(term_matrix(one))),
               unname(as.matrix(term_matrix(plain))), tolerance = 1e-12)
  expect_equal(term_value(one), term_value(plain), tolerance = 1e-12)
})

test_that("a subformula takes a term, whose penalty the term reports", {
  spec <- nl(~ a * exp(-r * x), a ~ ridge(~g), start = list(a = 2, r = 1.3))
  built <- term_build(spec, dd)

  # the interpreter's convention: the parametric block first (here the
  # intercept, which is the population value), then the penalized departures
  ent <- term_penalties(built)
  expect_length(ent, 1L)
  expect_identical(ent[[1L]]$name, "a::ridge(~g)")
  rsub <- term_build(ridge(~g), dd)
  expect_identical(term_coef_names(built)[ent[[1L]]$index],
                   paste("nl.a", term_coef_names(rsub), sep = "."))
  expect_true(term_smooth(built))

  # the jacobian carries the chain through the whole sub-design
  Z <- cbind(1, as.matrix(term_matrix(rsub)))
  k <- ncol(Z)
  f <- function(v) as.numeric(Z %*% v[seq_len(k)]) * exp(-v[k + 1L] * dd$x)
  b <- c(2, 0.3, -0.1, 1.3)[seq_len(k + 1L)]
  moved <- term_refresh(built, b)
  expect_equal(term_value(moved), f(b), tolerance = 1e-10)
  expect_equal(unname(as.matrix(term_matrix(moved))),
               numDeriv::jacobian(f, b), tolerance = 1e-6)

  res <- check_term(spec, dd, verbose = FALSE)
  expect_true(all(res$status == "OK"),
              info = paste(res$check[res$status != "OK"], collapse = ", "))
})

test_that("a smooth in a subformula predicts through its own blueprint", {
  d2 <- dd
  d2$z <- runif(n)
  spec <- nl(~ a * exp(-r * x), a ~ s(z, k = 6), start = list(a = 2, r = 1.3))
  built <- term_build(spec, d2)
  ent <- term_penalties(built)
  expect_length(ent, 1L)
  expect_identical(ent[[1L]]$name, "a::s(z, k = 6)")
  expect_true(term_smooth(built))

  # reapplied on a subset of the rows, the contribution is the fitted one
  cf <- c(1.5, stats::rnorm(length(term_coef_names(built)) - 2L, sd = 0.1), 1.1)
  keep <- seq_len(20)
  expect_equal(term_value(built, coef = cf,
                          newdata = d2[keep, , drop = FALSE]),
               term_value(built, coef = cf)[keep], tolerance = 1e-10)
})

test_that("a sparse sub-design stays sparse through the jacobian", {
  d3 <- data.frame(x = rep(seq(0, 3, length.out = 12), 15),
                   g3 = factor(rep(sprintf("g%02d", 1:15), each = 12)))
  spec <- nl(~ a * exp(-r * x), a ~ random(~1 | g3),
             start = list(a = 2, r = 1.3))
  built <- term_build(spec, d3)
  expect_s4_class(term_matrix(built), "Matrix")
  # and the values agree with the dense arithmetic
  ent <- term_penalties(built)
  expect_length(ent, 1L)
  rsub <- term_build(random(~1 | g3), d3)
  Z <- cbind(1, as.matrix(term_matrix(rsub)))
  k <- ncol(Z)
  b <- c(2, stats::rnorm(k - 1L, sd = 0.1), 1.3)
  f <- function(v) as.numeric(Z %*% v[seq_len(k)]) * exp(-v[k + 1L] * d3$x)
  moved <- term_refresh(built, b)
  expect_equal(unname(as.matrix(term_matrix(moved))),
               numDeriv::jacobian(f, b), tolerance = 1e-6)
})

test_that("a submodel must be a fixed design", {
  expect_error(term_build(nl(~ a * exp(-r * x), a ~ gas(p = 1, q = 1)), dd),
               "structural")
  expect_error(term_build(nl(~ a * exp(-r * x), a ~ seg(x)), dd),
               "moves with its coefficients")
})

test_that("term_value answers on other rows, where the block cannot", {
  spec <- nl(~ a * exp(-r * x), start = list(a = 2, r = 1.3))
  built <- term_build(spec, dd)
  cf <- c(2.1, 1.2)
  nd <- data.frame(x = c(0, 0.5, 1.5, 3), g = factor(c("a", "b", "a", "b")))
  got <- term_value(built, coef = cf, newdata = nd)
  # the link is the identity here, so the contribution is the function
  expect_equal(got, cf[1L] * exp(-cf[2L] * nd$x), tolerance = 1e-12)
  # on the fitting rows it is the ordinary answer
  expect_equal(term_value(built, coef = cf, newdata = dd),
               term_value(built, coef = cf), tolerance = 1e-12)
  # and a parameter with a submodel reapplies its levels rather than
  # rebuilding them: a subset that drops one still gets the same value
  bs <- term_build(nl(~ a * exp(-r * x), subformulas = list(a = ~g),
                      start = list(r = 1.3)), dd)
  cs <- c(2, 0.3, 1.1)
  keep <- which(dd$g == "a")
  expect_equal(term_value(bs, coef = cs,
                          newdata = droplevels(dd[keep, , drop = FALSE])),
               term_value(bs, coef = cs)[keep], tolerance = 1e-12)
})


test_that("the block's derivative along a direction is the contraction's adjoint", {
  # term_block_deriv() answers per ENTRY of the block where
  # term_block_contract() answers per coefficient, and neither computes the
  # other: the gradient of a marginal criterion needs the contraction, its
  # Hessian needs the direction. The identity that ties them shares no code
  # with either and must hold exactly:
  #     sum_ij A_ij (dX/dbeta . v)_ij  ==  v' term_block_contract(A)
  set.seed(31)
  n <- 200L
  m <- 5L
  id <- factor(rep(seq_len(m), each = n / m))
  d <- data.frame(x = stats::runif(n, 0.2, 3), id = id)
  d$y <- 3 * exp(-0.7 * d$x) + stats::rnorm(n, sd = 0.15)

  calls <- list(
    quote(nl(~ a * exp(-r * x), start = list(a = 3, r = 0.7))),
    quote(nl(~ a * exp(-r * x), a ~ 0 + id)),
    quote(nl(~ a * exp(-r * x) + c, start = list(a = 3, r = 0.7, c = 0.1))),
    quote(nl(~ a * exp(-r * x), links = list(r = linkfunctions7::log_link()),
             start = list(a = 3, r = 0.7)))
  )
  for (call in calls) {
    tm <- term_build(eval(call), d)
    cf <- term_coef_start(tm)
    tm <- term_refresh(tm, cf)
    p <- length(term_coef_names(tm))
    set.seed(5)
    v <- stats::rnorm(p)
    A <- matrix(stats::rnorm(n * p), n, p)
    ex <- term_block_deriv(tm, cf, v)
    expect_identical(dim(ex), c(n, p))
    # the adjoint identity, which is exact arithmetic and not an approximation
    expect_equal(sum(A * ex), sum(v * term_block_contract(tm, cf, A)),
                 tolerance = 1e-12)
    # and against the block differenced along v, which shares nothing with
    # either: it is the block itself, evaluated twice
    h <- 1e-6
    bf <- (as.matrix(term_matrix(term_refresh(tm, cf + h * v))) -
             as.matrix(term_matrix(term_refresh(tm, cf - h * v)))) / (2 * h)
    expect_lt(max(abs(ex - bf)), 1e-6 * max(1, max(abs(bf))))
  }

  # a term whose block does not move answers zeros, of the block's shape
  lp <- term_build(linpar(~ x), d)
  z <- term_block_deriv(lp, v = c(0.3, -0.2))
  expect_identical(dim(z), dim(as.matrix(term_matrix(lp))))
  expect_true(all(z == 0))

  # and a direction of the wrong length is refused where it is given
  tm <- term_refresh(term_build(nl(~ a * exp(-r * x),
                                   start = list(a = 3, r = 0.7)), d),
                     c(3, 0.7))
  expect_error(term_block_deriv(tm, c(3, 0.7), v = 1), "length 2")
})


# ---------------------------------------------------------------------------
# term_coef_start(target =): a nonlinear term estimating its own parameters
# ---------------------------------------------------------------------------

nl_logistic_data <- function(seed = 456, ni = 30, nt = 20) {
  set.seed(seed)
  u1 <- stats::rnorm(ni, 0, 5); u2 <- stats::rnorm(ni, 0, 1.2)
  u3 <- stats::rnorm(ni, 0, 0.3)
  d <- data.frame(id = factor(rep(seq_len(ni), each = nt)),
                  time = rep(0:(nt - 1), ni))
  p1 <- abs(50 + u1[d$id]); p2 <- abs(10 + u2[d$id]); p3 <- abs(2 + u3[d$id])
  d$y <- p1 / (1 + exp(-(d$time - p2) / p3)) + stats::rnorm(nrow(d), 0, 1.4)
  d
}
nl_logistic_links <- function() {
  list(phi = linkfunctions7::log_link(),
       theta = linkfunctions7::identity_link(),
       sigma = linkfunctions7::log_link())
}
nl_theta_at <- function(b, coef) {
  bp <- b@blueprint
  vapply(bp$params, function(p) {
    e <- if (is.null(bp$Z[[p]])) coef[bp$index[[p]]] else
      as.numeric(as.matrix(bp$Z[[p]]) %*% coef[bp$index[[p]]])
    linkfunctions7::linkinv(bp$links[[p]], e)[1L]
  }, numeric(1))
}

test_that("without a target the start is what term_build computed", {
  d <- nl_logistic_data()
  b <- term_build(nl(~ phi / (1 + exp(-(time - theta) / sigma)),
                     links = nl_logistic_links()), d)
  expect_identical(term_coef_start(b), b@blueprint$coef)
  expect_true(all(term_coef_start(b) == 0))
})

test_that("a target estimates the parameters the caller did not pin", {
  d <- nl_logistic_data()
  lk <- nl_logistic_links()
  b <- term_build(nl(~ phi / (1 + exp(-(time - theta) / sigma)),
                     phi ~ 1 + lasso(~id), theta ~ 1 + lasso(~id),
                     sigma ~ 1 + lasso(~id), links = lk), d)
  th <- nl_theta_at(b, term_coef_start(b, target = d$y))
  expect_equal(unname(th[["phi"]]), 50, tolerance = 0.15)
  expect_equal(unname(th[["theta"]]), 10, tolerance = 0.15)
  expect_equal(unname(th[["sigma"]]), 2, tolerance = 0.25)
  # and the value it produces is far closer to the response than zero is
  sse_at <- function(v) sum((d$y - term_value(term_refresh(b, v), v))^2)
  expect_lt(sse_at(term_coef_start(b, target = d$y)),
            sse_at(term_coef_start(b)) / 10)
})

test_that("start= wins over a target, and pins only what it names", {
  d <- nl_logistic_data()
  lk <- nl_logistic_links()
  b <- term_build(nl(~ phi / (1 + exp(-(time - theta) / sigma)),
                     links = lk, start = list(phi = 50, theta = 10, sigma = 2)),
                  d)
  expect_identical(term_coef_start(b, target = d$y), b@blueprint$coef)

  b2 <- term_build(nl(~ phi / (1 + exp(-(time - theta) / sigma)),
                      links = lk, start = list(theta = 10)), d)
  th <- nl_theta_at(b2, term_coef_start(b2, target = d$y))
  expect_equal(unname(th[["theta"]]), 10)          # pinned, exactly
  expect_equal(unname(th[["phi"]]), 50, tolerance = 0.15)
  expect_equal(unname(th[["sigma"]]), 2, tolerance = 0.25)
})

test_that("the start is deterministic and ignores the caller's seed", {
  d <- nl_logistic_data()
  b <- term_build(nl(~ phi / (1 + exp(-(time - theta) / sigma)),
                     links = nl_logistic_links()), d)
  set.seed(1); a1 <- term_coef_start(b, target = d$y)
  set.seed(99); a2 <- term_coef_start(b, target = d$y)
  a3 <- term_coef_start(b, target = d$y)
  expect_identical(a1, a2)
  expect_identical(a1, a3)
})

test_that("a constant reaches every column of a full-rank development", {
  set.seed(2)
  d <- data.frame(x = seq(0.2, 5, length.out = 200),
                  g = factor(rep(1:4, each = 50)))
  d$y <- 3 * exp(-0.8 * d$x) + stats::rnorm(200, 0, 0.05)
  lk <- list(a = linkfunctions7::log_link(), r = linkfunctions7::log_link())
  b <- term_build(nl(~ a * exp(-r * x), a ~ 0 + lasso(~g), links = lk), d)
  v <- term_coef_start(b, target = d$y)
  ia <- grep("[.]a[.]", term_coef_names(b))
  expect_length(ia, 4L)
  # the development carries no intercept, so the constant is every column
  expect_equal(max(v[ia]) - min(v[ia]), 0, tolerance = 1e-10)
  expect_equal(exp(v[ia][1L]), 3, tolerance = 0.05)
})

test_that("an opaque function gets no symbolic separation and still works", {
  set.seed(3)
  d <- data.frame(x = seq(0.2, 5, length.out = 200))
  d$y <- 3 * exp(-0.8 * d$x) + stats::rnorm(200, 0, 0.05)
  lk <- list(a = linkfunctions7::log_link(), r = linkfunctions7::log_link())
  fn <- function(x, theta) theta$a * exp(-theta$r * x)
  b <- term_build(nl(fn, params = c("a", "r"), x = x, links = lk), d)
  v <- term_coef_start(b, target = d$y)
  expect_equal(exp(v[[1L]]), 3, tolerance = 0.05)
  expect_equal(exp(v[[2L]]), 0.8, tolerance = 0.05)
})

test_that("a target that says nothing leaves the built coefficients", {
  d <- nl_logistic_data()
  b <- term_build(nl(~ phi / (1 + exp(-(time - theta) / sigma)),
                     links = nl_logistic_links()), d)
  expect_identical(term_coef_start(b, target = rep(NA_real_, nrow(d))),
                   b@blueprint$coef)
  expect_identical(term_coef_start(b, target = d$y[1:3]), b@blueprint$coef)
})

test_that("the affine parameters are read off the expression", {
  expect_identical(.nl_affine_params(quote(a * exp(-r * x)), c("a", "r")), "a")
  expect_identical(.nl_affine_params(quote(a + b * x), c("a", "b")),
                   c("a", "b"))
  expect_identical(
    .nl_affine_params(quote(phi / (1 + exp(-(time - theta) / sigma))),
                      c("phi", "theta", "sigma")), "phi")
  # d0 and a are affine jointly, e and b are not
  expect_setequal(
    .nl_affine_params(quote(d0 + (a - d0) / (1 + exp(-(x - e) / b))),
                      c("a", "d0", "e", "b")), c("a", "d0"))
})
