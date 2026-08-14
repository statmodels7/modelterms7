# The penalized quartet. Each term's penalty is pinned against the
# penalties7 constructor called directly: the same object, so no tolerance
# has to be chosen.

set.seed(7)
dd <- data.frame(x1 = rnorm(8), x2 = rnorm(8),
                 g = factor(rep(c("a", "b"), 4)))
dd$R <- matrix(rnorm(24), 8, 3, dimnames = list(NULL, c("r1", "r2", "r3")))

test_that("the four terms attach their penalty and read smoothness from it", {
  cases <- list(
    list(ctor = ridge, ref = penalties7::ridge_penalty, smooth = TRUE),
    list(ctor = lasso, ref = penalties7::lasso_penalty, smooth = FALSE),
    list(ctor = scad, ref = penalties7::scad_penalty, smooth = FALSE),
    list(ctor = mcp, ref = penalties7::mcp_penalty, smooth = FALSE)
  )
  # FOUR columns: the intercept is dropped, so a two-level factor is coded
  # full rank. The separable branch recycled a shorter vector without a word,
  # which is how this read three of them for as long as the ridge was on it.
  beta <- c(0.6, -1.1, 2.0, -0.3)
  for (case in cases) {
    built <- term_build(case$ctor(~ x1 + x2 + g), dd)
    pen <- term_penalty(built)
    ref <- case$ref(n_coef = term_npar(built))
    expect_identical(pen@params, ref@params)
    th <- lapply(pen@params_bounds, function(b) {
      if (all(is.finite(b))) mean(b) else if (is.finite(b[1])) b[1] + 1 else 1
    })
    expect_identical(penalties7::penalty_value(pen, beta, th),
                     penalties7::penalty_value(ref, beta, th))
    expect_identical(term_smooth(built), case$smooth)
  }
})

test_that("a formula input removes the intercept and keeps the blueprint discipline", {
  built <- term_build(ridge(~ x1 + x2), dd)
  ref <- stats::model.matrix(~ x1 + x2 - 1, dd)
  expect_equal(unname(term_matrix(built)), unname(ref), ignore_attr = TRUE)
  expect_identical(term_coef_names(built), c("ridge.x1", "ridge.x2"))
  res <- check_term(ridge(~ x1 + x2 + g), dd, verbose = FALSE)
  expect_true(all(res$status == "OK"),
              info = paste(res$check[res$status != "OK"], collapse = ", "))
})

test_that("a matrix column of the data builds, predicts and validates", {
  spec <- interpret_formula(~ lasso(R), dd)$terms[["lasso(R)"]]
  built <- term_build(spec, dd)
  expect_identical(term_coef_names(built),
                   c("lasso.r1", "lasso.r2", "lasso.r3"))
  expect_equal(unname(term_matrix(built)), unname(dd$R), ignore_attr = TRUE)

  # prediction re-evaluates the expression in the new data, so a subset of
  # the data carries the matching rows
  sub <- dd[c(2, 5, 7), , drop = FALSE]
  expect_equal(unname(term_predict(built, sub)), unname(dd$R[c(2, 5, 7), ]),
               ignore_attr = TRUE)

  for (lb in c("ridge(R)", "lasso(R)", "scad(R)", "mcp(R)")) {
    f <- stats::as.formula(paste("~", lb))
    sp <- interpret_formula(f, dd)$terms[[lb]]
    res <- check_term(sp, dd, verbose = FALSE)
    expect_true(all(res$status == "OK"),
                info = paste(lb, paste(res$check[res$status != "OK"],
                                       collapse = ", ")))
  }
})

test_that("a free-standing matrix builds but refuses a mismatched prediction", {
  M <- matrix(rnorm(16), 8, 2)
  built <- term_build(ridge(M), dd)
  expect_identical(term_coef_names(built), c("ridge.1", "ridge.2"))
  # the expression 'M' does not resolve in newdata, and the build-time rows
  # are deliberately not reused
  expect_error(term_predict(built, dd[1:3, , drop = FALSE]),
               "column of the data")
  # and a matrix whose rows do not match the data cannot build at all
  expect_error(term_build(ridge(M[1:5, ]), dd), "5 rows")
})

test_that("by is reserved and inputs are validated", {
  expect_error(ridge(~x1, by = dd$g), "reserved")
  expect_error(lasso(y ~ x1), "one-sided")
  expect_error(scad(matrix("a", 2, 2)), "numeric")
  expect_error(mcp(~x1, label = ""), "non-empty")
})

test_that("the formula interpreter routes the penalized terms", {
  out <- interpret_formula(x1 ~ x2 + ridge(~g, label = "rg") + lasso(R), dd)
  expect_named(out$terms,
               c("linpar", 'ridge(~g, label = "rg")', "lasso(R)"))
  built <- term_build(out$terms[["lasso(R)"]], dd)
  expect_identical(term_npar(built), 3L)
})

test_that("enet() carries the elastic-net penalty and both hyperparameters", {
  spec <- interpret_formula(~ enet(R), dd)$terms[["enet(R)"]]
  b <- term_build(spec, dd)
  expect_identical(term_coef_names(b), c("enet.r1", "enet.r2", "enet.r3"))
  pen <- term_penalty(b)
  expect_identical(pen@params, c("lambda", "alpha"))
  expect_identical(pen@n_coef, 3L)
  # it mixes a kink with a quadratic, so the penalized objective is not
  # differentiable and the term says so
  expect_false(term_smooth(b))

  # the value is the elastic net of the block's coefficients
  cf <- c(0.9, -1.4, 0.2)
  th <- list(lambda = 1.1, alpha = 0.35)
  got <- penalties7::penalty_value(pen, cf, th) -
    penalties7::penalty_value(pen, rep(0, 3), th)
  expect_equal(got, 1.1 * (0.35 * sum(abs(cf)) + 0.65 * sum(cf^2) / 2),
               tolerance = 1e-12)

  # and edf counts the nonzero coefficients, as it does for the lasso
  expect_identical(edf(b, coef = c(0.9, -1.4, 0)), 2)

  res <- check_term(spec, dd, verbose = FALSE)
  expect_true(all(res$status == "OK"),
              info = paste(res$check[res$status != "OK"], collapse = ", "))
})

test_that("the interpreter routes enet beside the other four", {
  out <- interpret_formula(x1 ~ x2 + enet(R), dd)
  expect_named(out$terms, c("linpar", "enet(R)"))
  b <- term_build(out$terms[["enet(R)"]], dd)
  expect_output(print(b), "enet")
})


test_that("a term reports every penalty it carries", {
  # term_penalty() answers for one penalty over a whole design block, which is
  # every term shipped here. term_penalties() generalizes it in the two
  # directions a model layer needs: several penalties on one term, and
  # penalties over parameters that are not coefficients of a block -- the
  # persistence of a gas(), the nonlinear parameters of nl(), the break-point
  # of seg(). The base method answers from term_penalty(), so nothing here
  # needs a method of its own and nothing downstream sees a change.
  d <- data.frame(x = stats::rnorm(30), z = stats::rnorm(30),
                  g = factor(rep(1:3, 10)))
  for (tm in list(ridge(~x), lasso(~ x + z), s(x, k = 5), random(~ 1 | g))) {
    b <- term_build(tm, d)
    ps <- term_penalties(b)
    expect_length(ps, 1L)
    expect_identical(ps[[1L]]$name, "")
    expect_identical(ps[[1L]]$index, seq_len(term_npar(b)))
    expect_identical(ps[[1L]]$penalty, term_penalty(b))
  }
  # an unpenalized term carries none
  expect_length(term_penalties(term_build(linpar(~x), d)), 0L)
  # and a structural term answers rather than raising, which is what lets a
  # caller enumerate over every term without knowing which kind it has
  expect_length(term_penalties(gas(p = 1, q = 1)), 0L)

  # the name is the entry's WITHIN the term and never the term's own: two
  # ridge() terms in one formula are two terms with two hyperparameters, and
  # it is the caller that knows what it called each one
  a <- term_build(ridge(~x), d)
  b2 <- term_build(ridge(~z), d)
  expect_identical(term_penalties(a)[[1L]]$name,
                   term_penalties(b2)[[1L]]$name)
})


test_that("a sparse matrix input is never densified", {
  # a penalized block is exactly where a sparse design turns up -- indicators
  # over many levels are what a lasso is for -- and densifying it at the
  # constructor would undo the whole saving, by a factor of 1/density.
  set.seed(9)
  n <- 400L; k <- 20L
  R <- Matrix::sparseMatrix(i = sample(n, n), j = sample(k, n, TRUE), x = 1,
                            dims = c(n, k))
  dat <- data.frame(y = rnorm(n))
  dat$R <- R
  for (ctor in list(ridge, lasso, enet, scad, mcp)) {
    spec <- ctor(R)
    expect_s4_class(spec@input, "dgCMatrix")
    built <- term_build(spec, dat)
    X <- term_matrix(built)
    expect_s4_class(X, "dgCMatrix")
    expect_identical(dim(X), c(n, k))
    # and it stays sparse at new data, where a densifying predict would
    # spend what the build was careful not to
    expect_s4_class(term_predict(built, dat[1:50, , drop = FALSE]),
                    "dgCMatrix")
  }
  # the values are the ones a dense build gives
  bs <- term_build(lasso(R), dat)
  Rd <- as.matrix(R)
  dat2 <- data.frame(y = dat$y); dat2$R <- Rd
  bd <- term_build(lasso(Rd), dat2)
  expect_equal(as.matrix(term_matrix(bs)), term_matrix(bd),
               ignore_attr = TRUE)
  expect_true(all(check_term(lasso(R), dat, verbose = FALSE)$status == "OK"))

  # a logical Matrix -- the commonest sparse input of all -- is carried to
  # double rather than rejected
  L <- Matrix::sparseMatrix(i = sample(n, n), j = sample(k, n, TRUE),
                            dims = c(n, k))
  expect_true(methods::is(lasso(L)@input, "dMatrix"))

  # a standardized sparse block reads its spreads without densifying
  bstd <- term_build(lasso(R, standardize = TRUE), dat)
  expect_s4_class(term_matrix(bstd), "dgCMatrix")
  expect_equal(unname(bstd@blueprint$standardize),
               unname(apply(as.matrix(R), 2L, stats::sd)))
})


# --- standardization --------------------------------------------------------

test_that("standardize attaches the column spreads as a diagonal map", {
  set.seed(11)
  d <- data.frame(x1 = rnorm(60, mean = 100, sd = 3),
                  x2 = rnorm(60) * 1000)
  d$k <- 4                                  # a constant column
  for (ctor in list(ridge, lasso, enet, scad, mcp)) {
    built <- term_build(ctor(~ x1 + x2 + k, standardize = TRUE), d)
    s <- built@blueprint$standardize
    expect_identical(names(s), term_coef_names(built))
    expect_equal(unname(s[1:2]), c(stats::sd(d$x1), stats::sd(d$x2)))
    # a constant column has no spread to divide by and keeps its own scale
    expect_identical(unname(s[[3L]]), 1)
    map <- term_penalty(built)@map
    expect_true(methods::is(map, "diagonalMatrix"))
    expect_equal(as.numeric(Matrix::diag(map)), unname(s))
  }
  # and without the argument there is no map at all
  expect_null(term_penalty(term_build(lasso(~ x1 + x2), d))@map)
  expect_null(term_build(lasso(~ x1 + x2), d)@blueprint$standardize)
})

test_that("the standardized penalty is the penalty on a standardized design", {
  # the whole claim: rho(S beta_x) with S = diag(sd) is what a penalty on the
  # columns divided by their own spread comes to, so the two routes agree at
  # machine precision and the design is never touched
  set.seed(12)
  d <- data.frame(x1 = rnorm(50, 10, 2), x2 = rnorm(50) * 500)
  built <- term_build(lasso(~ x1 + x2, standardize = TRUE), d)
  s <- unname(built@blueprint$standardize)
  X <- term_matrix(built)
  Z <- sweep(X, 2L, s, "/")
  expect_equal(unname(apply(Z, 2L, stats::sd)), c(1, 1))

  beta_x <- c(0.7, -0.002)
  th <- list(lambda = 1.4)
  pen_std <- term_penalty(built)
  pen_raw <- penalties7::lasso_penalty(n_coef = 2L)
  # beta_z = s beta_x is the same fitted function on the standardized design
  expect_equal(as.numeric(X %*% beta_x), as.numeric(Z %*% (s * beta_x)))
  expect_equal(penalties7::penalty_value(pen_std, beta_x, th),
               penalties7::penalty_value(pen_raw, s * beta_x, th))
})

test_that("standardizing SCAD and MCP composes BOTH hyperparameters", {
  # the trap: substituting s b is not "multiply lambda by s". Measured on
  # the published piecewise forms, transcribed here and sharing no code with
  # penalties7, the naive substitution is out by 4.4 for MCP and 11.2 for
  # SCAD, while the exact relations carry an overall factor of s^2 as well.
  rho_scad <- function(b, lam, a) {
    t <- abs(b)
    ifelse(t <= lam, lam * t,
      ifelse(t <= a * lam,
        (2 * a * lam * t - t^2 - lam^2) / (2 * (a - 1)),
        lam^2 * (a + 1) / 2))
  }
  rho_mcp <- function(b, lam, gam) {
    t <- abs(b)
    ifelse(t <= gam * lam, lam * t - t^2 / (2 * gam), gam * lam^2 / 2)
  }
  b <- c(-2.3, -0.4, 0.15, 1.7, 3.9)
  s <- c(0.5, 0.8, 1.25, 2, 3)
  lam <- 1.3; a <- 3.7; gam <- 2.6
  D <- Matrix::Diagonal(x = s)

  # what the map computes IS rho(s b), exactly
  expect_equal(penalties7::penalty_value(penalties7::scad_penalty(map = D), b,
                                         list(lambda = lam, a = a)),
               sum(rho_scad(s * b, lam, a)))
  expect_equal(penalties7::penalty_value(penalties7::mcp_penalty(map = D), b,
                                         list(lambda = lam, gamma = gam)),
               sum(rho_mcp(s * b, lam, gam)))

  # SCAD is NOT a SCAD at rescaled parameters: the slope near zero asks for
  # lambda*s and the first breakpoint for lambda/s, and only an overall
  # factor reconciles them
  expect_equal(rho_scad(s * b, lam, a), s^2 * rho_scad(b, lam / s, a))
  expect_gt(max(abs(rho_scad(s * b, lam, a) - rho_scad(b, lam * s, a))), 10)

  # MCP does have a member of its own family, and its second parameter moves
  # with the SQUARE of the scale, not the scale
  expect_equal(rho_mcp(s * b, lam, gam), rho_mcp(b, lam * s, gam / s^2))
  expect_equal(rho_mcp(s * b, lam, gam), s^2 * rho_mcp(b, lam / s, gam))
  expect_gt(max(abs(rho_mcp(s * b, lam, gam) - rho_mcp(b, lam * s, gam / s))), 4)
})

test_that("a standardized term keeps its proximal table", {
  # the route a compiled coordinate descent takes: without a table it would
  # fall back on the general operator and pay an R call per coordinate
  set.seed(13)
  d <- data.frame(x1 = rnorm(40), x2 = rnorm(40) * 100)
  for (ctor in list(lasso, enet, scad, mcp)) {
    pen <- term_penalty(term_build(ctor(~ x1 + x2, standardize = TRUE), d))
    expect_true(penalties7::has_prox(pen))
    th <- lapply(pen@params_bounds, function(b) {
      if (all(is.finite(b))) mean(b) else if (is.finite(b[1])) b[1] + 1 else 1
    })
    sp <- penalties7::penalty_prox_spec(pen, th, rep(1e-4, 2))
    expect_false(is.null(sp))
    expect_equal(penalties7::prox_apply(sp, c(1.5, -0.4)),
                 penalties7::penalty_prox(pen, c(1.5, -0.4), 1e-4, th),
                 tolerance = 1e-12)
  }
})

test_that("a term with no columns rejects standardize by name", {
  # a deviation is a parameter of a recursion and a random effect's columns
  # are grouping indicators: neither has a spread to divide by, and neither
  # constructor takes the argument, so the request is an error rather than
  # something silently ignored
  expect_error(gas(p = 1, q = 1, by = g, standardize = TRUE), "standardize")
  expect_error(random(~ 1 | g, standardize = TRUE), "standardize")
  expect_error(lasso(~x, standardize = "yes"), "TRUE or FALSE")
  expect_error(lasso(~x, standardize = NA), "TRUE or FALSE")
})

test_that("print says what a term was standardized by", {
  set.seed(14)
  d <- data.frame(x1 = rnorm(30) * 7)
  built <- term_build(lasso(~ x1, standardize = TRUE), d)
  expect_output(print(built), "standardized by")
  expect_output(print(lasso(~x1, standardize = TRUE)), "standardized")
  # and an unstandardized term prints no such line
  out <- utils::capture.output(print(term_build(lasso(~x1), d)))
  expect_false(any(grepl("standardized", out)))
})
