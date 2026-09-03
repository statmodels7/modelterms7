# Grouped random effects: intercepts, slopes, and the effect distributions.

set.seed(21)
dd <- data.frame(y = rnorm(12), x = rnorm(12),
                 g = factor(rep(c("a", "b", "c"), 4)))

.centered <- function(d, ...) distributions7::fixed(d, ...)

.mv_centered <- function(p, ..., inverted = FALSE) {
  ctor <- if (inverted) {
    distributions7::mvgaussian2_distrib
  } else {
    distributions7::mvgaussian1_distrib
  }
  do.call(distributions7::fixed,
          c(list(ctor(p, ...)),
            stats::setNames(as.list(rep(0, p)), paste0("mu", seq_len(p)))))
}

test_that("random intercepts build the indicator block with a gaussian default", {
  built <- term_build(random(~ 1 | g), dd)
  ref <- stats::model.matrix(~ 0 + g, dd)
  # the block is SPARSE by construction: a row belongs to one group, so the
  # density is 1/m whatever the data, and the dense form was the whole cost
  # of a random effect. It is compared against the dense reference by value
  Z <- term_matrix(built)
  expect_s4_class(Z, "sparseMatrix")
  expect_equal(unname(as.matrix(Z)), unname(ref), ignore_attr = TRUE)
  expect_identical(term_coef_names(built),
                   c("random.a", "random.b", "random.c"))

  # the hyperparameter IS the standard deviation of the effects: a reader of
  # a mixed model gets the variance component and not a precision to invert
  pen <- term_penalty(built)
  expect_identical(pen@params, "sigma")
  expect_true(term_smooth(built))

  res <- check_term(random(~ 1 | g), dd, verbose = FALSE)
  expect_true(all(res$status == "OK"),
              info = paste(res$check[res$status != "OK"], collapse = ", "))
})

test_that("the simple random effect IS a ridge, at the matching hyperparameter", {
  # they remain the same model and not the same object: the ridge reports a
  # precision and this reports a standard deviation, so the two are pinned
  # against each other rather than one being built from the other
  built <- term_build(random(~ 1 | g), dd)
  pen <- term_penalty(built)
  b <- c(0.4, -1.1, 2.2)
  s <- 1.7
  ridge <- penalties7::ridge_penalty(n_coef = 3)
  expect_equal(
    penalties7::penalty_gradient(pen, b, list(sigma = s)),
    penalties7::penalty_gradient(ridge, b, list(lambda = 1 / s^2)),
    tolerance = 1e-12)
  expect_equal(
    unname(penalties7::penalty_hessian(pen, b, list(sigma = s))),
    unname(penalties7::penalty_hessian(ridge, b, list(lambda = 1 / s^2))),
    tolerance = 1e-12)
  # the values differ by the constant that makes sigma estimable at all
  expect_equal(penalties7::penalty_value(pen, b, list(sigma = s)) -
               penalties7::penalty_value(ridge, b, list(lambda = 1 / s^2)),
               3 * log(s * sqrt(2 * pi)) -
               penalties7::penalty_value(ridge, rep(0, 3),
                                         list(lambda = 1 / s^2)),
               tolerance = 1e-12)
})

test_that("random slopes interact the within-group design with the groups", {
  built <- term_build(random(~ x | g), dd)
  expect_identical(term_npar(built), 6L)
  expect_identical(term_coef_names(built),
                   c("random.a.(Intercept)", "random.a.x",
                     "random.b.(Intercept)", "random.b.x",
                     "random.c.(Intercept)", "random.c.x"))

  # reference construction: indicator times within-group column, group-major
  Z <- as.matrix(term_matrix(built))
  G <- stats::model.matrix(~ 0 + g, dd)
  expect_equal(unname(Z[, 1]), unname(G[, 1]), ignore_attr = TRUE)
  expect_equal(unname(Z[, 2]), unname(G[, 1] * dd$x), ignore_attr = TRUE)
  expect_equal(unname(Z[, 6]), unname(G[, 3] * dd$x), ignore_attr = TRUE)
  # the block reproduces the fixed design when the group effects are summed
  expect_equal(rowSums(Z[, c(1, 3, 5)]), rep(1, 12))
  expect_equal(rowSums(Z[, c(2, 4, 6)]), dd$x)

  res <- check_term(random(~ x | g), dd, verbose = FALSE)
  expect_true(all(res$status == "OK"),
              info = paste(res$check[res$status != "OK"], collapse = ", "))

  # and the slope alone drops the intercept by the formula convention
  b0 <- term_build(random(~ 0 + x | g), dd)
  expect_identical(term_npar(b0), 3L)
  expect_identical(term_coef_names(b0),
                   c("random.a.x", "random.b.x", "random.c.x"))
})

test_that("the correlated default is one blockwise multivariate gaussian", {
  bc <- term_build(random(~ x | g), dd)
  ent <- term_penalties(bc)
  expect_length(ent, 1L)
  expect_identical(ent[[1L]]$name, "")
  expect_identical(term_penalty(bc)@params,
                   c("sigma_log_L1", "sigma_log_L2", "sigma_L2.1"))
  expect_identical(term_penalty(bc)@block, 2L)

  # and it IS the penalty penalties7 builds from the same centered family
  ref <- penalties7::distrib_penalty(
    .mv_centered(2, sigma = parameters7::log_cholesky(
      2)), n_coef = 6)
  beta <- stats::rnorm(6)
  th <- stats::setNames(as.list(c(0.1, -0.2, 0.3)), ref@params)
  expect_equal(penalties7::penalty_value(term_penalty(bc), beta, th),
               penalties7::penalty_value(ref, beta, th))
  expect_true(term_smooth(bc))

  # edf runs on the block
  H <- as.matrix(crossprod(term_matrix(bc)))
  e <- edf(bc, coef = beta, hessian = H, theta = th)
  expect_true(is.finite(e) && e > 0 && e < 6)
})

test_that("correlated = FALSE is one penalty PER within-group column", {
  # an intercept and a slope are quantities of different units, so they get a
  # standard deviation each rather than sharing one
  bu <- term_build(random(~ x | g, correlated = FALSE), dd)
  ent <- term_penalties(bu)
  expect_length(ent, 2L)
  expect_identical(vapply(ent, function(e) e$name, ""),
                   c("(Intercept)", "x"))
  # the coefficients are ordered group by group, so column j is the stride
  # j, d+j, 2d+j
  expect_identical(ent[[1L]]$index, c(1L, 3L, 5L))
  expect_identical(ent[[2L]]$index, c(2L, 4L, 6L))
  for (e in ent) expect_identical(e$penalty@params, "sigma")
  # term_penalty() answers for one penalty over the whole block and there is
  # not one here
  expect_null(term_penalty(bu))
})

test_that("a multivariate prior carries the dependence, and says which matrix", {
  st <- .mv_centered(2, parameters7::ar1(2), inverted = TRUE)
  built <- term_build(random(~ x | g, distrib = st), dd)
  # the free name says how the matrix is built AND which matrix it is
  expect_identical(term_penalty(built)@params,
                   c("omega_log_scale", "omega_z_rho"))
  expect_identical(term_penalty(built)@block, 2L)

  expect_error(term_build(random(~ 1 | g, distrib = st), dd), "1 columns")
  res <- check_term(random(~ x | g, distrib = st), dd, verbose = FALSE)
  expect_true(all(res$status == "OK"),
              info = paste(res$check[res$status != "OK"], collapse = ", "))
})

test_that("a univariate prior over several columns is a template", {
  ph <- .centered(distributions7::pseudohuber_distrib(), mu = 0)
  b1 <- term_build(random(~ 1 | g, distrib = ph), dd)
  expect_identical(term_penalty(b1)@params, c("sigma", "nu"))

  b2 <- term_build(random(~ x | g, distrib = ph), dd)
  ent <- term_penalties(b2)
  expect_length(ent, 2L)
  for (e in ent) expect_identical(e$penalty@params, c("sigma", "nu"))
  # each copy has its OWN hyperparameters: four in all, not two
  expect_length(unlist(lapply(ent, function(e) e$penalty@params)), 4L)

  # and a list gives one per column explicitly
  lap <- .centered(distributions7::laplace_distrib(), mu = 0)
  b3 <- term_build(random(~ x | g, distrib = list(ph, lap)), dd)
  ent3 <- term_penalties(b3)
  expect_identical(ent3[[1L]]$penalty@params, c("sigma", "nu"))
  expect_identical(ent3[[2L]]$penalty@params, "sigma")
  expect_error(term_build(random(~ x | g, distrib = list(ph)), dd),
               "one distribution per column")
})

test_that("the kinks of the effects' prior are DERIVED, not declared here", {
  # the default used to be numeric(0), which overrode penalties7's own
  # derivation, so a Laplace prior declared none and a fitting layer sent its
  # block to the scheme that cannot solve it
  lap <- .centered(distributions7::laplace_distrib(), mu = 0)
  built <- term_build(random(~ 1 | g, distrib = lap), dd)
  expect_identical(penalties7::penalty_kinks(term_penalty(built),
                                             list(sigma = 1)), 0)
  expect_false(term_smooth(built))

  # a smooth prior has none, and says so for itself
  ph <- .centered(distributions7::pseudohuber_distrib(), mu = 0)
  expect_true(term_smooth(term_build(random(~ 1 | g, distrib = ph), dd)))
})

test_that("a hyperparameter is held by name, qualified by its column", {
  ph <- .centered(distributions7::pseudohuber_distrib(), mu = 0)
  b1 <- term_build(random(~ 1 | g, distrib = ph, hyper = c(nu = 2)), dd)
  expect_equal(term_hyper(b1)[[1L]]$nu, 2)

  b2 <- term_build(random(~ x | g, distrib = ph,
                          hyper = c("nu.(Intercept)" = 2, "sigma.x" = 0.5)), dd)
  h <- term_hyper(b2)
  expect_equal(h[["(Intercept)"]]$nu, 2)
  expect_equal(h[["x"]]$sigma, 0.5)

  # an UNQUALIFIED name is an error listing what there is, not a value
  # recycled over every column
  expect_error(term_build(random(~ x | g, distrib = ph, hyper = c(nu = 2)), dd),
               "nu.\\(Intercept\\)")
  expect_error(term_build(random(~ 1 | g, hyper = c(lambda = 4)), dd),
               "not a hyperparameter")
})

test_that("holding inside the distribution and holding by name agree", {
  # the two differ in what is REPORTED and not in the fit: one removes the
  # hyperparameter from the model, the other keeps it and holds it
  ph <- .centered(distributions7::pseudohuber_distrib(), mu = 0)
  ph2 <- .centered(distributions7::pseudohuber_distrib(), mu = 0, nu = 2)
  pa <- term_penalty(term_build(random(~ 1 | g, distrib = ph2), dd))
  pb <- term_penalty(term_build(random(~ 1 | g, distrib = ph,
                                       hyper = c(nu = 2)), dd))
  expect_identical(pa@params, "sigma")
  expect_identical(pb@params, c("sigma", "nu"))
  b <- c(0.3, -0.8, 1.1)
  expect_equal(penalties7::penalty_value(pa, b, list(sigma = 0.8)),
               penalties7::penalty_value(pb, b, list(sigma = 0.8, nu = 2)))
  expect_equal(penalties7::penalty_gradient(pa, b, list(sigma = 0.8)),
               penalties7::penalty_gradient(pb, b, list(sigma = 0.8, nu = 2)))
})

test_that("a prior on the effects is centered, and a free mean is rejected", {
  # a free mean in the effects is confounded with the intercept of the
  # equation the term sits in, which is a flat direction and not a model
  expect_error(
    term_build(random(~ 1 | g,
                      distrib = distributions7::gaussian1_distrib()), dd),
    "free location")
  # a location HELD is identified whatever its value and is not policed: it
  # shrinks the effects towards that value, which is a modelling statement.
  # Nor could the value be policed in general -- where the prior is a
  # transformation of another family the parameter is the mean on the
  # ORIGINAL scale, and holding a gamma's mean at one is what centers its
  # logarithm, at a value zero would put outside the parameter's own domain.
  expect_no_error(
    term_build(random(~ 1 | g,
                      distrib = .centered(distributions7::gaussian1_distrib(),
                                          mu = 1)), dd))
  lg <- .centered(distributions7::transformation(
    distributions7::gamma2_distrib(), distributions7::log_transform()),
    mu = 1)
  expect_identical(
    term_penalty(term_build(random(~ 1 | g, distrib = lg), dd))@params,
    "sigma2")
  # and its free-location form is still refused, the interpretation string
  # carrying the parent's scale rather than being the bare word
  expect_error(
    term_build(random(~ 1 | g, distrib = distributions7::transformation(
      distributions7::gamma2_distrib(),
      distributions7::log_transform())), dd),
    "free location")
})

test_that("correlation is refused for a family that cannot express it", {
  # the question is a PROPERTY -- a location block as long as the dimension
  # and a matrix parameter -- and not a list of two admitted names
  expect_error(
    term_build(random(~ x | g,
                      distrib = distributions7::dirichlet_distrib(2)), dd),
    "no matrix parameter")
  expect_error(
    term_build(random(~ 1 | g,
                      distrib = .centered(distributions7::poisson_distrib(),
                                          mu = 1)), dd),
    "discrete")
})

test_that("levels are pinned at build and unknown levels are rejected", {
  dc <- data.frame(x = rnorm(6), g = rep(c("u", "v"), 3))
  built <- term_build(random(~ x | g), dc)
  sub <- data.frame(x = 0.5, g = "u")
  expect_identical(dim(term_predict(built, sub)), c(1L, 4L))
  expect_error(term_predict(built, data.frame(x = 1, g = "w")),
               "not present at build time")
})

test_that("degenerate inputs are rejected", {
  expect_error(random(~ g), "grouping bar")
  expect_error(random(y ~ 1 | g), "one-sided")
  expect_error(random(~ 1 | g, correlated = NA), "TRUE or FALSE")
  expect_error(term_build(random(~ 1 | g), data.frame(g = rep("a", 4))),
               "at least two levels")
})

test_that("the removed arguments are reported by name", {
  # a removed argument lands in the dots, where it would be swallowed; it is
  # named instead, with the spelling that replaces it
  expect_error(random(~ x | g, precision = parameters7::ar1(2)),
               "mvgaussian2_distrib")
  expect_error(random(~ 1 | g, kinks = 0), "distrib_kinks")
  expect_error(random(~ 1 | g, standardize = TRUE), "no argument")
  # and two ways of saying the same thing is an error, not a silent winner
  ph <- .centered(distributions7::pseudohuber_distrib(), mu = 0)
  expect_error(random(~ x | g, distrib = ph, correlated = FALSE), "Say it once")
})

test_that("the formula interpreter routes random()", {
  out <- interpret_formula(y ~ x + random(~ x | g), dd)
  expect_named(out$terms, c("linpar", "random(~x | g)"))
  built <- term_build(out$terms[["random(~x | g)"]], dd)
  expect_identical(term_npar(built), 6L)
})


test_that("a multivariate Student t prior is admitted", {
  # it carries a location block and a scale matrix, and it answers its mixed
  # response-parameter block, which is what the admissibility rule asks
  mvt <- do.call(distributions7::fixed,
                 list(distributions7::mvstudent_t1_distrib(2),
                      mu1 = 0, mu2 = 0))
  built <- term_build(random(~ x | g, distrib = mvt), dd)
  expect_identical(term_penalty(built)@block, 2L)
  expect_identical(term_penalty(built)@params,
                   c("sigma_log_L1", "sigma_log_L2", "sigma_L2.1", "nu"))
  res <- check_term(random(~ x | g, distrib = mvt), dd, verbose = FALSE)
  expect_true(all(res$status == "OK"),
              info = paste(res$check[res$status != "OK"], collapse = ", "))
  # and its degrees of freedom are held like any other hyperparameter
  h <- term_build(random(~ x | g, distrib = mvt, hyper = c(nu = 5)), dd)
  expect_equal(term_hyper(h)[[1L]]$nu, 5)
})


# The covariance label, and the grouping variable's type -------------------

test_that("the second bar is read as a covariance label", {
  # R nests bars to the left, so the last position is the grouping variable
  # and the middle one the label: brms's convention, out of the parser.
  expect_identical(term_tag(random(~ 1 | u | g)), "u")
  expect_identical(term_tag(random(~ 1 + x | u | g)), "u")
  # and a term without one says so, as does any other kind of term
  expect_identical(term_tag(random(~ 1 | g)), NA_character_)
  expect_identical(term_tag(linpar(~ x)), NA_character_)
  expect_identical(term_tag(s(x)), NA_character_)
})

test_that("the label is a name and is never evaluated", {
  # This is what the two-bar form did before it had a meaning: with no such
  # column it stopped with "object 'u' not found"; with a factor column it
  # built twice the columns with a warning about '|' applied to factors; with
  # a NUMERIC column it built twice the columns and said nothing at all. All
  # three now build the block the one-bar formula builds.
  ref <- term_build(random(~ 1 + x | g), dd)
  for (u in list(NULL, factor(rep(c("p", "q"), 6)), as.numeric(1:12))) {
    d <- dd
    if (!is.null(u)) d$u <- u
    b <- expect_silent(term_build(random(~ 1 + x | u | g), d))
    expect_identical(term_coef_names(b), term_coef_names(ref))
    expect_identical(as.matrix(term_matrix(b)), as.matrix(term_matrix(ref)))
    # the block is the one-bar formula's; the PENALTY is not, and that is the
    # label's whole effect -- it belongs to the class and not to this term
    expect_null(term_penalty(b))
  }
})

test_that("the label must be a symbol, and three bars are refused", {
  expect_error(random(~ 1 | "u" | g), "must be a name")
  expect_error(random(~ 1 | u + v | g), "must be a name")
  expect_error(random(~ 1 | f(u) | g), "must be a name")
  expect_error(random(~ 1 | u | v | g), "three bars")
})

test_that("a continuous grouping variable is refused", {
  # Without the check, factor() turns a continuous covariate into one level
  # per distinct value: measured, random(~ 1 | x) at 60 observations built 60
  # columns and reported nothing.
  d <- data.frame(x = rnorm(60), cont = rnorm(60),
                  num = rep(1:10, each = 6),
                  dbl = as.numeric(rep(1:10, each = 6)),
                  chr = rep(letters[1:10], each = 6),
                  fac = factor(rep(1:10, each = 6)))
  expect_error(term_build(random(~ 1 | cont), d), "is continuous")
  # what stays legal is every ordinary way of writing a grouping variable,
  # integer codes included
  for (g in c("fac", "chr", "num", "dbl")) {
    b <- term_build(random(stats::as.formula(paste("~ 1 |", g))), d)
    expect_identical(term_npar(b), 10L)
  }
  # and the explicit spelling says the values really are labels
  expect_identical(term_npar(term_build(random(~ 1 | factor(cont)), d)), 60L)
})

test_that("a labelled term declares no penalty, and is otherwise unchanged", {
  # The coefficients of a labelled term share a covariance block with those of
  # the other terms carrying the label, and that block may span another
  # equation. The prior over it is the class's, so declaring one here would be
  # a second answer to the same question. An empty list is the right answer.
  d <- dd
  d$u <- as.numeric(1:12)
  a <- term_build(random(~ 1 + x | g), d)
  b <- term_build(random(~ 1 + x | u | g), d)
  expect_identical(term_npar(a), term_npar(b))
  expect_length(term_penalties(a), 1L)
  expect_length(term_penalties(b), 0L)
  expect_null(term_penalty(b))
  expect_identical(term_predict(b, d), term_matrix(b))
  res <- check_term(random(~ 1 + x | u | g), d, verbose = FALSE)
  expect_true(all(res$status == "OK"),
              info = paste(res$check[res$status != "OK"], collapse = ", "))
})


test_that("a labelled term reports its grouping, and reads it back", {
  # a fitting layer needs three things to stack two labelled blocks into one
  # covariance: the grouping expression, its levels and how many columns one
  # level carries. The levels are what tell two groupings apart -- the
  # expression alone cannot, since `id` under two subsets writes the same.
  b <- term_build(random(~ 1 + x | u | g), dd)
  gr <- term_group(b)
  expect_identical(gr$dim, 2L)
  expect_identical(gr$levels, levels(dd$g))
  expect_identical(deparse1(gr$expr), "g")
  # an unlabelled random effect answers too: the grouping is not the label's
  expect_identical(term_group(term_build(random(~ 1 | g), dd))$dim, 1L)
  # and no other kind of term has one
  expect_null(term_group(term_build(linpar(~ x), dd)))
  expect_null(term_group(term_build(s(x), dd)))
})

test_that("what a labelled term refuses, it refuses by name", {
  # `hyper` names a hyperparameter of the term's own penalty, and it has none
  expect_error(term_build(random(~ 1 | u | g, hyper = c(sigma = 1)), dd),
               "belong to the class")
  # the two scopes ask for opposite things at once
  expect_error(term_build(random(~ 1 + x | u | g, correlated = FALSE), dd),
               "two things at once")
  # a univariate prior has no correlations to carry, and margins of different
  # families are a copula rather than an elliptical law
  expect_error(
    term_build(random(~ 1 | u | g,
                      distrib = distributions7::fixed(
                        distributions7::student_t1_distrib(), mu = 0)), dd),
    "copula")
  # a list gives one prior per column, which makes the columns independent
  expect_error(
    term_build(random(~ 1 + x | u | g,
                      distrib = list(.centered(distributions7::gaussian1_distrib(),
                                               mu = 0),
                                     .centered(distributions7::gaussian1_distrib(),
                                               mu = 0))), dd),
    "list")
  # a multivariate one names the class's joint prior and is admitted; its
  # dimension is the CLASS's total and is checked where the class is known
  expect_s7_class(term_build(random(~ 1 + x | u | g,
                                    distrib = .mv_centered(4)), dd),
                  RandomTerm)
  # and none of this reaches an unlabelled term
  expect_no_error(term_build(random(~ 1 | g, hyper = c(sigma = 1)), dd))
  expect_no_error(term_build(random(~ 1 + x | g, correlated = FALSE), dd))
})
