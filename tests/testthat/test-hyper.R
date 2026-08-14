test_that("a constructor holds a hyperparameter and leaves the rest free", {
  # WHICH hyperparameters are estimated is a property of the TERM, since the
  # term is where the penalty is named. NULL, the default, means estimated.
  expect_identical(term_hyper(lasso(~x)), list())
  expect_equal(term_hyper(lasso(~x, lambda = 3)),
               stats::setNames(list(list(lambda = 3)), ""))
  # part of them: alpha held, lambda still to be estimated
  h <- term_hyper(enet(~x, alpha = 0.3))
  expect_identical(names(h[[1L]]), "alpha")
  expect_equal(h[[1L]]$alpha, 0.3)
  expect_equal(term_hyper(scad(~x, lambda = 2, a = 3.7))[[1L]],
               list(lambda = 2, a = 3.7))
  expect_equal(term_hyper(mcp(~x, gamma = 3))[[1L]], list(gamma = 3))
  expect_equal(term_hyper(ridge(~x, lambda = 0.5))[[1L]], list(lambda = 0.5))
})

test_that("the arguments carry the penalty's own names", {
  # a Gaussian prior is parametrized by its SCALE, so `ridge(lambda = 2)` is
  # a reasonable thing to write and reaches nothing. R reports an unused
  # argument, which says it was not read and not what to write instead.
  expect_error(ridge(~x, sigma = 2), "no argument 'sigma'")
  expect_error(ridge(~x, sigma = 2), "lambda")
  expect_error(scad(~x, gamma = 3), "lambda, a")
  expect_error(mcp(~x, a = 3), "lambda, gamma")
  expect_error(lasso(~x, alpha = 0.5), "no argument 'alpha'")
})

test_that("a held value must lie strictly inside the penalty's bounds", {
  expect_error(enet(~x, alpha = 1.5), "strictly inside")
  expect_error(enet(~x, alpha = 0), "strictly inside")
  expect_error(scad(~x, a = 2), "strictly inside")
  expect_error(mcp(~x, gamma = 1), "strictly inside")
  expect_error(lasso(~x, lambda = 0), "strictly inside")
  expect_error(lasso(~x, lambda = c(1, 2)), "single finite number")
  expect_error(lasso(~x, lambda = "a"), "single finite number")
})

test_that("a smooth holds its smoothing parameter, one per margin", {
  expect_identical(term_hyper(s(x)), list())
  expect_equal(term_hyper(s(x, lambda = 2))[[1L]], list(lambda = 2))
  expect_equal(term_hyper(te(x, z, lambda = c(1, 5)))[[1L]],
               list(lambda1 = 1, lambda2 = 5))
  # named, so one margin is held and the other estimated
  expect_equal(term_hyper(te(x, z, lambda = c(lambda2 = 5)))[[1L]],
               list(lambda2 = 5))
  # an isotropic tensor product has ONE, so one number is right there
  expect_equal(term_hyper(te(x, z, anisotropic = FALSE, lambda = 3))[[1L]],
               list(lambda = 3))
  # and one number for an anisotropic one is a different model, not a
  # shorthand for this one
  expect_error(te(x, z, lambda = 2), "one per margin")
  expect_error(s(x, lambda = 0), "strictly positive")
  expect_error(te(x, z, lambda = c(lambda3 = 1)), "no smoothing parameter")
})

test_that("a random effect is checked against the penalty it builds", {
  dd <- data.frame(y = stats::rnorm(20), x = stats::rnorm(20),
                   g = factor(rep(1:4, 5)))
  built <- term_build(random(~ 1 | g, hyper = c(lambda = 6)), dd)
  expect_equal(term_hyper(built)[[1L]], list(lambda = 6))
  # which names there are depends on what the term was given, so the check
  # is at the build, the first point at which the penalty exists
  expect_error(term_build(random(~ 1 | g, hyper = c(sigma = 1)), dd),
               "no hyperparameter 'sigma'")
  expect_error(random(~ 1 | g, hyper = 0.4), "must be named")
})

test_that("a held value travels with the entry, through a structural term", {
  # a term that copies its sub-terms' entries propagates what they hold
  # without knowing hyperparameters exist
  set.seed(9)
  dd <- data.frame(y = stats::rnorm(60), x = stats::rnorm(60),
                   id = factor(rep(1:6, 10)))
  b <- term_build(nl(~ a * exp(-r * x), a ~ 0 + ridge(~ id, lambda = 0.7),
                     start = list(a = 1, r = 0.5)), dd)
  h <- term_hyper(b)
  expect_true(length(h) >= 1L)
  expect_equal(h[[1L]]$lambda, 0.7)
  ent <- term_penalties(b)
  expect_equal(ent[[1L]]$fixed$lambda, 0.7)
})
