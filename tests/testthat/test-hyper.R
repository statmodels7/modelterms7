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
  expect_error(lasso(~x, lambda = "a"), "must be a number")
  # several values are a GRID and not a held value, so the length is what
  # tells the two states apart and nothing rejects it for being more than one
  expect_silent(lasso(~x, lambda = c(1, 2)))
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
  # the effects are gaussian by default and the hyperparameter is their
  # standard deviation, which is the name the penalty carries
  built <- term_build(random(~ 1 | g, hyper = c(sigma = 6)), dd)
  expect_equal(term_hyper(built)[[1L]], list(sigma = 6))
  # which names there are depends on what the term was given, so the check
  # is at the build, the first point at which the penalty exists
  expect_error(term_build(random(~ 1 | g, hyper = c(lambda = 1)), dd),
               "not a hyperparameter")
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

test_that("a term says how fine its own grid is, per hyperparameter", {
  # HOW FINELY is the term's answer for the same reason as WHICH: a block of
  # four columns and one of four hundred want different grids, and the
  # criterion is put to every hyperparameter of the model at once, the
  # smooth ones included, so it cannot know which it is looking at.
  #
  # The default is a NUMBER ON THE SIGNATURE and not a criterion's that
  # nothing prints, so the term always answers.
  expect_equal(term_grid(lasso(~x))[[1L]], list(lambda = 25L))
  expect_equal(term_grid(enet(~x))[[1L]], list(lambda = 25L, alpha = 5L))
  expect_equal(term_grid(lasso(~x, n_lambda = 50))[[1L]], list(lambda = 50L))
  expect_equal(term_grid(enet(~x, n_lambda = 40, n_alpha = 12))[[1L]],
               list(lambda = 40L, alpha = 12L))
  expect_equal(term_grid(scad(~x, n_a = 6))[[1L]], list(lambda = 25L, a = 6L))
  expect_equal(term_grid(mcp(~x, n_gamma = 7))[[1L]],
               list(lambda = 25L, gamma = 7L))
  # the two axes differ because what they span differs: lambda descends the
  # size of the kink over four decades, a shape spans its useful range
  expect_gt(eval(formals(enet)$n_lambda), eval(formals(enet)$n_alpha))

  # a grid on a hyperparameter the penalty does not carry, and a grid that
  # is not a grid
  expect_error(lasso(~x, n_alpha = 5), "no argument 'n_alpha'")
  expect_error(mcp(~x, n_gamma = 1), "at least 2")
  expect_error(mcp(~x, n_gamma = 2.5), "whole number")

  # and it travels with the entry, as the held values do
  dd <- data.frame(x = stats::rnorm(20), z = stats::rnorm(20))
  b <- term_build(enet(~ x + z, n_lambda = 9), dd)
  expect_equal(term_penalties(b)[[1L]]$n_values,
               list(lambda = 9L, alpha = 5L))
})

test_that("a term says how far down its own path reaches", {
  expect_equal(term_path_min(lasso(~x))[[1L]], 1e-4)
  expect_equal(term_path_min(lasso(~x, min_ratio = 1e-6))[[1L]], 1e-6)
  # one number per term and not one per hyperparameter: only the sweep by
  # kink size uses it
  expect_error(lasso(~x, min_ratio = 1), "in (0, 1)", fixed = TRUE)
  expect_error(lasso(~x, min_ratio = 0), "in (0, 1)", fixed = TRUE)
  expect_error(lasso(~x, min_ratio = c(0.1, 0.2)), "single number")
  # and it travels with the entry
  dd <- data.frame(x = stats::rnorm(20), z = stats::rnorm(20))
  b <- term_build(scad(~ x + z, min_ratio = 1e-3), dd)
  expect_equal(term_penalties(b)[[1L]]$min_ratio, 1e-3)
})

test_that("several values are a grid the path visits, not a held value", {
  # one argument, three states, settled per hyperparameter: NULL builds the
  # grid, one number holds, several ARE the grid
  expect_identical(term_values(lasso(~x)), list())
  expect_identical(term_hyper(lasso(~x, lambda = 3))[[1L]], list(lambda = 3))
  expect_identical(term_values(lasso(~x, lambda = 3)), list())

  v <- c(0.1, 1, 10)
  expect_identical(term_values(lasso(~x, lambda = v))[[1L]], list(lambda = v))
  # and a written-out grid is NOT held: what the caller fixed is where to
  # look, not the answer
  expect_identical(term_hyper(lasso(~x, lambda = v)), list())

  # independently per hyperparameter, which is the case the contract exists
  # for: a grid on one and a built grid on the other
  e <- enet(~x, lambda = v)
  expect_identical(term_values(e)[[1L]], list(lambda = v))
  expect_identical(term_hyper(e), list())
  e2 <- enet(~x, lambda = v, alpha = 0.5)
  expect_identical(term_values(e2)[[1L]], list(lambda = v))
  expect_identical(term_hyper(e2)[[1L]], list(alpha = 0.5))

  # sorted and deduplicated, a path being walked in one direction
  expect_identical(term_values(mcp(~x, gamma = c(4, 2, 4, 3)))[[1L]],
                   list(gamma = c(2, 3, 4)))

  # every value is checked against the bounds, as a held one is
  expect_error(enet(~x, alpha = c(0.2, 1.5)), "strictly inside")
  expect_error(scad(~x, a = c(3, 1)), "strictly inside")
  expect_error(lasso(~x, lambda = c(1, NA)), "several numbers")
  expect_error(lasso(~x, lambda = character(0)), "several numbers")

  # a penalty with no kink has no path to visit them on, and says so rather
  # than taking the vector and doing nothing with it. The question is put to
  # the PENALTY at a probe value, so a random effect under a Gaussian prior
  # is covered by the same line as a ridge.
  expect_error(ridge(~x, lambda = c(1, 2)), "no path")
  dg <- data.frame(y = stats::rnorm(30), g = factor(rep(1:6, each = 5)))
  expect_error(term_build(random(~ 1 | g, hyper = list(sigma = c(1, 2))), dg),
               "no path")

  # and it travels with the entry, as the held values and the grid size do
  dd <- data.frame(x = stats::rnorm(20), z = stats::rnorm(20))
  b <- term_build(scad(~ x + z, lambda = c(0.5, 2)), dd)
  expect_equal(term_penalties(b)[[1L]]$values, list(lambda = c(0.5, 2)))
  expect_identical(term_penalties(b)[[1L]]$fixed, list())
})

test_that("a term says how its own hyperparameters are combined", {
  # HOW they are combined with EACH OTHER is not a property any one of them
  # has, so it is one word per term. It belongs to the term because a penalty
  # with a kink is fitted by a scheme of its own; a criterion is asked of
  # every hyperparameter of the model, the smooth ones included, and most of
  # those are read at the mode rather than swept.
  expect_equal(term_search(enet(~x))[[1L]], "grid")
  expect_equal(term_search(enet(~x, search = "cyclic"))[[1L]], "cyclic")
  expect_equal(term_search(scad(~x, search = "grid"))[[1L]], "grid")
  expect_equal(term_search(mcp(~x, search = "cyclic"))[[1L]], "cyclic")

  expect_error(enet(~x, search = "product"), "must be")
  expect_error(enet(~x, search = c("grid", "cyclic")), "must be")
  expect_error(enet(~x, search = 1), "must be")

  # the two with one hyperparameter or none do not carry it: there is
  # nothing to combine, and the argument is reported rather than swallowed
  expect_error(lasso(~x, search = "grid"), "no argument 'search'")
  expect_error(ridge(~x, search = "grid"), "no argument 'search'")

  # and it travels with the entry, as the grid size and the depth do
  dd <- data.frame(x = stats::rnorm(20), z = stats::rnorm(20))
  b <- term_build(enet(~ x + z, search = "cyclic"), dd)
  expect_equal(term_penalties(b)[[1L]]$search, "cyclic")
})

test_that("a term does not take an argument it would not read", {
  # a ridge is twice differentiable everywhere, so its lambda is estimated at
  # the mode and there is no path: a grid size and a depth describe nothing
  # here, and they were accepted only so the five constructors read alike
  expect_error(ridge(~x, n_lambda = 10), "no argument 'n_lambda'")
  expect_error(ridge(~x, min_ratio = 1e-3), "no argument 'min_ratio'")
  # while the four that do have a path keep both
  expect_equal(term_grid(lasso(~x, n_lambda = 10))[[1L]], list(lambda = 10L))
  expect_equal(term_path_min(lasso(~x, min_ratio = 1e-3))[[1L]], 1e-3)
})
