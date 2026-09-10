# A smooth's penalty is a choice, and there are two of them: WHICH penalty
# covers the coordinates, and whether the levels of a factor `by` are
# smoothed together or apart. Both are answered through term_penalties(),
# and the default answer is the one that shipped.

set.seed(11)
n <- 240
dd <- data.frame(x = sort(runif(n)),
                 g = factor(rep(c("a", "b", "c"), length.out = n)),
                 w = rnorm(n))
dd$y <- sin(2 * pi * dd$x) + rnorm(n, sd = 0.2)

lp <- penalties7::lasso_penalty

test_that("the default declares one penalty over the whole block", {
  b <- term_build(s(x, bspline_smooth(k = 10)), dd)
  ent <- term_penalties(b)
  expect_length(ent, 1L)
  expect_identical(ent[[1L]]$name, "")
  expect_identical(ent[[1L]]$index, seq_len(term_npar(b)))
  # and it is the term's own penalty, not a second one built beside it
  expect_identical(ent[[1L]]$penalty, term_penalty(b))
  # which is what the base method answers: the method delegates rather than
  # repeating it, so the two cannot drift
  expect_identical(ent, S7::method(term_penalties, model_term)(b))
})

test_that("a factory covers the penalized coordinates and not the free ones", {
  # A ROUGHNESS MATRIX CARRIES A ZERO ROW for the linear column and leaves
  # it free by its own arithmetic; a separable penalty has no such row and
  # would shrink it, so the free columns are outside the entry.
  b <- term_build(s(x, bspline_smooth(k = 10, penalty = lp)), dd)
  ent <- term_penalties(b)
  expect_length(ent, 1L)
  expect_identical(ent[[1L]]$index, 2:9)
  expect_identical(as.integer(ent[[1L]]$penalty@n_coef), 8L)
  expect_true(penalties7::has_prox(ent[[1L]]$penalty))
  # the block itself is the one the roughness matrix would have produced:
  # the factory penalizes the coordinates, it does not build them
  ref <- term_build(s(x, bspline_smooth(k = 10)), dd)
  expect_identical(term_matrix(b), term_matrix(ref))
  expect_identical(term_coef_names(b), term_coef_names(ref))
  # there is no single penalty over the whole block, and the property says
  # so rather than carrying one of the entries
  expect_null(term_penalty(b))

  # dropping the null space leaves nothing free, so the factory covers all
  d0 <- term_build(s(x, bspline_smooth(k = 10, null_space = "drop",
                                       penalty = lp)), dd)
  expect_identical(term_penalties(d0)[[1L]]$index, seq_len(term_npar(d0)))
})

test_that("by_hyper = 'level' declares one penalty per level", {
  b <- term_build(s(x, bspline_smooth(k = 10), by = g, by_hyper = "level"), dd)
  ent <- term_penalties(b)
  expect_length(ent, 3L)
  expect_identical(vapply(ent, function(e) e$name, character(1)),
                   levels(dd$g))
  # group-major: each level's columns are adjacent, and the levels partition
  # the block
  expect_identical(ent[[1L]]$index, 1:9)
  expect_identical(ent[[2L]]$index, 10:18)
  expect_identical(ent[[3L]]$index, 19:27)
  expect_identical(sort(unlist(lapply(ent, function(e) e$index))),
                   seq_len(term_npar(b)))
  # one smoothing parameter each, and each is the SAME matrix: what differs
  # between the levels is the value estimated, not the penalty
  expect_true(all(vapply(ent, function(e) e$penalty@params, character(1)) ==
                    "lambda"))
  expect_identical(ent[[1L]]$penalty@P, ent[[2L]]$penalty@P)

  # shared is the default and is one entry over all of them
  sh <- term_build(s(x, bspline_smooth(k = 10), by = g), dd)
  expect_length(term_penalties(sh), 1L)
  expect_identical(term_matrix(b), term_matrix(sh))
})

test_that("a factory and a per-level by compose", {
  b <- term_build(s(x, bspline_smooth(k = 10, penalty = lp), by = g,
                    by_hyper = "level"), dd)
  ent <- term_penalties(b)
  expect_length(ent, 3L)
  # each level's own penalized coordinates, the free column of each left out
  expect_identical(ent[[1L]]$index, 2:9)
  expect_identical(ent[[2L]]$index, 11:18)
  expect_identical(ent[[3L]]$index, 20:27)
  expect_true(all(vapply(ent, function(e) as.integer(e$penalty@n_coef),
                         integer(1)) == 8L))

  # shared over the same block is ONE penalty over every level's penalized
  # coordinates: repeating a separable penalty blockwise is the identity
  # operation, so this is the same penalty under one smoothing parameter
  sh <- term_build(s(x, bspline_smooth(k = 10, penalty = lp), by = g), dd)
  e1 <- term_penalties(sh)
  expect_length(e1, 1L)
  expect_identical(e1[[1L]]$index,
                   sort(unlist(lapply(ent, function(e) e$index))))
  expect_identical(as.integer(e1[[1L]]$penalty@n_coef), 24L)
})

test_that("the hyperparameters are the penalty's own, checked at the build", {
  # WHICH THERE ARE is not known when the term is written: a factory names
  # them at a coefficient count the data settle. The shape is checked at the
  # constructor and the names here.
  expect_silent(s(x, bspline_smooth(penalty = penalties7::heavy_penalty),
                  hyper = c(nu = 3)))
  b <- term_build(s(x, bspline_smooth(k = 10,
                                      penalty = penalties7::heavy_penalty),
                    hyper = c(nu = 3)), dd)
  expect_identical(term_penalties(b)[[1L]]$penalty@params, c("sigma", "nu"))
  expect_identical(term_penalties(b)[[1L]]$fixed, list(nu = 3))
  # a name that penalty has not is an error naming what it carries
  expect_error(term_build(s(x, bspline_smooth(k = 10,
                                              penalty = lp),
                            hyper = c(nu = 3)), dd),
               "not a hyperparameter")
  # under by_hyper = "level" a held value is qualified by the level
  bl <- term_build(s(x, bspline_smooth(k = 10), by = g, by_hyper = "level",
                     hyper = c(lambda.b = 4)), dd)
  ent <- term_penalties(bl)
  expect_identical(ent[[1L]]$fixed, list())
  expect_identical(ent[[2L]]$fixed, list(lambda = 4))
  expect_identical(ent[[3L]]$fixed, list())
  # and an unqualified one is an error rather than a value recycled over
  # every level
  expect_error(term_build(s(x, bspline_smooth(k = 10), by = g,
                            by_hyper = "level", hyper = c(lambda = 4)), dd),
               "not a hyperparameter")
})

test_that("te() splits its smoothing parameters by level too", {
  dt <- data.frame(x = dd$x, z = stats::runif(nrow(dd)), g = dd$g)
  # ANISOTROPIC is a SET of smoothing parameters, one per margin, so per
  # level is one set each rather than one parameter each
  a <- term_build(te(x, z, smooths = bspline_smooth(k = 4), by = g,
                     by_hyper = "level"), dt)
  ea <- term_penalties(a)
  expect_length(ea, 3L)
  expect_identical(vapply(ea, function(e) e$name, character(1)), levels(dt$g))
  expect_identical(ea[[1L]]$penalty@params, c("lambda1", "lambda2"))
  # the levels partition the block, in level-major order
  expect_identical(sort(unlist(lapply(ea, function(e) e$index))),
                   seq_len(term_npar(a)))
  expect_identical(max(ea[[1L]]$index) + 1L, min(ea[[2L]]$index))

  # isotropic is ONE parameter, so per level is one each
  i <- term_build(te(x, z, smooths = bspline_smooth(k = 4), by = g,
                     by_hyper = "level", anisotropic = FALSE), dt)
  expect_identical(term_penalties(i)[[1L]]$penalty@params, "lambda")

  # and shared is one entry over every level, as it was
  s0 <- term_build(te(x, z, smooths = bspline_smooth(k = 4), by = g), dt)
  expect_length(term_penalties(s0), 1L)
  expect_identical(term_matrix(a), term_matrix(s0))
})

test_that("an id reaches one level's smoothing parameter", {
  # SHARING AND SPLITTING ARE THE SAME MACHINERY read in two directions:
  # `by_hyper = "level"` gives the levels a smoothing parameter each, and
  # `id` then shares ONE of them with another term. Neither needed writing
  # for the other -- the entries carry qualified names, so `check_ids()`
  # sees exactly the names `hyper` does.
  b <- term_build(s(x, bspline_smooth(k = 6), by = g, by_hyper = "level",
                    id = c(lambda.a = "sp")), dd)
  ent <- term_penalties(b)
  expect_identical(ent[[1L]]$ids, c(lambda = "sp"))
  # empty, and named-empty at that: subsetting a named vector by an empty
  # intersection keeps the names attribute, which identical() distinguishes
  # from character(0). What the test is about is that nothing is shared.
  expect_length(ent[[2L]]$ids, 0L)
  expect_length(ent[[3L]]$ids, 0L)
  expect_identical(term_ids(b), list(a = c(lambda = "sp")))
  # a name no level answers for lists the ones that exist
  expect_error(term_build(s(x, bspline_smooth(k = 6), by = g,
                            by_hyper = "level", id = c(lambda.z = "sp")), dd),
               "lambda.a, lambda.b, lambda.c", fixed = TRUE)
})

test_that("the combinations that contradict each other are rejected", {
  # a numeric `by` has no levels to give a smoothing parameter each
  expect_error(term_build(s(x, bspline_smooth(k = 10), by = w,
                            by_hyper = "level"), dd),
               "no factor 'by'")
  expect_error(term_build(s(x, bspline_smooth(k = 10), by_hyper = "level"), dd),
               "no factor 'by'")
  # a tensor product's penalty is a sum over the margins, and a factory does
  # not say how it composes with that sum
  expect_error(te(x, w, smooths = bspline_smooth(k = 4, penalty = lp)),
               "not read by te")
  # an unknown value is rejected by name
  expect_error(s(x, by_hyper = "nope"), "arg")
})

test_that("a factory returning the wrong thing is caught where it is called", {
  expect_error(term_build(s(x, bspline_smooth(k = 10,
                                              penalty = function(n_coef) 1)),
                          dd),
               "must give a penalties7 penalty")
  # a penalty of the wrong width would be evaluated at a coefficient vector
  # of another length and R would recycle it in silence
  expect_error(term_build(
    s(x, bspline_smooth(k = 10,
                        penalty = function(n_coef) lp(n_coef = 3L))), dd),
    "covers 3 coefficients")
})
