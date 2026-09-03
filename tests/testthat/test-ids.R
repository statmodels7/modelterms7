test_that("a single label needs no name where the penalty carries one", {
  # the answer is keyed by entry name, and a term with one penalty over
  # the whole of itself names that entry "", as term_hyper() does
  one <- function(v) stats::setNames(list(v), "")
  expect_identical(term_ids(ridge(~ x, id = "L")), one(c(lambda = "L")))
  expect_identical(term_ids(lasso(~ x, id = "L")), one(c(lambda = "L")))
  expect_identical(term_ids(s(x, id = "L")), one(c(lambda = "L")))
  # an isotropic tensor product carries lambda alone, so the short form works
  expect_identical(term_ids(te(x, z, anisotropic = FALSE, id = "L")),
                   one(c(lambda = "L")))
})

test_that("a single label is refused where the penalty carries several", {
  # what it would mean depends on which was meant, and the reply says which
  # there are rather than picking one by position
  expect_error(enet(~ x, id = "L"), "carries 2: lambda, alpha")
  expect_error(scad(~ x, id = "L"), "carries 2: lambda, a")
  expect_error(mcp(~ x, id = "L"), "carries 2: lambda, gamma")
  expect_error(te(x, z, id = "L"), "lambda1, lambda2")
})

test_that("named labels reach the hyperparameters they name", {
  one <- function(v) stats::setNames(list(v), "")
  expect_identical(term_ids(enet(~ x, id = c(alpha = "A"))),
                   one(c(alpha = "A")))
  expect_identical(term_ids(enet(~ x, id = c(lambda = "L", alpha = "A"))),
                   one(c(lambda = "L", alpha = "A")))
  expect_identical(term_ids(te(x, z, id = c(lambda1 = "A"))),
                   one(c(lambda1 = "A")))
})

test_that("a name the penalty does not carry is reported where it is written", {
  expect_error(mcp(~ x, id = c(a = "A")), "is not a hyperparameter")
  expect_error(s(x, id = c(sigma = "A")), "It carries: lambda")
})

test_that("id must be a non-empty string", {
  expect_error(ridge(~ x, id = 1), "must be a non-empty string")
  expect_error(ridge(~ x, id = ""), "must be a non-empty string")
  expect_error(ridge(~ x, id = NA_character_), "must be a non-empty string")
})

test_that("a term sharing nothing answers with an empty list", {
  expect_identical(term_ids(ridge(~ x)), list())
  expect_identical(term_ids(s(x)), list())
  expect_identical(term_ids(linpar(~ x)), list())
  expect_identical(term_ids(enet(~ x)), list())
})

test_that("sharing and holding are independent", {
  # the two arguments say different things about the same hyperparameter and
  # neither reads the other
  tm <- s(x, lambda = 2, id = "L")
  expect_identical(term_hyper(tm),
                   stats::setNames(list(list(lambda = 2)), ""))
  expect_identical(term_ids(tm),
                   stats::setNames(list(c(lambda = "L")), ""))
})

test_that("a random effect's labels are checked where its penalty exists", {
  set.seed(11)
  n <- 120
  d <- data.frame(y = rnorm(n), x = rnorm(n),
                  g = factor(sample(letters[1:6], n, TRUE)))
  # the constructor cannot check: the names are the effects' distribution's,
  # qualified by the entry, and neither is known before the group is read
  expect_silent(random(~ 1 | g, id = c(zeta = "R")))
  expect_error(term_build(random(~ 1 | g, id = c(zeta = "R")), d),
               "is not a hyperparameter")

  b <- term_build(random(~ 1 | g, id = "R"), d)
  expect_identical(term_ids(b), stats::setNames(list(c(sigma = "R")), ""))
  expect_identical(term_penalties(b)[[1L]]$ids, c(sigma = "R"))

  # one entry per within-group column when they do not correlate, so the
  # label carries the column as `hyper` does
  b2 <- term_build(random(~ 1 + x | g, correlated = FALSE,
                          id = c("sigma.(Intercept)" = "R")), d)
  ids <- term_ids(b2)
  expect_identical(names(ids), "(Intercept)")
  expect_identical(ids[[1L]], c(sigma = "R"))
})

test_that("a label reaches a term written inside a subformula", {
  set.seed(12)
  n <- 200
  d <- data.frame(x = rep(seq(0, 10, length.out = 25), 8),
                  g = factor(rep(seq_len(8), each = 25)))
  d$y <- 1 + 0.3 * d$x + 1.2 * pmax(d$x - 5, 0) + rnorm(n, sd = 0.5)
  b <- term_build(seg(x, beta ~ 0 + ridge(~ g, id = "S"),
                      gamma1 ~ 0 + ridge(~ g, id = "S")), d)
  got <- term_ids(b)
  expect_length(got, 2L)
  expect_true(all(vapply(got, function(v) identical(v, c(lambda = "S")),
                         logical(1))))
  # and the entries carry it, which is what the fitting layer reads
  for (e in term_penalties(b)) expect_identical(e$ids, c(lambda = "S"))
})

test_that("a term is unaffected by an id it was not given", {
  # the negative control: every quantity a term reports is the same with the
  # argument absent as it was before the argument existed
  set.seed(13)
  n <- 150
  d <- data.frame(y = rnorm(n), x = rnorm(n),
                  g = factor(sample(letters[1:5], n, TRUE)))
  a <- term_build(s(x, k = 6), d)
  b <- term_build(s(x, k = 6, id = "L"), d)
  expect_identical(term_matrix(a), term_matrix(b))
  expect_identical(term_npar(a), term_npar(b))
  expect_identical(term_hyper(a), term_hyper(b))
  expect_identical(term_penalty(a)@params, term_penalty(b)@params)
  expect_identical(term_ids(a), list())
})

test_that("a labelled random term has no hyperparameter of its own to share", {
  set.seed(14)
  n <- 150
  d <- data.frame(y = rnorm(n), x = rnorm(n),
                  g = factor(sample(letters[1:5], n, TRUE)))
  # the two bars already share the whole block, and the prior over it belongs
  # to the class -- so there is nothing here for id to name
  expect_error(term_build(random(~ 1 | u | g, id = "R"), d),
               "the covariance block it shares is")
  # and the same term without a label takes it
  expect_silent(term_build(random(~ 1 | g, id = "R"), d))
})
