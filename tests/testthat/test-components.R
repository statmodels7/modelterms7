test_that("a term whose columns do not divide answers with nothing", {
  dd <- data.frame(x = seq(0.2, 3, length.out = 40),
                   g = factor(rep(c("a", "b"), 20)))
  dd$y <- 1 + 0.5 * dd$x
  for (tm in list(linpar(~ x), s(x, k = 6), random(~ 1 | g), ridge(~ x + g),
                  lasso(~ x + g))) {
    expect_identical(term_components(term_build(tm, dd)), list())
  }
})

test_that("the components partition the columns of a term written in its own parameters", {
  dd <- data.frame(x = seq(0.2, 3, length.out = 40),
                   g = factor(rep(c("a", "b"), 20)))
  dd$y <- 2 * exp(-1.3 * dd$x)
  built <- list(
    nl = term_build(nl(~ a * exp(-r * x), start = list(a = 2, r = 1.3)), dd),
    nl_dev = term_build(nl(~ a * exp(-r * x), a ~ 0 + g,
                           start = list(r = 1.3)), dd),
    seg = term_build(seg(x, psi = 1.5), dd),
    jump = term_build(jump(x, psi = 1.5), dd),
    jseg = term_build(jseg(x, psi = 1.5), dd))
  for (nm in names(built)) {
    cp <- term_components(built[[nm]])
    expect_gt(length(cp), 0)
    idx <- unlist(lapply(cp, `[[`, "index"), use.names = FALSE)
    # every column belongs to exactly one own parameter: a partition, not a
    # cover, so a consumer that reports the components reports the whole
    # block and reports nothing twice
    expect_identical(sort(as.integer(idx)),
                     seq_len(term_npar(built[[nm]])))
    expect_identical(names(cp),
                     unname(vapply(cp, `[[`, character(1), "name")))
  }
})

test_that("sub_index divides a developed parameter among its sub-terms", {
  dd <- data.frame(x = seq(0.2, 3, length.out = 60),
                   id = factor(rep(letters[1:6], each = 10)),
                   g = factor(rep(c("a", "b"), 30)))
  dd$y <- 2 * exp(-1.3 * dd$x)
  built <- list(
    term_build(nl(~ a * exp(-r * x), a ~ 0 + g, start = list(r = 1.3)), dd),
    term_build(seg(x, psi ~ random(~ 1 | id), psi = 1.5), dd),
    term_build(jump(x, psi ~ 0 + id, psi = 1.5), dd))
  for (b in built) {
    nms <- term_coef_names(b)
    for (cp in term_components(b)) {
      if (!length(cp$subs)) {
        expect_identical(cp$sub_index, list())
        next
      }
      expect_identical(length(cp$sub_index), length(cp$subs))
      expect_identical(unlist(cp$sub_index, use.names = FALSE),
                       as.integer(cp$index))
      for (i in seq_along(cp$subs)) {
        si <- cp$sub_index[[i]]
        expect_identical(length(si), as.integer(term_npar(cp$subs[[i]])))
        # THE STRUCTURAL CHECK, and the one a count alone would not make:
        # the term's own coefficient names at those columns end in the
        # sub-term's, so the columns really are that sub-term's and not the
        # next one's shifted by an equal count
        expect_true(all(endsWith(nms[si], term_coef_names(cp$subs[[i]]))))
      }
    }
  }
})

test_that("a parameter developed by two kinds at once reports both", {
  dd <- data.frame(x = seq(0.2, 3, length.out = 60),
                   id = factor(rep(letters[1:6], each = 10)))
  dd$y <- 2 + 0.5 * dd$x
  b <- term_build(seg(x, psi ~ random(~ 1 | id), psi = 1.5), dd)
  cp <- term_components(b)[["psi1"]]
  expect_identical(length(cp$subs), 2L)
  expect_true(S7::S7_inherits(cp$subs[[1L]], LinparTerm))
  expect_true(S7::S7_inherits(cp$subs[[2L]], RandomTerm))
})

test_that("a structural term divides its parameter vector", {
  set.seed(15)
  n <- 60
  dd <- data.frame(t = seq_len(n), y = rnorm(n),
                   z = as.numeric(scale(runif(n))),
                   g = factor(rep(c("a", "b"), each = n / 2)))
  b <- term_build(gas(p = 1, q = 1, omega ~ random(~ 1 | g), by = g,
                      time = t), dd)
  cp <- term_components(b)
  expect_identical(names(cp), c("omega", "alpha1", "pacf1"))
  # a structural term has no design columns: what is divided is the vector
  # term_params() reports, which its state and its variance matrix share
  idx <- unlist(lapply(cp, `[[`, "index"), use.names = FALSE)
  expect_identical(sort(as.integer(idx)),
                   seq_along(term_params(b)))
  expect_identical(length(cp$omega$subs), 2L)
  expect_identical(cp$alpha1$subs, list())
  expect_identical(unlist(cp$omega$sub_index, use.names = FALSE),
                   as.integer(cp$omega$index))
  # and the positions really are that sub-term's: the term's own parameter
  # names at those positions end in the sub-term's coefficient names
  nms <- term_params(b)
  for (i in seq_along(cp$omega$subs)) {
    expect_true(all(endsWith(nms[cp$omega$sub_index[[i]]],
                             term_coef_names(cp$omega$subs[[i]]))))
  }
})

test_that("a scalar structural term divides into one position each", {
  set.seed(15)
  n <- 60
  dd <- data.frame(t = seq_len(n), y = rnorm(n))
  b <- term_build(gas(p = 1, q = 1, time = t), dd)
  cp <- term_components(b)
  expect_identical(vapply(cp, function(z) length(z$index), integer(1)),
                   c(omega = 1L, alpha1 = 1L, pacf1 = 1L))
  expect_true(all(vapply(cp, function(z) length(z$subs) == 0L, logical(1))))
})
