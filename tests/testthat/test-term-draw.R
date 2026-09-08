# Drawing a structural term's own parameters. What is asserted is the shape,
# the scale it is on and the two properties the charts are there to guarantee:
# whatever comes out, a loading is positive and a persistence is stationary.

test_that("a draw has term_start's shape and spreads around it", {
  set.seed(1)
  g <- gas(p = 2, q = 2)
  z <- term_draw(g)
  expect_identical(names(z), term_params(g))
  expect_length(z, term_npar(g))
  expect_true(all(is.finite(z)))

  # it is centred on the start rather than on zero, which for a loading is
  # the whole point: zero on a log chart is a loading of one
  m <- rowMeans(replicate(400, term_draw(g)))
  expect_equal(unname(m), unname(term_start(g)), tolerance = 0.1)
  expect_lt(term_start(g)[["alpha1"]], -2)
})

test_that("the width is half of sd", {
  set.seed(2)
  g <- gas(p = 1, q = 1)
  for (s in c(0.5, 1, 3)) {
    v <- replicate(2000, term_draw(g, sd = s)[["omega"]])
    expect_equal(stats::sd(v), s / 2, tolerance = 0.06)
  }
  # a width of zero is the start itself
  expect_identical(term_draw(g, sd = 0), term_start(g))
})

test_that("a draw is admissible on every chart it rides", {
  set.seed(3)
  g <- gas(p = 2, q = 2)
  lk <- term_links(g)
  z <- replicate(200, term_draw(g, sd = 4))
  for (j in seq_len(ncol(z))) {
    psi <- vapply(term_params(g),
                  function(q) linkfunctions7::linkinv(lk[[q]], z[q, j]),
                  numeric(1))
    # the loadings are positive and the persistences stationary, at any
    # coordinate whatever, which is what the charts are for
    expect_true(all(psi[c("alpha1", "alpha2")] > 0))
    expect_true(all(abs(psi[c("pacf1", "pacf2")]) < 1))
  }
})

test_that("a regime's levels come out ordered, by the chart alone", {
  set.seed(4)
  r <- regime(k = 3)
  lk <- term_links(r)
  nm <- term_params(r)
  for (i in seq_len(200)) {
    z <- term_draw(r, sd = 3)
    psi <- vapply(nm, function(q) linkfunctions7::linkinv(lk[[q]], z[[q]]),
                  numeric(1))
    gaps <- psi[grepl("gap", nm, fixed = TRUE)]
    if (length(gaps)) expect_true(all(gaps > 0))
  }
})

test_that("a developed term draws every coordinate, not only the scalars", {
  set.seed(5)
  dd <- data.frame(id = factor(rep(1:6, each = 8)),
                   t = rep(1:8, 6), y = stats::rnorm(48))
  b <- term_build(gas(p = 1, q = 1, omega ~ 1 + random(~1 | id),
                      by = id, time = t), dd)
  z <- term_draw(b)
  expect_identical(names(z), term_params(b))
  expect_true(all(is.finite(z)))
  # the development's own coordinates are drawn rather than left where
  # term_start() puts them, which for a deviation is zero
  dev <- grep("random", names(z), fixed = TRUE)
  expect_length(dev, 6L)
  expect_true(all(term_start(b)[dev] == 0))
  expect_true(all(z[dev] != 0))
})

test_that("an ordinary term has no parameters of its own to draw", {
  expect_error(term_draw(linpar(~1)), "term_draw")
  expect_error(term_draw(s(x, k = 5)), "term_draw")
})

test_that("the draw reads the caller's stream and sets no seed", {
  g <- gas(p = 1, q = 1)
  set.seed(6)
  a <- term_draw(g)
  set.seed(6)
  expect_identical(term_draw(g), a)
  set.seed(6)
  expect_false(isTRUE(all.equal(term_draw(g), term_draw(g))))
})

test_that("'sd' is checked where it is given", {
  g <- gas(p = 1, q = 1)
  expect_error(term_draw(g, sd = -1), "non-negative")
  expect_error(term_draw(g, sd = c(1, 2)), "one finite")
  expect_error(term_draw(g, sd = NA_real_), "one finite")
})
