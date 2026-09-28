# term_starts(): one start by default, n_start of them for regime(), the first
# always term_start()'s.

test_that("the base method returns term_start() alone", {
  g <- gas(p = 1, q = 1)
  st <- term_starts(g)
  expect_length(st, 1L)
  expect_identical(st[[1]], term_start(g))
})

test_that("regime() returns one start unless asked for more", {
  y <- c(1, 2, 3, 10, 11, 12)
  expect_length(term_starts(regime(k = 2), target = y), 1L)
  expect_identical(term_starts(regime(k = 2), target = y)[[1]],
                   term_start(regime(k = 2), target = y))
})

test_that("n_start starts, the first the quantile start", {
  y <- c(1, 2, 3, 10, 11, 12, 20, 21)
  tm <- regime(k = 3, n_start = 4)
  st <- term_starts(tm, target = y)
  expect_length(st, 4L)
  expect_identical(st[[1]], term_start(tm, target = y))
  for (z in st) expect_identical(names(z), term_params(tm))
  # the first level stays, the gaps move and stay positive on their chart,
  # the log-ratios are drawn afresh
  lv <- vapply(st, function(z) z[["level1"]], numeric(1))
  expect_identical(unique(lv), st[[1]][["level1"]])
  alr <- vapply(st, function(z) z[["alr1.1"]], numeric(1))
  expect_length(unique(alr), 4L)
  gp <- vapply(st, function(z) z[["gap2"]], numeric(1))
  expect_length(unique(gp), 4L)
})

test_that("the starts are the same every time and leave the stream alone", {
  y <- c(1, 2, 3, 10, 11, 12)
  tm <- regime(k = 2, n_start = 3)
  set.seed(9)
  a <- term_starts(tm, target = y)
  u1 <- runif(1)
  set.seed(9)
  u2 <- runif(1)
  b <- term_starts(tm, target = y)
  expect_identical(a, b)
  expect_identical(u1, u2)
  # and a session with no stream yet is left without one
  had <- exists(".Random.seed", envir = globalenv())
  if (had) {
    old <- get(".Random.seed", envir = globalenv())
    rm(".Random.seed", envir = globalenv())
    on.exit(assign(".Random.seed", old, envir = globalenv()), add = TRUE)
  }
  term_starts(tm, target = y)
  expect_false(exists(".Random.seed", envir = globalenv()))
})

test_that("n_start must be a whole number of at least one", {
  expect_error(regime(k = 2, n_start = 0), "n_start")
  expect_error(regime(k = 2, n_start = 1.5), "n_start")
  expect_error(regime(k = 2, n_start = NA), "n_start")
  expect_error(regime(k = 2, n_start = c(2, 3)), "n_start")
  expect_identical(regime(k = 2, n_start = 5)@n_start, 5L)
})
