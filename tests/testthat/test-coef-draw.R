make_data <- function(n_id = 12, seed = 1) {
  set.seed(seed)
  nj <- 3 + rpois(n_id, 60)
  id <- rep(seq_len(n_id), nj)
  x <- unlist(lapply(nj, function(m) sort(runif(m))))
  data.frame(id = id, x = x, g = factor(id))
}

psi_of <- function(tm, cf) {
  p <- term_refresh(tm, cf)@blueprint$psi
  apply(p, 2, mean)
}

# EVERY position and not their average: under a partition each group has its
# own, and a mean over groups pinned at both limits sits comfortably in the
# middle while not one of them is where it should be.
psi_cells <- function(tm, cf) {
  as.numeric(term_refresh(tm, cf)@blueprint$psi)
}

test_that("a term with no opinion hands the draw back untouched", {
  dd <- make_data()
  b <- term_build(linpar(~ x), dd)
  d <- term_coef_draw(b, c(1, 2))
  expect_identical(d$coef, c(1, 2))
  expect_length(d$index, 0L)
  expect_length(d$scale, 0L)
})

test_that("a break-point is drawn inside the covariate, where a plain draw is not", {
  dd <- make_data()
  shapes <- list(
    "seg" = seg(x, npsi = 1),
    "seg, three" = seg(x, npsi = 3),
    "seg, developed" = seg(x, by = ~ random(~ 1 | id), npsi = 1),
    "seg, partition" = seg(x, psi ~ 0 + g, npsi = 1),
    "jump" = jump(x, npsi = 1),
    "jump, two" = jump(x, npsi = 2),
    "jump, partition" = jump(x, by = ~ 0 + g, npsi = 1),
    "jseg" = jseg(x, npsi = 1))
  for (nm in names(shapes)) {
    tm <- term_build(shapes[[nm]], dd)
    lim <- tm@blueprint$lim
    set.seed(3)
    R <- 60
    drawn <- unlist(replicate(R, psi_cells(tm, term_coef_draw(
      tm, stats::rnorm(term_npar(tm)))$coef), simplify = FALSE))
    # NOT ONE of them is pinned against a confinement limit
    expect_true(all(drawn > lim[1L] + 1e-9), info = nm)
    expect_true(all(drawn < lim[2L] - 1e-9), info = nm)
    # THE NEGATIVE CONTROL, without which the assertion above is a
    # statement about the confinement rather than about the draw: the
    # plain normal this replaces does end up pinned.
    set.seed(3)
    plain <- unlist(replicate(R, psi_cells(tm, stats::rnorm(term_npar(tm))),
                              simplify = FALSE))
    expect_true(any(plain <= lim[1L] + 1e-9 | plain >= lim[2L] - 1e-9),
                info = nm)
  }
})

test_that("the position is written exactly where the construction holds it", {
  dd <- make_data()
  set.seed(5)
  for (kind in c("seg", "jump")) {
    tm <- term_build(get(kind)(x, npsi = 1), dd)
    ix <- tm@blueprint$index[["psi1"]]
    gap <- replicate(40, {
      cf <- term_coef_draw(tm, stats::rnorm(term_npar(tm)))$coef
      aim <- if (kind == "seg") cf[ix] else
        -cf[ix] / cf[tm@blueprint$index[["delta1"]]]
      abs(seg_psi(term_refresh(tm, cf)) - aim)
    })
    expect_lt(max(gap), 1e-12)
  }
})

test_that("the draw spreads from the covariate and not from where the term sits", {
  dd <- make_data()
  tm <- term_build(seg(x, npsi = 1), dd)
  lim <- tm@blueprint$lim
  # A CALLER MAY HAND OVER A TERM SITTING ANYWHERE: a simulation builds its
  # specification against a placeholder response, and a break-point chosen
  # on a least-squares profile of noise lands against a limit as readily as
  # anywhere. The draw must not inherit that.
  moved <- seg_relocate(tm, lim[1L])
  expect_equal(as.numeric(seg_psi(moved)), lim[1L], tolerance = 1e-9)
  set.seed(11)
  drawn <- replicate(60, psi_of(moved, term_coef_draw(
    moved, stats::rnorm(term_npar(moved)))$coef))
  q0 <- as.numeric(stats::quantile(dd$x, 0.5, names = FALSE))
  expect_equal(mean(drawn), q0, tolerance = 0.05)
  expect_true(all(drawn > lim[1L] + 1e-9))
})

test_that("sd = 0 places the break-points at the covariate's interior quantiles", {
  dd <- make_data()
  tm <- term_build(seg(x, npsi = 3), dd)
  cf <- term_coef_draw(tm, stats::rnorm(term_npar(tm)), sd = 0)$coef
  expect_equal(as.numeric(seg_psi(term_refresh(tm, cf))),
               as.numeric(stats::quantile(dd$x, c(1, 2, 3) / 4, names = FALSE)),
               tolerance = 1e-9)
})

test_that("the reported scale is the spread the positions get, and sd scales it", {
  dd <- make_data()
  tm <- term_build(seg(x, npsi = 1), dd)
  set.seed(2)
  d <- term_coef_draw(tm, stats::rnorm(term_npar(tm)))
  expect_length(d$index, 1L)
  got <- replicate(400, psi_of(tm, term_coef_draw(
    tm, stats::rnorm(term_npar(tm)))$coef))
  expect_equal(stats::sd(got), d$scale[[1L]], tolerance = 0.15)
  wide <- replicate(400, psi_of(tm, term_coef_draw(
    tm, stats::rnorm(term_npar(tm)), sd = 2)$coef))
  expect_equal(stats::sd(wide) / stats::sd(got), 2, tolerance = 0.2)
})

test_that("the width is the derived one", {
  dd <- make_data()
  for (K in 1:3) {
    tm <- term_build(seg(x, npsi = K), dd)
    lim <- tm@blueprint$lim
    d <- term_coef_draw(tm, stats::rnorm(term_npar(tm)))
    expect_equal(d$scale[[1L]], (lim[2L] - lim[1L]) / (6 * (K + 1)),
                 tolerance = 1e-12)
  }
})

test_that("a bad width is refused and a bad coefficient count with it", {
  dd <- make_data()
  tm <- term_build(seg(x, npsi = 1), dd)
  expect_error(term_coef_draw(tm, stats::rnorm(term_npar(tm)), sd = -1),
               "finite non-negative")
  expect_error(term_coef_draw(tm, c(1, 2)), "one per column")
})
