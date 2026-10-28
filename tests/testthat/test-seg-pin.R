# A continuous term with a held coefficient: the free break-points stay on
# their side of a held position, and the held coefficients keep their slot.

pin_data <- function(n = 120L, seed = 1L) {
  set.seed(seed)
  data.frame(x = sort(runif(n, 0, 10)))
}

# where a free break-point can sit beside a held one at h: the k-th
# observation above it, and the k-th at or below it less twice the tolerance
pin_obs <- function(d, h, k = 2L) {
  tie <- sqrt(.Machine$double.eps) * diff(range(d$x))
  xs <- sort(d$x)
  c(left = sort(xs[xs - h <= tie], decreasing = TRUE)[k] - 2 * tie,
    right = xs[xs - h > tie][k])
}

test_that("without a pin the confinement interval is the term's own", {
  d <- pin_data()
  b <- term_build(seg(x, npsi = 2, psi = c(3, 7)), d)
  lim <- as.numeric(stats::quantile(d$x, c(0.05, 0.95), names = FALSE))
  m <- seg_limits(b)
  expect_identical(dimnames(m), list(c("psi1", "psi2"), c("lower", "upper")))
  expect_equal(m[1, ], c(lower = lim[1], upper = lim[2]))
  expect_equal(m[2, ], c(lower = lim[1], upper = lim[2]))
})

test_that("a free break-point is confined to the side of the held one that it starts on", {
  d <- pin_data()
  b <- term_build(seg(x, npsi = 2, psi = c(3, 7)), d)
  lim <- as.numeric(stats::quantile(d$x, c(0.05, 0.95), names = FALSE))
  o <- pin_obs(d, 5)
  # psi2 held at 5, psi1 starts below it
  p <- seg_pin(b, held = c(seg.psi2 = 5), start = c(3, 7))
  m <- seg_limits(p)
  expect_equal(unname(m[1, ]), c(lim[1], o[["left"]]))
  expect_equal(unname(m[2, ]), lim)
  # two observations lie between the free and the held break-point at the
  # ends of the interval, and the third does not
  expect_equal(sum(d$x > o[["left"]] & d$x <= 5), 2L)
  expect_equal(sum(d$x > 5 & d$x <= o[["right"]]), 2L)
  # the same held position with the free one starting above it
  p <- seg_pin(b, held = c(seg.psi2 = 5), start = c(7, 8))
  expect_equal(unname(seg_limits(p)[1, ]), c(o[["right"]], lim[2]))
  # a free break-point that starts too close takes the nearest interval
  p <- seg_pin(b, held = c(seg.psi2 = 5), start = c(5.2, 7))
  expect_equal(unname(seg_limits(p)[1, ]), c(o[["right"]], lim[2]))
  # a longer separation moves both ends out
  p <- seg_pin(b, held = c(seg.psi2 = 5), start = c(3, 7), min_obs = 5L)
  expect_equal(unname(seg_limits(p)[1, 2]), pin_obs(d, 5, k = 5L)[["left"]])
  # a held position beyond the confinement limits is clamped, and there is no
  # room on the side that is not left
  p <- seg_pin(b, held = c(seg.psi2 = 100), start = c(3, 7))
  expect_equal(unname(seg_limits(p)[1, ]),
               c(lim[1], pin_obs(d, lim[2])[["left"]]))
})

test_that("the position of a free break-point is clamped inside its interval", {
  d <- pin_data()
  b <- term_build(seg(x, npsi = 2, psi = c(3, 7)), d)
  p <- seg_pin(b, held = c(seg.psi2 = 5), start = c(3, 7))
  o <- pin_obs(d, 5)
  cf <- b@build_info$coef
  cf[b@build_info$index$psi1] <- 8
  cf[b@build_info$index$psi2] <- 5
  r <- term_refresh(p, cf)
  expect_equal(unname(r@build_info$psi[1, ]), c(o[["left"]], 5))
  # its column of the block is zero there, as at a confinement limit, and the
  # coefficients are not exchanged
  expect_true(all(term_matrix(r)[, "seg.psi1"] == 0))
  expect_identical(r@build_info$coef, cf)
  # without the pin the crossed break-points exchange their coefficients
  r0 <- term_refresh(b, cf)
  expect_equal(unname(r0@build_info$psi[1, ]), c(5, 8))
  expect_false(identical(r0@build_info$coef, cf))
})

test_that("a held change of slope keeps its slot when the positions cross", {
  d <- pin_data()
  b <- term_build(seg(x, npsi = 2, psi = c(3, 7)), d)
  p <- seg_pin(b, held = c(seg.gamma2 = 1.5), start = c(3, 7))
  cf <- b@build_info$coef
  cf[b@build_info$index$gamma1] <- -2
  cf[b@build_info$index$gamma2] <- 1.5
  cf[b@build_info$index$psi1] <- 8
  cf[b@build_info$index$psi2] <- 4
  # no position is held, so the intervals are the term's own
  expect_identical(seg_limits(p), seg_limits(b))
  expect_identical(term_refresh(p, cf)@build_info$coef, cf)
  expect_false(identical(term_refresh(b, cf)@build_info$coef, cf))
})

test_that("free break-points on one side of a held one still exchange their labels", {
  d <- pin_data()
  b <- term_build(seg(x, npsi = 3, psi = c(2, 4, 8)), d)
  p <- seg_pin(b, held = c(seg.psi3 = 8), start = c(2, 4, 8))
  idx <- b@build_info$index
  cf <- b@build_info$coef
  cf[idx$psi1] <- 4
  cf[idx$psi2] <- 2
  cf[idx$gamma1] <- 1
  cf[idx$gamma2] <- -1
  cf[idx$psi3] <- 8
  cf[idx$gamma3] <- 0.5
  r <- term_refresh(p, cf)@build_info$coef
  expect_equal(unname(r[c(idx$psi1, idx$psi2, idx$psi3)]), c(2, 4, 8))
  expect_equal(unname(r[c(idx$gamma1, idx$gamma2, idx$gamma3)]), c(-1, 1, 0.5))
})

test_that("a smoothed term is pinned like a sharp one", {
  d <- pin_data()
  b <- term_build(seg(x, npsi = 2, psi = c(3, 7),
                      smoothed = numericals7::smooth_probit()), d)
  p <- seg_pin(b, held = c(seg.psi2 = 5), start = c(3, 7))
  expect_lt(seg_limits(p)[1, "upper"], 5)
  cf <- b@build_info$coef
  cf[b@build_info$index$psi1] <- 8
  cf[b@build_info$index$psi2] <- 5
  expect_identical(term_refresh(p, cf)@build_info$coef, cf)
})

test_that("seg_pin rejects what it cannot pin", {
  d <- pin_data()
  b <- term_build(seg(x, npsi = 2, psi = c(3, 7)), d)
  expect_error(seg_pin(b, held = c(seg.nothing = 5)), "distinct coefficients")
  expect_error(seg_pin(b, held = 5), "named numeric")
  expect_error(seg_pin(b, held = c(seg.psi2 = 5), min_obs = 0), "'min_obs'")
  expect_error(seg_pin(b, held = c(seg.psi2 = 5), min_obs = 1.5), "'min_obs'")
  expect_error(seg_pin(b, held = c(seg.psi2 = 5), start = 1), "'start'")
  expect_error(seg_pin(term_build(jump(x, npsi = 2, psi = c(3, 7)), d),
                       held = c(jump.psi2 = 5)), "seg_hold")
  expect_error(seg_pin(1, held = c(seg.psi2 = 5)), "break-point term")
  # a separation so long that the held position leaves no room on either side
  expect_error(seg_pin(b, held = c(seg.psi2 = 5), min_obs = 80L), "no room")
})
