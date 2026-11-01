# scs() with a factor or a numeric `by`: the factor term against the separate
# curves with the same coefficients, the numeric term against z f(x), the
# block and its derivatives against central differences, the penalties, the
# versions searched one level at a time, and the rejections.

set.seed(41)
by_d <- data.frame(x = runif(150), g = gl(3, 50, labels = c("a", "b", "c")),
                   z = rnorm(150))
by_sm <- basis7::bspline_smooth(k = 7)

by_fd <- function(tm, g, h = 1e-5) {
  p <- length(g); n <- length(term_value(tm))
  e <- function(j) replace(numeric(p), j, h)
  v <- stats::rnorm(p); u <- stats::rnorm(p); A <- matrix(stats::rnorm(n * p), n)
  M <- function(gg) term_matrix(term_refresh(tm, gg))
  list(J = sapply(seq_len(p), function(j) (term_value(tm, coef = g + e(j)) -
                                            term_value(tm, coef = g - e(j))) / (2 * h)),
       v = v, u = u, A = A,
       D1 = (M(g + h * v) - M(g - h * v)) / (2 * h),
       D2 = (term_block_deriv(tm, coef = g + h * u, v = v) -
             term_block_deriv(tm, coef = g - h * u, v = v)) / (2 * h),
       C = sapply(seq_len(p), function(k) sum(A * (M(g + e(k)) - M(g - e(k)))) / (2 * h)))
}

test_that("a factor by is the separate curves of the levels", {
  set.seed(42)
  for (sh in list(shape_increasing(), shape_concave(),
                  shape_increasing_convex(value = c(at = 0.4, value = 1)),
                  shape_convex(slope = c(at = 0.6, value = 0)))) {
    tm <- term_build(scs(x, sh, by_sm, by = g), by_d)
    one <- term_build(scs(x, sh, by_sm), by_d)
    nc <- length(term_coef_names(one))
    expect_length(term_coef_names(tm), 3L * nc)
    gam <- stats::rnorm(3L * nc, sd = 0.7)
    sep <- sapply(1:3, function(l) term_value(one, coef = gam[(l - 1L) * nc + seq_len(nc)]))
    expect_equal(term_value(term_refresh(tm, gam)),
                 sep[cbind(seq_len(nrow(by_d)), as.integer(by_d$g))], tolerance = 1e-13)
  }
  tm <- term_build(scs(x, shape_increasing(), by_sm, by = g), by_d)
  expect_identical(head(term_coef_names(tm), 2L), c("scs(x).a.g1", "scs(x).a.g2"))
})

test_that("a numeric by is z times a shaped coefficient with a free level", {
  set.seed(43)
  tm <- term_build(scs(x, shape_increasing(), by_sm, by = z), by_d)
  expect_identical(tm@label, "scs(x):z")
  expect_identical(term_coef_names(tm)[1:2], c("scs(x):z.const", "scs(x):z.g1"))
  gam <- stats::rnorm(length(term_coef_names(tm)), sd = 0.7)
  f1 <- term_value(tm, coef = gam, newdata = transform(by_d, z = 1))
  expect_equal(term_value(term_refresh(tm, gam)), by_d$z * f1, tolerance = 1e-13)
  grid <- data.frame(x = seq(0, 1, length.out = 301), z = 1)
  expect_true(all(diff(term_value(tm, coef = gam, newdata = grid)) >= -1e-12))
  # the level moves f by a constant
  f2 <- term_value(tm, coef = replace(gam, 1L, gam[1L] + 2), newdata = grid)
  expect_equal(f2 - term_value(tm, coef = gam, newdata = grid), rep(2, 301), tolerance = 1e-13)
  # a value at a point fixes the level: one coefficient fewer, f(a) = v
  tv <- term_build(scs(x, shape_increasing(value = c(at = 0.5, value = 0.3)), by_sm, by = z), by_d)
  expect_length(term_coef_names(tv), length(term_coef_names(tm)) - 1L)
  gv <- stats::rnorm(length(term_coef_names(tv)))
  expect_equal(term_value(tv, coef = gv, newdata = data.frame(x = 0.5, z = 1)), 0.3,
               tolerance = 1e-13)
})

test_that("the block and its derivatives agree with central differences", {
  set.seed(44)
  terms <- list(scs(x, shape_increasing(), by_sm, by = g),
                scs(x, shape_convex(slope = c(at = 0.3, value = 1)), by_sm, by = g),
                scs(x, shape_increasing_concave(), by_sm, by = z),
                scs(x, shape_concave(slope = c(at = 0.3, value = 0)), by_sm, by = z))
  for (t0 in terms) {
    tm <- term_build(t0, by_d)
    g <- stats::rnorm(length(term_coef_names(tm)), sd = 0.6)
    r <- term_refresh(tm, g)
    fd <- by_fd(tm, g)
    rel <- function(a, b) max(abs(a - b)) / max(abs(b))
    expect_lt(rel(term_matrix(r), fd$J), 1e-8)
    expect_lt(rel(term_block_deriv(tm, coef = g, v = fd$v), fd$D1), 1e-8)
    expect_lt(rel(term_block_deriv2(tm, coef = g, v = fd$v, u = fd$u), fd$D2), 1e-8)
    expect_lt(rel(term_block_contract(tm, coef = g, A = fd$A), fd$C), 1e-8)
    expect_equal(term_predict(r, by_d), term_matrix(r), tolerance = 1e-14,
                 ignore_attr = TRUE)
  }
})

test_that("the penalty is one block per level, shared or per level", {
  tm <- term_build(scs(x, shape_increasing(), by_sm, by = g), by_d)
  one <- term_build(scs(x, shape_increasing(), by_sm), by_d)
  P1 <- penalties7::penalty_hessian(term_penalty(one), numeric(ncol(term_matrix(one))),
                                    c(lambda = 1))
  P3 <- penalties7::penalty_hessian(term_penalty(tm), numeric(ncol(term_matrix(tm))),
                                    c(lambda = 1))
  expect_equal(as.matrix(P3), as.matrix(kronecker(diag(3), as.matrix(P1))),
               ignore_attr = TRUE)
  lv <- term_build(scs(x, shape_increasing(), by_sm, by = g, by_hyper = "level",
                       hyper = c(lambda.b = 2)), by_d)
  ent <- term_penalties(lv)
  expect_identical(vapply(ent, `[[`, "", "name"), c("a", "b", "c"))
  nc <- ncol(term_matrix(one))
  expect_identical(ent[[2L]]$index, nc + seq_len(nc))
  expect_equal(unlist(ent[[2L]]$fixed), c(lambda = 2))
})

test_that("a located shape under a factor by is searched one level at a time", {
  tm <- term_build(scs(x, shape_unimodal(), by_sm, by = g), by_d)
  expect_identical(term_variant_stages(tm), 3L)
  expect_identical(term_variant_stages(term_build(scs(x, shape_unimodal(), by_sm, by = z), by_d)), 1L)
  expect_identical(term_variant_stages(term_build(scs(x, shape_unimodal(mode = 0.5), by_sm, by = g), by_d)), 1L)
  expect_identical(term_variant_stages(s(x)), 1L)
  mo <- c(a = 0.2, b = 0.55, c = 0.8)
  pil <- -(by_d$x - mo[as.character(by_d$g)])^2
  xi <- tm@build_info$xi
  rng <- .shape_index_range(tm@shape, length(xi))
  near <- vapply(mo, function(m) rng[which.min(abs(xi[rng] - m))], 1L)
  v1 <- term_variants(tm, pilot = pil, stage = 1L)
  expect_identical(names(v1), sprintf("a: M = %d", intersect(near[["a"]] + (-2):2, rng)))
  expect_identical(unname(v1[[1L]]@spec$index_by[c("b", "c")]), unname(near[c("b", "c")]))
  cur <- v1[[1L]]
  v2 <- term_variants(tm, pilot = pil, stage = 2L, current = cur)
  expect_false(any(vapply(v2, function(v) v@spec$index_by[["b"]], 1L) == cur@spec$index_by[["b"]]))
  expect_true(all(vapply(v2, function(v) v@spec$index_by[["a"]], 1L) == cur@spec$index_by[["a"]]))
  # each level is built with its own index
  b2 <- term_build(v2[[1L]], by_d)
  expect_identical(length(b2@build_info$pieces) >= 2L, TRUE)
  # the pilot is s() with the same by
  pt <- term_pilot(tm)
  expect_true(S7::S7_inherits(pt, SmoothTerm))
  expect_identical(deparse(pt@by), "g")
})

test_that("a by that cannot be built is rejected", {
  expect_error(term_build(scs(x, by = z, by_hyper = "level"), by_d), "no factor 'by'")
  d2 <- by_d; d2$g[3] <- NA
  expect_error(term_build(scs(x, by = g), d2), "missing or infinite")
  tm <- term_build(scs(x, by = g), by_d)
  expect_error(term_predict(tm, data.frame(x = 0.5, g = "d")), "absent at build")
  tz <- term_build(scs(x, by = z), by_d)
  expect_error(term_predict(tz, data.frame(x = 0.5, z = "a")), "was numeric")
})

test_that("edf() reads a penalty repeated in blocks", {
  set.seed(45)
  for (t0 in list(scs(x, shape_increasing(), by_sm, by = g),
                  s(x, smoother = by_sm, by = g))) {
    tm <- term_build(t0, by_d)
    k <- length(term_coef_names(tm))
    A <- matrix(stats::rnorm(k * k), k); H <- crossprod(A) / k + diag(k)
    cf <- if (S7::S7_inherits(tm, ScsTerm)) tm@build_info$coef else numeric(k)
    S <- as.matrix(penalties7::penalty_hessian(term_penalty(tm), cf, list(lambda = 0.7)))
    expect_equal(edf(tm, coef = cf, hessian = H, theta = list(lambda = 0.7)),
                 sum(diag(solve(H + S, H))), tolerance = 1e-10)
  }
})

test_that("a tied left slope is one coefficient for every level", {
  set.seed(46)
  for (sh in list(shape_convex(left_slope_by = "shared"),
                  shape_increasing_convex(left_slope_by = "shared"))) {
    tm <- term_build(scs(x, sh, by_sm, by = g), by_d)
    bp <- tm@build_info
    expect_length(term_coef_names(tm), 1L + 3L * (bp$nc - 1L))
    expect_identical(term_coef_names(tm)[1:2], c("scs(x).g1", "scs(x).a.g2"))
    gam <- stats::rnorm(length(term_coef_names(tm)), sd = 0.6)
    sh0 <- sh; sh0@shape_params$left_slope_by <- "level"
    untied <- term_build(scs(x, sh0, by_sm, by = g), by_d)
    expect_equal(term_value(term_refresh(tm, gam)),
                 term_value(untied, coef = gam[bp$ix]), tolerance = 1e-13)
    fd <- by_fd(tm, gam)
    rel <- function(a, b) max(abs(a - b)) / max(abs(b))
    r <- term_refresh(tm, gam)
    expect_lt(rel(term_matrix(r), fd$J), 1e-8)
    expect_lt(rel(term_block_deriv(tm, coef = gam, v = fd$v), fd$D1), 1e-8)
    expect_lt(rel(term_block_deriv2(tm, coef = gam, v = fd$v, u = fd$u), fd$D2), 1e-8)
    expect_lt(rel(term_block_contract(tm, coef = gam, A = fd$A), fd$C), 1e-8)
  }
  # the penalty is the levels' penalties on the coefficients
  tm <- term_build(scs(x, shape_increasing_convex(left_slope_by = "shared"), by_sm, by = g), by_d)
  one <- term_build(scs(x, shape_increasing_convex(), by_sm), by_d)
  P1 <- as.matrix(penalties7::penalty_hessian(term_penalty(one), numeric(ncol(term_matrix(one))),
                                              list(lambda = 1)))
  C <- tm@build_info$C
  expect_equal(as.matrix(penalties7::penalty_hessian(term_penalty(tm), tm@build_info$coef,
                                                     list(lambda = 1))),
               crossprod(C, kronecker(diag(3), P1) %*% C), ignore_attr = TRUE)
  expect_error(term_build(scs(x, shape_increasing_convex(left_slope_by = "shared"),
                              by_sm, by = g, by_hyper = "level"), by_d), "once per")
  expect_length(term_penalties(term_build(scs(x, shape_convex(left_slope_by = "shared"),
                                              by_sm, by = g, by_hyper = "level"), by_d)), 3L)
})

test_that("the declared settings of a shape may differ by level", {
  set.seed(47)
  lv <- factor(c("a", "b", "c"))
  at <- c(a = 0.2, b = 0.5, c = 0.8); v <- c(a = 0, b = 1, c = 2)
  tm <- term_build(scs(x, shape_increasing(value = list(at = at, value = v)), by_sm, by = g),
                   by_d)
  gam <- stats::rnorm(length(term_coef_names(tm)))
  for (l in names(at)) {
    expect_equal(term_value(tm, coef = gam, newdata = data.frame(x = at[[l]], g = factor(l, levels(lv)))),
                 v[[l]], tolerance = 1e-13)
  }
  ms <- term_build(scs(x, shape_increasing(min_slope = c(a = 0, b = 0.5, c = 1)), by_sm, by = g),
                   by_d)
  expect_identical(length(ms@build_info$pieces), 3L)
  # each level is the term without by built with that level's settings
  for (l in c("a", "c")) {
    one <- term_build(scs(x, shape_increasing(min_slope = c(a = 0, b = 0.5, c = 1)[[l]]), by_sm), by_d)
    k <- ncol(term_matrix(one)); j <- (match(l, c("a", "b", "c")) - 1L) * k + seq_len(k)
    g2 <- stats::rnorm(3L * k)
    nd <- data.frame(x = by_d$x, g = factor(l, levels(lv)))
    expect_equal(term_value(ms, coef = g2, newdata = nd), term_value(one, coef = g2[j], newdata = nd),
                 tolerance = 1e-13)
  }
  um <- term_build(scs(x, shape_unimodal(mode = c(a = 0.3, b = 0.5, c = 0.7),
                                         valley = c(a = FALSE, b = TRUE, c = FALSE)), by_sm, by = g),
                   by_d)
  expect_identical(vapply(um@build_info$pieces, function(pc) pc$shape@shape_params$sign, 1),
                   c(1, -1, 1))
  expect_output(print(shape_increasing(min_slope = c(a = 0, b = 0.2))), "min_slope = (a = 0.0, b = 0.2)",
                fixed = TRUE)
})

test_that("a location may be one for every level", {
  tm <- term_build(scs(x, shape_unimodal(mode_by = "shared"), by_sm, by = g), by_d)
  expect_identical(term_variant_stages(tm), 1L)
  mo <- c(a = 0.2, b = 0.55, c = 0.8)
  pil <- -(by_d$x - mo[as.character(by_d$g)])^2
  v <- term_variants(tm, pilot = pil)
  expect_true(all(grepl("^M = ", names(v))))
  b1 <- term_build(v[[1L]], by_d)
  expect_identical(length(b1@build_info$pieces), 1L)
  expect_output(print(shape_sigmoid(inflection_by = "shared")), "one location for every level")
})

test_that("settings by level that cannot be resolved are rejected", {
  expect_error(term_build(scs(x, shape_increasing(min_slope = c(a = 0, b = 1))), by_d),
               "need a\n  factor 'by'|need a factor 'by'")
  expect_error(term_build(scs(x, shape_increasing(min_slope = c(a = 0, b = 1)), by = g), by_d),
               "no value for the level 'c'")
  expect_error(term_build(scs(x, shape_increasing(min_slope = c(a = 0, b = 1, c = 0, z = 2)),
                              by = g), by_d), "'z'")
  expect_error(term_build(scs(x, shape_unimodal(mode_by = "shared")), by_d), "needs a factor 'by'")
  expect_error(term_build(scs(x, shape_unimodal(mode_by = "shared"), by = z), by_d),
               "needs a factor 'by'")
  expect_error(term_build(scs(x, shape_convex(slope = c(at = 0.5, value = 0),
                                              left_slope_by = "shared"), by = g), by_d),
               "removes the left")
  expect_error(shape_increasing(min_slope = c(0, 1)), "named by the levels")
  expect_error(shape_unimodal(mode_by = "both"), "\"level\" or \"shared\"")
  expect_error(shape_increasing(value = list(at = 0)), "'at' and")
  expect_error(shape_map(shape_increasing(min_slope = c(a = 0, b = 1)), seq(0, 1, length.out = 5)),
               "values by level")
})

test_that("a given location under a by has its two starts, except under a numeric by", {
  tm <- term_build(scs(x, shape_unimodal(mode = c(a = 0.3, b = 0.5, c = 0.7)), by_sm, by = g), by_d)
  expect_identical(names(term_variants(tm)), c("start: sides", "start: spline"))
  tz <- term_build(scs(x, shape_unimodal(mode = 0.5), by_sm, by = z), by_d)
  expect_null(term_pilot(tz))
  expect_length(term_variants(tz), 1L)
})
