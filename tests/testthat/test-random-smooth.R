rs_data <- function(n = 120, m = 6, seed = 1) {
  set.seed(seed)
  data.frame(x = runif(n), z = runif(n),
             g = factor(rep(letters[seq_len(m)], length.out = n)))
}

test_that("random() builds a smooth random effect from s() or te()", {
  tm <- random(~ s(x, smoother = basis7::bspline_smooth(k = 6)) | g)
  expect_true(S7::S7_inherits(tm, RandomSmoothTerm))
  expect_true(S7::S7_inherits(tm, RandomTerm))
  expect_true(tm@intercept)
  expect_false(random(~ 0 + s(x) | g)@intercept)
  expect_false(random(~ s(x) - 1 | g)@intercept)
  expect_true(random(~ 1 + s(x) | g)@intercept)
  expect_true(S7::S7_inherits(random(~ te(x, z) | g), RandomSmoothTerm))
  # an ordinary bar formula is unchanged
  expect_false(S7::S7_inherits(random(~ x | g), RandomSmoothTerm))
})

test_that("the left side of a smooth random effect is checked", {
  expect_error(random(~ s(x) + s(z) | g), "more than one smooth")
  expect_error(random(~ s(x) + z | g), "random\\(~ z \\| g\\)")
  expect_error(random(~ s(x, by = z) | g), "'by'")
  expect_error(random(~ s(x, hyper = c(lambda = 1)) | g), "'hyper' or 'id'")
  expect_error(random(~ s(x) | u | g), "covariance label")
})

test_that("the block, its names and its two priors", {
  dd <- rs_data()
  b <- term_build(random(~ s(x, smoother = basis7::bspline_smooth(k = 6)) | g),
                  dd)
  gr <- term_group(b)
  expect_identical(gr$names, c("(Intercept)", "lin", paste0("z", 1:4)))
  expect_identical(gr$dim, 6L)
  expect_identical(term_npar(b), 36L)
  expect_s4_class(term_matrix(b), "dgCMatrix")
  ent <- term_penalties(b)
  expect_length(ent, 2L)
  # the random level and the line share an unstructured covariance
  expect_identical(ent[[1L]]$index, as.integer(outer(1:2, (0:5) * 6, `+`)))
  expect_equal(ent[[1L]]$penalty@block, 2)
  # the penalized coordinates of every group have one standard deviation
  expect_identical(ent[[2L]]$name, "smooth")
  expect_identical(ent[[2L]]$index, as.integer(outer(3:6, (0:5) * 6, `+`)))
  expect_identical(ent[[2L]]$penalty@params, "sigma")
})

test_that("a Demmler-Reinsch smooth keeps its own coordinates", {
  dd <- rs_data()
  sm <- basis7::bspline_smooth(k = 6)
  b <- term_build(random(~ s(x, smoother = sm) | g), dd)
  X <- as.matrix(term_matrix(term_build(s(x, smoother = sm), dd)))
  W <- term_within(b, dd)
  expect_equal(unname(W[, -1L]), unname(X), tolerance = 1e-13)
  expect_equal(unname(W[, 1L]), rep(1, nrow(dd)))
})

test_that("the penalized coordinates are whitened", {
  dd <- rs_data()
  for (sm in list(basis7::bspline_smooth(k = 6, null_space = "shrink"),
                  basis7::bspline_smooth(k = 6, reparam = "none"),
                  basis7::pspline_smooth(k = 8, reparam = "orthonorm"))) {
    sb <- term_build(s(x, smoother = sm), dd)
    cp <- .random_smooth_components(sb)
    sp <- .random_smooth_split(cp, paste0("c", seq_len(nrow(cp[[1L]]))))
    Tp <- sp$T[, sp$nnull + seq_len(sp$npen), drop = FALSE]
    expect_equal(crossprod(Tp, cp[[1L]] %*% Tp), diag(1, sp$npen),
                 tolerance = 1e-8)
    if (sp$nnull) {
      Tn <- sp$T[, seq_len(sp$nnull), drop = FALSE]
      expect_lt(max(abs(cp[[1L]] %*% Tn)), 1e-10)
    }
  }
})

test_that("null_space decides what the null space contributes", {
  dd <- rs_data()
  nm <- function(ns) term_group(term_build(random(
    ~ s(x, smoother = basis7::bspline_smooth(k = 6, null_space = ns)) | g),
    dd))$names
  expect_identical(nm("keep"), c("(Intercept)", "lin", paste0("z", 1:4)))
  expect_identical(nm("drop"), c("(Intercept)", paste0("z", 1:4)))
  b <- term_build(random(
    ~ s(x, smoother = basis7::bspline_smooth(k = 6, null_space = "shrink")) |
      g), dd)
  ent <- term_penalties(b)
  # under shrink the line is penalized with the rest, so the first prior is
  # the random level alone
  expect_identical(ent[[1L]]$penalty@params, "sigma")
  expect_length(ent[[1L]]$index, 6L)
  expect_length(ent[[2L]]$index, 30L)
})

test_that("the null space is the operator's, whatever its form", {
  dd <- rs_data()
  four <- term_build(random(~ s(x, smoother = basis7::fourier_smooth(
    k = 7, lower = 0, upper = 1, null_space = "keep")) | g), dd)
  expect_identical(term_group(four)$names[1:3],
                   c("(Intercept)", "sin1", "cos1"))
  cyc <- term_build(random(~ s(x, smoother = basis7::cyclic_smooth(
    k = 6, lower = 0, upper = 1)) | g), dd)
  expect_identical(term_penalties(cyc)[[1L]]$penalty@params, "sigma")
  expect_length(term_penalties(cyc)[[1L]]$index, 6L)
  ex <- term_build(random(~ s(x, smoother = basis7::bspline_smooth(
    k = 6, order = basis7::linear_operator(c(0, 2)))) | g), dd)
  expect_match(term_group(ex)$names[2L], "exp")
})

test_that("a tensor product's null space is rotated out", {
  dd <- rs_data(n = 200)
  b <- term_build(random(~ te(x, z, smooths = basis7::bspline_smooth(k = 4)) |
                           g), dd)
  gr <- term_group(b)
  # x, z and xz: the product of the margins' lines, less the constant
  expect_identical(gr$names[1:4], c("(Intercept)", "null1", "null2", "null3"))
  ent <- term_penalties(b)
  expect_identical(ent[[2L]]$penalty@params, c("lambda1", "lambda2"))
  expect_identical(ent[[2L]]$penalty@readable, "sd")
  expect_identical(penalties7::penalty_rank(ent[[2L]]$penalty),
                   ent[[2L]]$penalty@n_coef)
})

test_that("an adaptive smoother keeps one profile for every group", {
  dd <- rs_data()
  b <- term_build(random(~ s(x, smoother = basis7::adaptive_smooth(
    k = 10, m = 3)) | g), dd)
  pen <- term_penalties(b)[[2L]]$penalty
  expect_identical(pen@params, c("lambda1", "lambda2", "lambda3"))
  expect_identical(penalties7:::additive_blocks(pen), 6L)
})

test_that("prediction reapplies the map and reads a new level as zero", {
  dd <- rs_data()
  b <- term_build(random(~ te(x, z, smooths = basis7::bspline_smooth(k = 4)) |
                           g), dd)
  expect_equal(as.matrix(term_predict(b, dd[1:7, ])),
               as.matrix(term_matrix(b))[1:7, ], tolerance = 1e-12)
  nd <- dd[1:3, ]
  nd$g <- "new"
  expect_error(term_predict(b, nd), "not present at build time")
  expect_true(all(as.matrix(term_predict(b, nd, unseen = "zero")) == 0))
  expect_identical(colnames(term_within(b, nd)), term_group(b)$names)
})

test_that("distrib and correlated reach the null-space prior", {
  dd <- rs_data()
  sm <- basis7::bspline_smooth(k = 6)
  ind <- term_build(random(~ s(x, smoother = sm) | g, correlated = FALSE), dd)
  expect_identical(vapply(term_penalties(ind), `[[`, "", "name"),
                   c("(Intercept)", "lin", "smooth"))
  t4 <- distributions7::fixed(distributions7::student_t1_distrib(), mu = 0,
                              nu = 4)
  # a univariate prior is a template, one copy per null-space column
  tb <- term_build(random(~ s(x, smoother = sm) | g, distrib = t4), dd)
  expect_identical(vapply(term_penalties(tb), `[[`, "", "name"),
                   c("(Intercept)", "lin", "smooth"))
  expect_error(term_build(random(~ 0 + s(x, smoother = basis7::bspline_smooth(
    k = 6, null_space = "drop")) | g, distrib = t4), dd), "has neither")
})
