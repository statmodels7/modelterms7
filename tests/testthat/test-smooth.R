# Smooth terms: the Demmler-Reinsch construction, its penalty, the by
# argument, and the tensor product.

set.seed(11)
n <- 120
dd <- data.frame(x = sort(runif(n)), z = runif(n),
                 g = factor(rep(c("a", "b", "c"), length.out = n)),
                 w = rnorm(n))
dd$y <- sin(2 * pi * dd$x) + rnorm(n, sd = 0.2)

test_that("s() separates the linear effect and penalizes only the deviation", {
  built <- term_build(s(x, k = 10), dd)
  Z <- term_matrix(built)
  cn <- term_coef_names(built)

  expect_identical(cn[1], "s(x).lin")
  expect_true(all(grepl("^s\\(x\\)\\.z", cn[-1])))

  # the Demmler-Reinsch block is empirically orthogonal to a constant and to
  # the covariate, which is what makes the split unbiased
  nl <- Z[, -1, drop = FALSE]
  expect_lt(max(abs(crossprod(cbind(1, dd$x), nl))), 1e-8)
  # and its crossproduct is diagonal
  C <- crossprod(nl)
  expect_lt(max(abs(C[upper.tri(C)])), 1e-8)

  # the penalty is the identity on the deviation and zero on the linear part
  pen <- term_penalty(built)
  P <- penalties7::penalty_matrix(pen, list(lambda = 1))
  expect_equal(unname(diag(P)), c(0, rep(1, ncol(Z) - 1)))
  expect_lt(max(abs(P - diag(diag(P)))), 1e-12)
  expect_identical(pen@n_coef, ncol(Z))
  expect_true(term_smooth(built))
})

test_that("the edf runs from the basis dimension down to the straight line", {
  built <- term_build(s(x, k = 10), dd)
  H <- crossprod(term_matrix(built))
  b <- rep(0, term_npar(built))
  e <- vapply(c(1e-8, 1, 1e8), function(lam)
    edf(built, coef = b, hessian = H, theta = list(lambda = lam)), numeric(1))
  expect_equal(e[1], term_npar(built), tolerance = 1e-4)
  expect_equal(e[3], 1, tolerance = 1e-4)   # the unpenalized linear effect
  expect_true(e[1] > e[2] && e[2] > e[3])
})

test_that("prediction reapplies the stored transform", {
  built <- term_build(s(x, k = 8), dd)
  res <- check_term(s(x, k = 8), dd, verbose = FALSE)
  expect_true(all(res$status == "OK"),
              info = paste(res$check[res$status != "OK"], collapse = ", "))

  # a fit recovers the function it was given, which is what says the
  # construction is usable and not merely well shaped
  Z <- term_matrix(built)
  P <- penalties7::penalty_matrix(term_penalty(built), list(lambda = 1))
  bhat <- solve(crossprod(Z) + 1e-4 * P, crossprod(Z, dd$y))
  fitted <- as.numeric(Z %*% bhat)
  expect_gt(stats::cor(fitted, sin(2 * pi * dd$x)), 0.98)

  # and predicting at new points inside the range follows the same curve
  nd <- data.frame(x = seq(min(dd$x), max(dd$x), length.out = 50))
  pred <- as.numeric(term_predict(built, nd) %*% bhat)
  expect_gt(stats::cor(pred, sin(2 * pi * nd$x)), 0.98)
})

test_that("linear = FALSE leaves the deviation alone", {
  built <- term_build(s(x, k = 8, linear = FALSE), dd)
  expect_false(any(grepl("lin$", term_coef_names(built))))
  P <- penalties7::penalty_matrix(term_penalty(built), list(lambda = 1))
  expect_equal(unname(diag(P)), rep(1, term_npar(built)))
})

test_that("by a factor gives one smooth per level with a shared parameter", {
  built <- term_build(s(x, by = g, k = 6), dd)
  one <- term_build(s(x, k = 6), dd)
  expect_identical(term_npar(built), 3L * term_npar(one))
  expect_true(all(grepl("^s\\(x\\)\\.(a|b|c)\\.", term_coef_names(built))))

  # each level's columns vanish off its own rows
  Z <- term_matrix(built)
  cols <- term_npar(one)
  for (l in seq_len(3)) {
    blk <- Z[, (l - 1) * cols + seq_len(cols), drop = FALSE]
    expect_true(all(blk[dd$g != levels(dd$g)[l], ] == 0))
  }
  # the penalty repeats blockwise, so one lambda governs every level
  expect_identical(term_penalty(built)@params, "lambda")
  res <- check_term(s(x, by = g, k = 6), dd, verbose = FALSE)
  expect_true(all(res$status == "OK"),
              info = paste(res$check[res$status != "OK"], collapse = ", "))
})

test_that("by a numeric is a varying-coefficient term", {
  built <- term_build(s(x, by = w, k = 6), dd)
  plain <- term_build(s(x, k = 6), dd)
  expect_identical(term_npar(built), term_npar(plain))
  expect_equal(term_matrix(built), dd$w * term_matrix(plain),
               ignore_attr = TRUE)
  res <- check_term(s(x, by = w, k = 6), dd, verbose = FALSE)
  expect_true(all(res$status == "OK"))
})

test_that("te() smooths each margin with a parameter of its own", {
  built <- term_build(te(x, z, k = 4), dd)
  # the product of the marginal dimensions less the centering constraint
  expect_identical(term_npar(built), 15L)

  pen <- term_penalty(built)
  # one smoothing parameter per margin: that is what anisotropic means
  expect_identical(pen@params, c("lambda1", "lambda2"))
  P <- penalties7::penalty_matrix(pen, list(lambda1 = 1, lambda2 = 1))
  expect_identical(dim(P), c(15L, 15L))

  # a roughness penalty is rank deficient: its null space holds the
  # surfaces of no curvature in either direction
  r <- penalties7::penalty_rank(pen)
  expect_lt(r, 15L)
  expect_gt(r, 0L)

  # and the two directions really are penalized apart
  b <- rnorm(15)
  v1 <- penalties7::penalty_value(pen, b, list(lambda1 = 1e3, lambda2 = 1e-3))
  v2 <- penalties7::penalty_value(pen, b, list(lambda1 = 1e-3, lambda2 = 1e3))
  expect_false(isTRUE(all.equal(v1, v2)))

  res <- check_term(te(x, z, k = 4), dd, verbose = FALSE)
  expect_true(all(res$status == "OK"),
              info = paste(res$check[res$status != "OK"], collapse = ", "))
})

test_that("anisotropic = FALSE sums the margins under one parameter", {
  iso <- term_build(te(x, z, k = 4, anisotropic = FALSE), dd)
  expect_identical(term_penalty(iso)@params, "lambda")
  # the isotropic matrix is the anisotropic one at equal parameters
  ani <- term_build(te(x, z, k = 4), dd)
  expect_equal(penalties7::penalty_matrix(term_penalty(iso), list(lambda = 1)),
               penalties7::penalty_matrix(term_penalty(ani),
                                          list(lambda1 = 1, lambda2 = 1)),
               ignore_attr = TRUE)
  expect_error(te(x, z, anisotropic = NA), "TRUE or FALSE")
})

test_that("a tensor smooth recovers an interaction surface", {
  set.seed(3)
  m <- 300
  d2 <- data.frame(x = runif(m), z = runif(m))
  truth <- function(x, z) sin(pi * x) * (z - 0.5)
  d2$y <- truth(d2$x, d2$z) + rnorm(m, sd = 0.05)

  built <- term_build(te(x, z, k = 5), d2)
  Z <- term_matrix(built)
  P <- penalties7::penalty_matrix(term_penalty(built),
                                  list(lambda1 = 1, lambda2 = 1))
  bhat <- solve(crossprod(Z) + 1e-3 * P + 1e-8 * diag(ncol(Z)),
                crossprod(Z, d2$y))
  expect_gt(stats::cor(as.numeric(Z %*% bhat), truth(d2$x, d2$z)), 0.98)
})

test_that("te() with by keeps one surface per level", {
  built <- term_build(te(x, z, k = 4, by = g), dd)
  expect_identical(term_npar(built), 45L)
  res <- check_term(te(x, z, k = 4, by = g), dd, verbose = FALSE)
  expect_true(all(res$status == "OK"))
})

test_that("te() is centered, so a design carrying an intercept has full rank", {
  # the tensor product contains the constant and the penalty's null space
  # contains it too, so without the constraint the design beside an intercept
  # is rank deficient by exactly one and nothing covers the deficiency
  built <- term_build(te(x, z, k = 5), dd)
  Z <- term_matrix(built)
  expect_identical(ncol(Z), 24L)
  expect_lt(max(abs(colSums(Z))), 1e-10)

  X <- cbind(1, Z)
  expect_identical(qr(X)$rank, ncol(X))
  sv <- svd(X)$d
  # the unconstrained block gives a condition number at the reciprocal of the
  # machine epsilon; this one is a design that can be solved
  expect_lt(sv[1] / sv[length(sv)], 1e6)

  # the direction removed was one of the penalty's null directions, so the
  # rank of the penalty does not move with the dimension
  marg <- lapply(list(dd$x, dd$z), function(v) {
    r <- range(v); pad <- diff(r) * 0.001 + .Machine$double.eps
    basis7::bspline_basis(lower = r[1] - pad, upper = r[2] + pad,
                          dimension = 5L, degree = 3L)
  })
  comps <- lapply(1:2, function(j) {
    Pj <- basis7::basis_gram(marg[[j]], order = 2L)
    Pj <- Pj / max(1, max(abs(Pj)))
    Reduce(kronecker, rev(lapply(1:2, function(i)
      if (i == j) Pj else diag(5L))))
  })
  expect_identical(penalties7::penalty_rank(term_penalty(built)),
                   penalties7::penalty_rank(penalties7::additive_penalty(comps)))
})

test_that("the penalized information of a centered tensor is definite", {
  # chol() is not a rank test: on the unconstrained block it succeeded or
  # failed by the luck of rounding while the smallest eigenvalue sat at the
  # rounding floor, so vcov(), confint() and the outer criterion were computed
  # on a singular matrix. The eigenvalue is what has to be asserted.
  built <- term_build(te(x, z, k = 5), dd)
  Z <- term_matrix(built)
  X <- cbind(1, Z)
  pen <- term_penalty(built)
  S <- matrix(0, ncol(X), ncol(X))
  idx <- 1L + seq_len(ncol(Z))
  S[idx, idx] <- penalties7::penalty_hessian(pen, numeric(ncol(Z)),
                                             list(lambda1 = 1, lambda2 = 1))
  ev <- eigen(crossprod(X) + S, symmetric = TRUE, only.values = TRUE)$values
  expect_gt(min(ev), sqrt(.Machine$double.eps) * max(ev))
  expect_silent(chol(crossprod(X) + S))
})

test_that("the smooths are routed by the interpreter and validated", {
  out <- interpret_formula(y ~ w + s(x) + te(x, z, k = 4), dd)
  expect_named(out$terms, c("linpar", "s(x)", "te(x, z, k = 4)"))
  expect_true(S7::S7_inherits(out$terms[["s(x)"]], SmoothTerm))

  expect_error(te(x), "at least two covariates")
  expect_error(s(x, k = 2), "at least 3")
  expect_error(s(x, k = 3, degree = 3), "must exceed")
  expect_error(s(x, label = ""), "non-empty")
  expect_error(term_build(s(nope), dd), "not found")
})

test_that("a factor by is where a smooth's block can be sparse", {
  # each row sits in the block of its own level and nowhere else, a density
  # of 1/m -- the shape .random_block() has. The basis itself is dense by
  # construction, so this is the only place a smooth admits the storage.
  set.seed(41)
  n <- 600L
  m <- 40L
  d <- data.frame(x = stats::runif(n), z = stats::runif(n),
                  g = factor(sample.int(m, n, TRUE)))

  a <- term_build(s(x, k = 8, by = g), d)
  b <- term_build(s(x, k = 8, by = g, sparse = TRUE), d)
  expect_true(is.matrix(term_matrix(a)))
  expect_true(methods::is(term_matrix(b), "sparseMatrix"))
  # the same block, only stored differently
  expect_identical(term_coef_names(b), term_coef_names(a))
  expect_equal(unname(as.matrix(term_matrix(b))), unname(term_matrix(a)))
  expect_lt(as.numeric(utils::object.size(term_matrix(b))),
            as.numeric(utils::object.size(term_matrix(a))) / 5)

  # the storage is in the blueprint, so prediction does not densify
  nd <- d[1:25, , drop = FALSE]
  expect_true(methods::is(term_predict(b, nd), "sparseMatrix"))
  expect_equal(unname(as.matrix(term_predict(b, nd))),
               unname(term_predict(a, nd)))

  # and a tensor product takes it on the same terms
  tb <- term_build(te(x, z, k = 4, by = g, sparse = TRUE), d)
  expect_true(methods::is(term_matrix(tb), "sparseMatrix"))
})

test_that("a smooth refuses sparse where there is nothing to build on", {
  set.seed(42)
  d <- data.frame(x = stats::runif(50), z = stats::runif(50),
                  g = factor(rep(1:5, 10)))
  # no `by` at all: the basis is dense by construction, the Demmler-Reinsch
  # rotation making it so
  expect_error(s(x, k = 6, sparse = TRUE), "nothing to build on")
  expect_error(te(x, z, k = 5, sparse = TRUE), "nothing to build on")
  # a NUMERIC by multiplies the basis, so the block is as dense as it is.
  # Known only at the build, which is where it is said.
  expect_error(term_build(s(x, k = 6, by = z, sparse = TRUE), d),
               "nothing to build on")
  expect_error(s(x, sparse = "yes"), "TRUE or FALSE")
  # the default is untouched
  expect_true(methods::is(term_matrix(term_build(s(x, k = 6, by = g), d)),
                          "matrix"))
})
