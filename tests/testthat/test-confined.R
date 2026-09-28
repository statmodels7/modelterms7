# A break-point against its confinement limit. The position is clamped, so the
# contribution does not move with the coefficients that place it: the block,
# which is the Jacobian for seg() and for every smoothed construction, has a
# zero column there, and its derivatives carry the same gate. Both are checked
# against a difference of term_value() and of the block, which share no
# arithmetic with the closed forms.

cf_data <- function(n = 160, seed = 4) {
  set.seed(seed)
  dd <- data.frame(x = sort(runif(n, 0, 10)), g = factor(rep(1:4, length.out = n)))
  dd$y <- 1 + 0.3 * dd$x + 2 * (dd$x > 6) + rnorm(n, sd = 0.3)
  dd
}

confined_terms <- function(dd) {
  sm <- numericals7::smooth_probit(h = 0.4)
  list(
    seg_sharp = term_build(seg(x, psi = 5), dd),
    seg_dev = term_build(seg(x, psi ~ 0 + g), dd),
    seg_smooth = term_build(seg(x, psi = 5, smoothed = sm), dd),
    jump_smooth = term_build(jump(x, psi = 5, smoothed = sm), dd),
    jseg_smooth = term_build(jseg(x, psi = 5, smoothed = sm), dd)
  )
}

push_out <- function(b) {
  cf <- b@blueprint$coef
  # the first coefficient that places the break-point: the whole break-point
  # where it is not developed, one group's where it is, the others staying
  # inside the data
  i <- b@blueprint$index[["psi1"]][1L]
  cf[i] <- b@blueprint$lim[1L] - 2
  cf
}

num_jac <- function(f, cf) {
  cols <- lapply(seq_along(cf), function(j) {
    h <- 1e-6 * max(1, abs(cf[j]))
    up <- cf; up[j] <- cf[j] + h
    dn <- cf; dn[j] <- cf[j] - h
    (f(up) - f(dn)) / (2 * h)
  })
  cols
}

test_that("a confined break-point's block is the Jacobian of the contribution", {
  dd <- cf_data()
  for (nm in names(tl <- confined_terms(dd))) {
    b <- tl[[nm]]
    cf <- push_out(b)
    X <- as.matrix(term_matrix(term_refresh(b, cf)))
    J <- do.call(cbind, num_jac(function(p) term_value(b, p), cf))
    expect_lt(max(abs(X - J)), 1e-6 * max(1, max(abs(J))), label = nm)
    i <- b@blueprint$index[["psi1"]][1L]
    expect_equal(max(abs(X[, i])), 0, label = nm)
    # and the position did leave the data, or the test would be empty
    pos <- .seg_positions(b@blueprint, cf, nrow(dd))$psi[, 1L]
    expect_true(any(pos == b@blueprint$lim[1L]), label = nm)
    expect_gt(max(abs(J)), 0, label = nm)
  }
})

test_that("a confined break-point's block derivatives match its block", {
  dd <- cf_data()
  for (nm in names(tl <- confined_terms(dd))) {
    b <- tl[[nm]]
    cf <- push_out(b)
    m <- length(cf)
    n <- nrow(dd)
    D <- num_jac(function(p) as.matrix(term_matrix(term_refresh(b, p))), cf)
    set.seed(3)
    A <- matrix(rnorm(n * m), n, m)
    v <- rnorm(m)
    ct <- term_block_contract(b, coef = cf, A = A)
    ct_ref <- vapply(seq_len(m), function(j) sum(A * D[[j]]), numeric(1))
    expect_lt(max(abs(ct - ct_ref)), 1e-5 * max(1, max(abs(ct_ref))), label = nm)
    dv <- term_block_deriv(b, coef = cf, v = v)
    dv_ref <- Reduce(`+`, Map(function(Dj, vj) Dj * vj, D, v))
    expect_lt(max(abs(dv - dv_ref)), 1e-5 * max(1, max(abs(dv_ref))), label = nm)
  }
})
