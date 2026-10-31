#' @include term_classes.R generics.R nonlinear.R shape.R
NULL

#' @title S7 Class for Shape-Constrained Smooth Terms
#' @name ScsTerm
#'
#' @description
#' `ScsTerm` is the subclass of [additive_term()] for a smooth of one
#' covariate constrained to a shape: [scs()] constructs it. The B-spline
#' coefficients are written through the map of a [shape_constraint],
#' \eqn{\beta = A\,b(\gamma)}, so the contribution to the predictor,
#' \eqn{f = X_c A_{-1}\,b_{-1}(\gamma)}, is nonlinear in the coefficients
#' \eqn{\gamma} and the block is the Jacobian at the current coefficients, as
#' for [nl()].
#'
#' @inheritParams additive_term
#' @inheritParams model_term
#' @param var The covariate expression, unevaluated.
#' @param by `NULL`: a `by` variable is not yet supported.
#' @param shape A [shape_constraint].
#' @param spec A named list of the settings that the build uses: the
#'   `smoother` and `default_label`.
#' @param sparse `NULL` or `FALSE`.
#'
#' @return An S7 object of class `ScsTerm`, inheriting from [additive_term()]
#'   and [model_term()].
#'
#' @seealso [scs()], the constructor.
#'
#' @examples
#' S7::S7_inherits(scs(x, shape_increasing()), ScsTerm)
#'
#' @export
ScsTerm <- S7::new_class(
  name = "ScsTerm",
  parent = additive_term,
  properties = list(
    var = S7::class_any,
    by = S7::class_any,
    shape = S7::class_any,
    spec = S7::class_list,
    sparse = S7::class_any
  )
)

#' Shape-Constrained Smooth of One Covariate
#'
#' @description
#' Constructs a smooth function of one covariate constrained to a shape:
#' monotone ([shape_increasing()], [shape_decreasing()]) or of constant
#' curvature ([shape_convex()], [shape_concave()]). The term enters the
#' formula of any distribution parameter, like [s()].
#'
#' @details
#' # The construction
#'
#' The covariate is expanded in the B-spline basis of `smoother`, whose
#' coefficients are written as \eqn{\beta = A\,b(\gamma)} through the map of
#' the shape (see [shape_constraint]): cumulative sums of slopes, or of second
#' divided differences, on the Greville abscissae of the basis, each
#' constrained one through the inverse of the shape's link. Every value of the
#' free coefficients \eqn{\gamma} gives a curve of the shape, so the constraint
#' holds exactly at every iterate of a fit.
#'
#' The level is removed against the intercept of the model: the first column of
#' \eqn{A} is dropped and the remaining columns of \eqn{X A} are centered over
#' the data on which the term is built, so the smooth averages to zero there.
#' The coefficients are therefore \eqn{\gamma_2, \dots, \gamma_K}, named
#' `<label>.g1`, `<label>.g2`, and so on.
#'
#' # The penalty
#'
#' The penalty is the quadratic penalty of first differences of the
#' constrained coordinates of \eqn{\gamma}, under one smoothing parameter
#' `lambda` selected with the rest of the model. Under a strong penalty a
#' monotone smooth approaches a straight line, and a convex or concave one the
#' curve whose coefficients have constant second divided differences, which is
#' a parabola in the interior and has a smaller curvature in the end panels of
#' the basis, where the knots are repeated. The penalty of the smoother is not
#' used: the smoother supplies
#' the basis only.
#'
#' # The block
#'
#' The contribution is nonlinear in \eqn{\gamma}, so the block is the
#' Jacobian \eqn{X_c A_{-1}\,\mathrm{diag}(b'(\gamma))} at the current
#' coefficients, refreshed by [term_refresh()], and [term_value()] gives the
#' contribution itself. Because \eqn{b} acts on one coordinate at a time, the
#' derivatives of the block are diagonal in \eqn{\gamma}.
#'
#' # Prediction outside the data
#'
#' Beyond the ends of the basis the smooth is continued by the straight line
#' with the slope at the end, which preserves the shape.
#'
#' # Inference
#'
#' Pointwise intervals from the curvature at the mode cover somewhat less than
#' their nominal level where the curve is flat (about 93 per cent for 95 per
#' cent in simulation); [statmodels7::statmod_bootstrap()] gives a reading that does not
#' rest on that curvature.
#'
#' @param x The covariate, a bare name or an expression in the data.
#' @param shape A [shape_constraint], by default [shape_increasing()].
#' @param smoother A \pkg{basis7} smoother whose basis is a B-spline:
#'   [basis7::bspline_smooth()] (the default) or [basis7::pspline_smooth()].
#'   Only its basis is used.
#' @param by Not yet supported; must be `NULL`.
#' @param by_hyper Not yet used.
#' @param hyper `NULL` to estimate the smoothing parameter, or a value at
#'   which to hold it, as for [s()].
#' @param id `NULL`, or a string shared with other terms whose smoothing
#'   parameter is the same, as for [s()].
#' @param label The term's label; `NULL` builds `scs(<x>)` from the call.
#' @param sparse `NULL` or `FALSE`: the block is dense.
#'
#' @return An unbuilt [ScsTerm].
#'
#' @references
#' Pya, N. and Wood, S. N. (2015). Shape constrained additive models.
#' *Statistics and Computing*, 25, 543--559.
#'
#' @seealso [shape_increasing()] and the other shapes, [s()] for the
#'   unconstrained smooth.
#'
#' @examples
#' set.seed(1)
#' dd <- data.frame(x = runif(100))
#' tm <- term_build(scs(x, shape_increasing(), basis7::bspline_smooth(k = 8)), dd)
#' term_coef_names(tm)
#'
#' # Every coefficient vector gives an increasing curve.
#' f <- term_value(tm, coef = rnorm(7))
#' all(diff(f[order(dd$x)]) >= -1e-12)
#'
#' @export
scs <- function(x, shape = shape_increasing(), smoother = basis7::bspline_smooth(),
                by = NULL, by_hyper = c("shared", "level"), hyper = NULL, id = NULL,
                label = NULL, sparse = NULL) {
  xe <- substitute(x)
  be <- substitute(by)
  by_hyper <- match.arg(by_hyper)
  if (!S7::S7_inherits(shape, shape_constraint)) {
    stop("'shape' must be a shape constraint, such as shape_increasing().",
         call. = FALSE)
  }
  if (!(S7::S7_inherits(smoother, basis7::BsplineSmoother) ||
        S7::S7_inherits(smoother, basis7::PsplineSmoother))) {
    stop(paste0("'smoother' must be bspline_smooth() or pspline_smooth(): a",
                " shape on the\n  coefficients is a shape of the curve only",
                " for a B-spline basis."), call. = FALSE)
  }
  if (!is.null(be)) {
    stop("'by' is not yet supported for scs(); leave it NULL.", call. = FALSE)
  }
  if (isTRUE(sparse)) {
    stop("'sparse' must be NULL or FALSE: the block of scs() is dense.",
         call. = FALSE)
  }
  default <- is.null(label)
  if (default) label <- sprintf("scs(%s)", paste(deparse(xe), collapse = ""))
  if (!is.character(label) || length(label) != 1L || is.na(label) ||
      !nzchar(label)) {
    stop("'label' must be a single non-empty character string.", call. = FALSE)
  }
  ScsTerm(label = label, var = xe, by = NULL, shape = shape, sparse = FALSE,
          spec = list(smoother = smoother, default_label = default),
          hyper = smooth_hyper(hyper, "lambda", label),
          ids = check_ids(id, "lambda", label),
          X = NULL, coef_names = character(0), build_info = list(),
          penalty = NULL)
}

.scs_greville <- function(b) {
  d <- b@basis_params$degree
  t <- c(rep(b@lower, d + 1L), b@basis_params$knots, rep(b@upper, d + 1L))
  vapply(seq_len(b@dimension), function(j) mean(t[(j + 1L):(j + d)]), 0)
}

# the centered design W = X A_{-1} - colMeans at the build's rows, at any x:
# inside the basis X A_{-1}; beyond its ends the line with the end slope
.scs_W <- function(bp, x) {
  lo <- bp$basis@lower; up <- bp$basis@upper
  xc <- pmin(pmax(x, lo), up)
  W <- basis7::basis_eval(bp$basis, xc) %*% bp$A1
  out <- which(x != xc)
  if (length(out)) {
    D <- basis7::basis_deriv(bp$basis, xc[out], order = 1L) %*% bp$A1
    W[out, ] <- W[out, , drop = FALSE] + (x[out] - xc[out]) * D
  }
  sweep(W, 2L, bp$center)
}

.scs_coef <- function(bp, coef) {
  if (is.null(coef)) return(bp$coef)
  coef <- as.numeric(coef)
  if (length(coef) != bp$ncoef) {
    stop(sprintf("'coef' must have length %d, the term's coefficients.",
                 bp$ncoef), call. = FALSE)
  }
  coef
}

.scs_x <- function(term, data) {
  x <- as.numeric(eval(term@var, data, baseenv()))
  if (length(x) != nrow(data)) {
    stop(sprintf("the covariate of '%s' must give one value per row.",
                 term@label), call. = FALSE)
  }
  if (anyNA(x) || !all(is.finite(x))) {
    stop(sprintf("the covariate of '%s' has missing or infinite values.",
                 term@label), call. = FALSE)
  }
  x
}

S7::method(term_build, ScsTerm) <- function(term, data, ...) {
  x <- .scs_x(term, data)
  sm <- term@spec$smoother
  b <- basis7::smoother_basis(sm, x)
  xi <- .scs_greville(b)
  m <- shape_map(term@shape, xi)
  K <- length(xi)
  A1 <- m$A[, -1L, drop = FALSE]
  cons <- m$cons - 1L
  lower <- m$lower[-1L]
  bp <- list(basis = b, xi = xi, A1 = A1, center = numeric(K - 1L), cons = cons,
             lower = lower, link = term@shape@link, ncoef = K - 1L, n = length(x),
             var = term@var, x = x)
  W <- .scs_W(bp, x)
  bp$center <- colMeans(W)
  W <- sweep(W, 2L, bp$center)
  bp$W <- W
  b0 <- shape_start(term@shape, x, NULL, xi)[-1L]
  coef <- .shape_gamma(bp$link, b0, cons, lower)
  bp$coef0 <- coef
  # first differences of the constrained coordinates; zero on the free slope
  # of a curvature shape, which the penalty does not see
  P <- matrix(0, K - 1L, K - 1L)
  if (length(cons) > 1L) {
    D <- diff(diag(length(cons)))
    P[cons, cons] <- crossprod(D)
  }
  term@coef_names <- sprintf("%s.g%d", term@label, seq_len(K - 1L))
  term@penalty <- penalties7::quadratic_penalty(P)
  term@build_info <- bp
  term@X <- matrix(0, length(x), K - 1L)
  term_refresh(term, coef)
}

S7::method(term_refresh, ScsTerm) <- function(term, coef, ...) {
  .assert_built(term)
  bp <- term@build_info
  coef <- .scs_coef(bp, coef)
  b0 <- .shape_b(bp$link, coef, bp$cons, bp$lower, 0L)
  b1 <- .shape_b(bp$link, coef, bp$cons, bp$lower, 1L)
  bp$coef <- coef
  bp$value <- drop(bp$W %*% b0)
  X <- bp$W * rep(b1, each = bp$n)
  colnames(X) <- term@coef_names
  term@X <- X
  term@build_info <- bp
  term
}

S7::method(term_value, ScsTerm) <- function(term, coef = NULL, newdata = NULL, ...) {
  .assert_built(term)
  bp <- term@build_info
  if (is.null(coef) && is.null(newdata)) return(bp$value)
  coef <- .scs_coef(bp, coef)
  W <- if (is.null(newdata)) bp$W else .scs_W(bp, .scs_x(term, newdata))
  drop(W %*% .shape_b(bp$link, coef, bp$cons, bp$lower, 0L))
}

S7::method(term_predict, ScsTerm) <- function(term, newdata, ...) {
  .assert_built(term)
  bp <- term@build_info
  W <- .scs_W(bp, .scs_x(term, newdata))
  X <- W * rep(.shape_b(bp$link, bp$coef, bp$cons, bp$lower, 1L), each = nrow(W))
  colnames(X) <- term@coef_names
  X
}

# X_ij = W_ij b'_j(gamma_j): the derivative in gamma_c touches column c only,
# so the three contractions are diagonal
S7::method(term_block_deriv, ScsTerm) <- function(term, coef = NULL, v, ...) {
  .assert_built(term)
  bp <- term@build_info
  coef <- .scs_coef(bp, coef)
  v <- as.numeric(v)
  if (length(v) != bp$ncoef) {
    stop(sprintf("'v' must have length %d, the term's coefficients.", bp$ncoef),
         call. = FALSE)
  }
  bp$W * rep(.shape_b(bp$link, coef, bp$cons, bp$lower, 2L) * v, each = bp$n)
}

S7::method(term_block_contract, ScsTerm) <- function(term, coef = NULL, A, ...) {
  .assert_built(term)
  bp <- term@build_info
  coef <- .scs_coef(bp, coef)
  A <- as.matrix(A)
  if (nrow(A) != bp$n || ncol(A) != bp$ncoef) {
    stop(sprintf("'A' must be %d by %d, the observations by the term's coefficients.",
                 bp$n, bp$ncoef), call. = FALSE)
  }
  .shape_b(bp$link, coef, bp$cons, bp$lower, 2L) * colSums(A * bp$W)
}

S7::method(term_block_deriv2, ScsTerm) <- function(term, coef = NULL, v, u, ...) {
  .assert_built(term)
  bp <- term@build_info
  coef <- .scs_coef(bp, coef)
  v <- as.numeric(v); u <- as.numeric(u)
  if (length(v) != bp$ncoef || length(u) != bp$ncoef) {
    stop(sprintf("'v' and 'u' must have length %d, the term's coefficients.",
                 bp$ncoef), call. = FALSE)
  }
  bp$W * rep(.shape_b(bp$link, coef, bp$cons, bp$lower, 3L) * v * u, each = bp$n)
}

S7::method(term_coef_start, ScsTerm) <- function(term, target = NULL, ...) {
  .assert_built(term)
  bp <- term@build_info
  if (is.null(target)) return(bp$coef0)
  x <- .scs_x_built(term)
  b0 <- tryCatch(shape_start(term@shape, x, as.numeric(target), bp$xi)[-1L],
                 error = function(e) NULL)
  if (is.null(b0)) return(bp$coef0)
  g <- .shape_gamma(bp$link, b0, bp$cons, bp$lower)
  if (!all(is.finite(g))) bp$coef0 else g
}

# the covariate at the build's rows, recovered from the stored design: the
# build keeps the basis and the centered design, not the data
.scs_x_built <- function(term) term@build_info$x

S7::method(print, ScsTerm) <- function(x, ...) {
  if (term_is_built(x)) {
    cat(sprintf("<ScsTerm> '%s' built: %d coefficient%s; shape %s\n",
                x@label, ncol(x@X), if (ncol(x@X) == 1L) "" else "s",
                x@shape@shape_name))
  } else {
    cat(sprintf("<ScsTerm> '%s' (specification; call term_build() with data)\n",
                x@label))
  }
  invisible(x)
}
