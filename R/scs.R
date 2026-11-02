#' @include term_classes.R generics.R nonlinear.R variants.R shape.R
#' @include block_deriv3.R block_deriv4.R
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
#' @param by The `by` expression, unevaluated, or `NULL`.
#' @param shape A [shape_constraint].
#' @param spec A named list of the settings that the build uses: the
#'   `smoother`, `default_label`, `by_hyper` and, for a located shape under a
#'   factor `by`, `index_by`, the index of each level.
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
#' monotone ([shape_increasing()], [shape_decreasing()]), of constant
#' curvature ([shape_convex()], [shape_concave()]), both
#' ([shape_increasing_convex()] and the other three combinations), or
#' unimodal or sigmoid ([shape_unimodal()], [shape_sigmoid()]). The term
#' enters the
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
#' A shape with a value at a point, `value = c(at = a, value = v)`, removes
#' the level by that value instead: the columns are taken relative to their
#' row at \eqn{a}, and the contribution is
#' \eqn{f(x) = v + (B(x) - B(a))\,A_{-1}\,b_{-1}(\gamma)}. A convex or concave
#' shape with a slope at a point, `slope = c(at = a, value = v)`, has one
#' coefficient fewer, because the slope at the left end is solved from
#' \eqn{f'(a) = v}; the line that this slope contributes is a fixed part of
#' the contribution. See [shape_increasing()] for when a value other than zero
#' changes the fit.
#'
#' # The penalty
#'
#' The penalty is the quadratic penalty of first differences of the
#' constrained coordinates of \eqn{\gamma}, under one smoothing parameter
#' `lambda` selected with the rest of the model. For a convex or concave
#' smooth these are the second divided differences, and the slope at the left
#' end stays free; for a combination the anchored slope is constrained and
#' is the first coordinate of the differences. Under a strong penalty a
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
#' # Positive and bounded functions
#'
#' The shapes do not include a positive or a bounded function. When a
#' parameter must be positive, or must lie between two values, the range is
#' set by the link of that parameter in the distribution, such as
#' [linkfunctions7::log_link()] or [linkfunctions7::bounded_link()], and [s()]
#' gives the smooth on the scale of the link. The whole predictor then lies in
#' the range at every iterate. A positive term beside an intercept would not be
#' identified, because a constant moved between the intercept and the term
#' leaves the likelihood unchanged.
#'
#' # Varying the shape by another variable
#'
#' A factor `by` gives one constrained curve per level, each with the shape
#' of `shape`: the block of level \eqn{l} is the block of one curve
#' multiplied by the indicator of that level, every level is centered over
#' all the rows, as for [s()], and the level means come from the factor
#' written in the formula. The coefficients are named
#' `<label>.<level>.g1`, `<label>.<level>.g2`, and so on. A value or a slope
#' at a point holds for every level. Under `by_hyper = "level"` each level
#' has its own smoothing parameter.
#'
#' A numeric `by` \eqn{z} gives a varying-coefficient term
#' \eqn{z\,f(x)}, with the shape on the coefficient \eqn{f}. The term is not
#' centered: the level of \eqn{f} is the coefficient of \eqn{z}, a free and
#' unpenalized coordinate named `<label>.const`, so writing \eqn{z} in the
#' formula as well makes it aliased. With a value at a point the level is
#' fixed by that value instead, and `value = c(at = a, value = 0)` sets the
#' effect of \eqn{z} to zero at \eqn{a}. The default label is `scs(x):z`.
#'
#' Under a factor `by` the shape may also differ between the levels: its
#' declared settings (a minimum slope, a location, a peak or a valley, the
#' point and the value of a value or a slope) may be given by level, and its
#' estimated ones (a location that is not given, the left slope of a
#' curvature shape) may be one for every level or one per level, through the
#' arguments `mode_by`, `inflection_by` and `left_slope_by` of the shape (see
#' [shape_increasing()] and [shape_unimodal()]). A location chosen by level
#' is searched by [statmodels7::statmod()] one level at a time from a pilot
#' fit of `s(x, by = g)` ([term_variant_stages()]).
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
#' @param by An optional factor or numeric variable, given as a bare name or
#'   an expression in the data, as in [s()].
#' @param by_hyper With a factor `by`, whether the levels share one smoothing
#'   parameter (`"shared"`, the default) or each has its own (`"level"`).
#'   Under `"level"` a held value or an `id` is qualified by the level,
#'   `hyper = c(lambda.a = 2)`.
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
  per_level <- identical(by_hyper, "level")
  if (per_level && is.null(be)) {
    stop(sprintf(paste0("'%s' sets by_hyper = \"level\" and has no 'by'. One",
                        " smoothing parameter per level\n  needs a factor 'by'."),
                 label), call. = FALSE)
  }
  # one smoothing parameter per level is named after the levels, which only
  # the build knows, so its names are checked there (as for s())
  ScsTerm(label = label, var = xe, by = be, shape = shape, sparse = FALSE,
          spec = list(smoother = smoother, default_label = default,
                      by_hyper = by_hyper),
          hyper = if (per_level) as_hyper(hyper, label)
                  else smooth_hyper(hyper, "lambda", label),
          ids = check_ids(id, if (per_level) NULL else "lambda", label),
          X = NULL, coef_names = character(0), build_info = list(),
          penalty = NULL)
}

.scs_greville <- function(b) {
  d <- b@basis_params$degree
  t <- c(rep(b@lower, d + 1L), b@basis_params$knots, rep(b@upper, d + 1L))
  vapply(seq_len(b@dimension), function(j) mean(t[(j + 1L):(j + d)]), 0)
}

# X A at any x, uncentered: inside the basis X A; beyond its ends the line
# with the end slope
.scs_raw <- function(bp, x, A) {
  lo <- bp$basis@lower; up <- bp$basis@upper
  xc <- pmin(pmax(x, lo), up)
  W <- basis7::basis_eval(bp$basis, xc) %*% A
  out <- which(x != xc)
  if (length(out)) {
    D <- basis7::basis_deriv(bp$basis, xc[out], order = 1L) %*% A
    W[out, ] <- W[out, , drop = FALSE] + (x[out] - xc[out]) * D
  }
  W
}

# the centered design of one piece, W = X A_{-1} - center: the column means at
# the build's rows, or the row at the point where a value is fixed
.scs_W <- function(bp, pc, x) sweep(.scs_raw(bp, x, pc$A1), 2L, pc$center)

# the part of a piece's contribution that no coefficient moves: the value
# fixed at a point, and the column a0 of the slope that a slope at a point
# fixes
.scs_fixed <- function(bp, pc, x) {
  out <- rep(pc$offset, length(x))
  if (!is.null(pc$a0)) out <- out + drop(.scs_raw(bp, x, matrix(pc$a0))) - pc$center0
  out
}

.scs_check_point <- function(p, what, b, label) {
  if (is.null(p)) return(invisible(NULL))
  if (p$at < b@lower || p$at > b@upper) {
    stop(sprintf(paste0(
      "the point at which '%s' fixes its %s (%s) is outside the interval of the",
      " basis,\n  [%s, %s]; the arguments 'lower' and 'upper' of the smoother",
      " widen it."), label, what, format(p$at), format(b@lower), format(b@upper)),
      call. = FALSE)
  }
  invisible(NULL)
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

# THE MULTIPLIERS OF THE LEVELS, one column per level: the indicators of a
# factor `by`, the variable of a numeric one, a column of ones without `by`.
# The contribution is sum_l G[, l] f_l(x).
.scs_G <- function(term, kind, levels, data) {
  if (identical(kind, "none")) return(matrix(1, nrow(data), 1L))
  by <- .smooth_by(term@by, data, levels = levels)
  if (!identical(by$kind, kind)) {
    stop(sprintf("the 'by' variable of '%s' was %s at the build and is not here.",
                 term@label, kind), call. = FALSE)
  }
  if (anyNA(by$value) || (identical(kind, "numeric") && !all(is.finite(by$value)))) {
    stop(sprintf("the 'by' variable of '%s' has missing or infinite values.",
                 term@label), call. = FALSE)
  }
  if (identical(kind, "factor")) {
    outer(as.integer(by$value), seq_along(levels), `==`) * 1
  } else {
    matrix(by$value)
  }
}

# The map of one shape at the Greville abscissae, with the level removed or,
# under a numeric `by` without a value at a point, kept: the coordinates, the
# constrained and penalized ones, and the coordinates of shape_start() that
# the term does not carry.
.scs_piece <- function(shape, b, xi, x, keep_level, label) {
  m <- shape_map(shape, xi)
  A1 <- m$A[, -1L, drop = FALSE]
  cons <- m$cons - 1L
  lower <- m$lower[-1L]
  # first differences of the coordinates the shape declares penalized; zero
  # on the free slope of a plain curvature shape, which the penalty does not
  # see (the anchored slope of a combination is in the walk)
  pen <- (if (is.null(m$pen)) m$cons else m$pen) - 1L
  pv <- shape@shape_params$value
  ps <- shape@shape_params[["slope"]]
  .scs_check_point(pv, "value", b, label)
  .scs_check_point(ps, "slope", b, label)
  # A SLOPE AT A POINT ELIMINATES THE FREE SLOPE. f'(a) = d(a)' A_{-1} b_{-1}
  # (the level drops out, the derivative of the basis summing to zero), and
  # the free slope b_2, column 1 of A_{-1}, enters with c_1 = 1; solving for
  # it leaves the columns A_{-1,j} - A_{-1,1} c_j / c_1 and a fixed column
  # a0 = A_{-1,1} v / c_1. A shape that constrains b_2 has no such coordinate.
  a0 <- NULL
  drop_b <- 1L
  if (!is.null(ps)) {
    if (1L %in% cons) {
      stop(sprintf(paste0("the shape of '%s' constrains every slope, so a",
                          " slope at a point cannot be fixed."), label),
           call. = FALSE)
    }
    cc <- drop(basis7::basis_deriv(b, ps$at, order = 1L) %*% A1)
    a0 <- A1[, 1L] * ps$value / cc[1L]
    A1 <- A1[, -1L, drop = FALSE] - outer(A1[, 1L], cc[-1L] / cc[1L])
    cons <- cons - 1L
    pen <- pen - 1L
    lower <- lower[-1L]
    drop_b <- 1:2
  }
  # A NUMERIC `by` KEEPS THE LEVEL. The term is z f(x), and the level of f is
  # the coefficient of z, which no other term carries (s() puts back its
  # constant column for the same reason): the level is a free coordinate in
  # front, unpenalized, and the term is not centered.
  if (keep_level) {
    A1 <- cbind(m$A[, 1L], A1)
    cons <- cons + 1L
    pen <- pen + 1L
    lower <- c(m$lower[1L], lower)
    drop_b <- setdiff(drop_b, 1L)
  }
  pc <- list(shape = shape, A1 = A1, a0 = a0, cons = cons, lower = lower,
             pen = pen, drop_b = drop_b, center0 = 0, offset = 0)
  # THE LEVEL IS REMOVED BY CENTERING OVER THE DATA, OR BY A VALUE AT A POINT:
  # f(x) = v + (B(x) - B(a)) A_{-1} b_{-1}, the level entering f(a) with
  # coefficient 1 because the basis sums to one. Under a factor `by` every
  # level is centered over all the rows, as s() does.
  full <- cbind(if (is.null(a0)) NULL else a0, A1)
  bp <- list(basis = b)
  cen <- if (keep_level) numeric(ncol(full))
         else if (is.null(pv)) colMeans(.scs_raw(bp, x, full))
         else drop(.scs_raw(bp, pv$at, full))
  if (!is.null(a0)) {
    pc$center0 <- cen[1L]
    cen <- cen[-1L]
  }
  pc$center <- cen
  if (!is.null(pv)) pc$offset <- pv$value
  b0 <- shape_start(shape, x, NULL, xi)
  pc$coef0 <- .shape_gamma(shape@link, b0[setdiff(seq_along(b0), drop_b)],
                           cons, lower)
  P <- matrix(0, ncol(A1), ncol(A1))
  if (length(pen) > 1L) {
    D <- diff(diag(length(pen)))
    P[pen, pen] <- crossprod(D)
  }
  pc$P <- P
  pc
}

# the coefficient vector from one column per level: the tied coordinate is
# the mean of the levels' values
.scs_collapse <- function(bp, Gm) {
  if (is.null(bp$C)) return(as.vector(Gm))
  c(mean(Gm[1L, ]), as.vector(Gm[-1L, , drop = FALSE]))
}

# b(gamma) and its derivatives level by level, one column per level
.scs_B <- function(bp, coef, order) {
  Gm <- matrix(coef[bp$ix], bp$nc, bp$m)
  out <- vapply(seq_len(bp$m), function(l) {
    pc <- bp$pieces[[bp$piece_of[l]]]
    .shape_b(bp$link, Gm[, l], pc$cons, pc$lower, order)
  }, numeric(bp$nc))
  matrix(out, bp$nc, bp$m)
}

# the contribution, sum_l G[, l] (W_l b_l + fixed_l)
.scs_value <- function(bp, Ws, fx, G, B0) {
  out <- 0
  for (l in seq_len(bp$m)) {
    k <- bp$piece_of[l]
    out <- out + G[, l] * (drop(Ws[[k]] %*% B0[, l]) + fx[[k]])
  }
  drop(out)
}

# a block whose columns of level l are G[, l] W_l diag(D[, l]), collapsed
# onto the coefficients where a coordinate is tied across the levels
.scs_block <- function(bp, Ws, G, D) {
  n <- nrow(G)
  X <- do.call(cbind, lapply(seq_len(bp$m), function(l) {
    Ws[[bp$piece_of[l]]] * (G[, l] * rep(D[, l], each = n))
  }))
  if (is.null(bp$C)) X else X %*% bp$C
}

# the designs and fixed parts at new rows, one per piece
.scs_at <- function(term, bp, data) {
  x <- .scs_x(term, data)
  list(W = lapply(bp$pieces, function(pc) .scs_W(bp, pc, x)),
       fixed = lapply(bp$pieces, function(pc) .scs_fixed(bp, pc, x)),
       G = .scs_G(term, bp$by_kind, bp$by_levels, data), x = x)
}

S7::method(term_build, ScsTerm) <- function(term, data, ...) {
  x <- .scs_x(term, data)
  sm <- term@spec$smoother
  b <- basis7::smoother_basis(sm, x)
  xi <- .scs_greville(b)
  by <- .smooth_by(term@by, data)
  kind <- if (is.null(by)) "none" else by$kind
  lev <- if (identical(kind, "factor")) by$levels else NULL
  if (identical(term@spec$by_hyper, "level") && !identical(kind, "factor")) {
    stop(sprintf(paste0("'%s' sets by_hyper = \"level\" and has no factor 'by'.",
                        " One smoothing\n  parameter per level needs levels to have",
                        " one each."), term@label), call. = FALSE)
  }
  G <- .scs_G(term, kind, lev, data)
  m <- ncol(G)
  if (identical(kind, "numeric") && isTRUE(term@spec$default_label)) {
    term@label <- paste0(term@label, ":", paste(deparse(term@by), collapse = ""))
  }
  sp <- term@shape@shape_params
  named <- .shape_level_names(term@shape)
  tie <- identical(sp$left_slope_by, "shared")
  if (!identical(kind, "factor")) {
    if (length(named)) {
      stop(sprintf(paste0("'%s' gives settings of its shape by level, which need a",
                          " factor 'by'."), term@label), call. = FALSE)
    }
    if (identical(sp$location_by, "shared") || tie) {
      stop(sprintf(paste0("'%s' shares a quantity of its shape across the levels",
                          " (\"shared\" in mode_by,\n  inflection_by or",
                          " left_slope_by), which needs a factor 'by'."), term@label),
           call. = FALSE)
    }
  } else if (length(setdiff(named, lev))) {
    stop(sprintf("'%s' names %s in the settings of its shape, which are not levels of '%s'.",
                 term@label, paste0("'", setdiff(named, lev), "'", collapse = ", "),
                 paste(deparse(term@by), collapse = "")), call. = FALSE)
  }
  # THE SHAPE OF EACH LEVEL: the settings given by level resolved for it and,
  # for a located shape with no location searched level by level, the index
  # that term_variants() set (term@spec$index_by). Levels with the same shape
  # share one piece.
  idx_by <- term@spec$index_by
  shapes <- if (identical(kind, "factor")) {
    lapply(lev, function(l) {
      sh <- .shape_level(term@shape, l)
      if (.scs_auto_location(term) && !identical(sp$location_by, "shared") &&
          !is.null(idx_by)) {
        sh@shape_params$index <- as.integer(idx_by[[l]])
      }
      .scs_resolve_shape(sh, xi)
    })
  } else list(.scs_resolve_shape(term@shape, xi))
  key <- vapply(shapes, function(sh)
    paste(class(sh)[1L], paste(deparse(sh@shape_params), collapse = "")), "")
  ukey <- unique(key)
  piece_of <- match(key, ukey)
  if (!identical(kind, "factor")) piece_of <- rep(1L, m)
  keep_level <- identical(kind, "numeric") && is.null(sp$value)
  if (tie && any(vapply(shapes, function(sh) !is.null(sh@shape_params[["slope"]]), logical(1)))) {
    stop(sprintf(paste0("'%s' fixes a slope at a point, which removes the left",
                        " slope that\n  left_slope_by = \"shared\" would tie across",
                        " the levels."), term@label), call. = FALSE)
  }
  pieces <- lapply(ukey, function(k) .scs_piece(shapes[[match(k, key)]], b, xi, x,
                                                keep_level, term@label))
  nc <- unique(vapply(pieces, function(pc) ncol(pc$A1), integer(1)))
  if (length(nc) != 1L) {
    stop(sprintf("the levels of '%s' have shapes with different numbers of coefficients.",
                 term@label), call. = FALSE)
  }
  # THE COEFFICIENTS OF EACH LEVEL are column l of matrix(coef[ix], nc, m).
  # Without a tie they are the levels' blocks one after another; with
  # left_slope_by = "shared" the first coordinate (the left slope of a
  # curvature shape) is one coefficient for every level, and the block of
  # the levels' columns is collapsed onto the coefficients by the 0/1 matrix C.
  if (tie) {
    ix <- rbind(1L, matrix(1L + seq_len((nc - 1L) * m), nc - 1L, m))
    ncoef <- 1L + (nc - 1L) * m
    C <- matrix(0, nc * m, ncoef)
    C[cbind(seq_len(nc * m), as.vector(ix))] <- 1
  } else {
    ix <- matrix(seq_len(nc * m), nc, m)
    ncoef <- nc * m
    C <- NULL
  }
  bp <- list(basis = b, xi = xi, pieces = pieces, piece_of = piece_of,
             link = term@shape@link, nc = nc, m = m, ncoef = ncoef, ix = ix, C = C,
             n = length(x), var = term@var, x = x, shape = pieces[[1L]]$shape,
             by = term@by, by_kind = kind, by_levels = lev, keep_level = keep_level,
             G = G)
  bp$W <- lapply(pieces, function(pc) .scs_W(bp, pc, x))
  bp$fixed <- lapply(pieces, function(pc) .scs_fixed(bp, pc, x))
  G0 <- matrix(unlist(lapply(piece_of, function(k) pieces[[k]]$coef0), use.names = FALSE),
               nc, m)
  coef <- .scs_collapse(bp, G0)
  bp$coef0 <- coef
  gn <- if (keep_level) c("const", sprintf("g%d", seq_len(nc - 1L)))
        else sprintf("g%d", seq_len(nc))
  lab <- function(g) as.character(t(outer(lev, g, function(a, h)
    sprintf("%s.%s.%s", term@label, a, h))))
  term@coef_names <- if (tie) c(sprintf("%s.%s", term@label, gn[1L]), lab(gn[-1L]))
  else if (identical(kind, "factor")) lab(gn)
  else sprintf("%s.%s", term@label, gn)
  Ps <- lapply(piece_of, function(k) pieces[[k]]$P)
  bp$penalties <- NULL
  if (identical(term@spec$by_hyper, "level")) {
    # ONE PENALTY PER LEVEL, under its own smoothing parameter, as s() under
    # by_hyper = "level": the names are checked here, where the levels exist.
    # A tied coordinate sits in every level's entry, and the entries' log
    # pseudo-determinants add up correctly only if none of them penalizes it.
    if (tie && any(vapply(Ps, function(P) any(P[1L, ] != 0), logical(1)))) {
      stop(sprintf(paste0("'%s' ties a penalized left slope across the levels,",
                          " and by_hyper = \"level\"\n  would penalize it once per",
                          " level; use by_hyper = \"shared\"."), term@label),
           call. = FALSE)
    }
    parts <- lapply(seq_len(m), function(l) {
      .penalty_entry(lev[l], ix[, l], penalties7::quadratic_penalty(Ps[[l]]))
    })
    parts <- .entry_hyper(parts, term@hyper, term@ids, term@label,
                          "this smooth's penalty")
    term@hyper <- do.call(c, c(list(list()), lapply(parts, function(en) {
      if (!length(en$fixed) || !nzchar(en$name)) return(en$fixed)
      stats::setNames(en$fixed, paste0(names(en$fixed), ".", en$name))
    })))
    bp$penalties <- parts
    term@penalty <- NULL
  } else if (tie) {
    S <- crossprod(C, as.matrix(Matrix::bdiag(Ps)) %*% C)
    term@penalty <- penalties7::quadratic_penalty(S)
  } else if (length(unique(Ps)) == 1L) {
    term@penalty <- penalties7::quadratic_penalty(Ps[[1L]], blocks = m)
  } else {
    term@penalty <- penalties7::quadratic_penalty(as.matrix(Matrix::bdiag(Ps)))
  }
  term@build_info <- bp
  term@X <- matrix(0, length(x), ncoef)
  term_refresh(term, coef)
}

S7::method(term_penalties, ScsTerm) <- function(term, ...) {
  ent <- term@build_info$penalties
  if (is.null(ent)) return(S7::method(term_penalties, model_term)(term, ...))
  ent
}

S7::method(term_refresh, ScsTerm) <- function(term, coef, ...) {
  .assert_built(term)
  bp <- term@build_info
  coef <- .scs_coef(bp, coef)
  bp$coef <- coef
  bp$value <- .scs_value(bp, bp$W, bp$fixed, bp$G, .scs_B(bp, coef, 0L))
  X <- .scs_block(bp, bp$W, bp$G, .scs_B(bp, coef, 1L))
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
  B0 <- .scs_B(bp, coef, 0L)
  if (is.null(newdata)) return(.scs_value(bp, bp$W, bp$fixed, bp$G, B0))
  at <- .scs_at(term, bp, newdata)
  .scs_value(bp, at$W, at$fixed, at$G, B0)
}

S7::method(term_predict, ScsTerm) <- function(term, newdata, ...) {
  .assert_built(term)
  bp <- term@build_info
  at <- .scs_at(term, bp, newdata)
  X <- .scs_block(bp, at$W, at$G, .scs_B(bp, bp$coef, 1L))
  colnames(X) <- term@coef_names
  X
}

# X_ij = G_il W_ij b'_j(gamma_j) for column j of level l: the derivative in
# gamma_c touches column c only, so the three contractions are diagonal
S7::method(term_block_deriv, ScsTerm) <- function(term, coef = NULL, v, ...) {
  .assert_built(term)
  bp <- term@build_info
  coef <- .scs_coef(bp, coef)
  v <- as.numeric(v)
  if (length(v) != bp$ncoef) {
    stop(sprintf("'v' must have length %d, the term's coefficients.", bp$ncoef),
         call. = FALSE)
  }
  .scs_block(bp, bp$W, bp$G, .scs_B(bp, coef, 2L) * matrix(v[bp$ix], bp$nc, bp$m))
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
  B2 <- .scs_B(bp, coef, 2L)
  out <- unlist(lapply(seq_len(bp$m), function(l) {
    B2[, l] * colSums(A[, bp$ix[, l], drop = FALSE] * bp$W[[bp$piece_of[l]]] * bp$G[, l])
  }), use.names = FALSE)
  if (is.null(bp$C)) out else drop(crossprod(bp$C, out))
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
  .scs_block(bp, bp$W, bp$G, .scs_B(bp, coef, 3L) * matrix((v * u)[bp$ix], bp$nc, bp$m))
}

# one order more: the k-th derivative of column c is G W_c b^(k+1)(gamma_c),
# so the contraction along k directions multiplies them coordinate by
# coordinate
.scs_block_deriv_k <- function(term, coef, dirs) {
  .assert_built(term)
  bp <- term@build_info
  coef <- .scs_coef(bp, coef)
  dirs <- lapply(dirs, as.numeric)
  if (any(lengths(dirs) != bp$ncoef)) {
    stop(sprintf("'%s' must have length %d, the term's coefficients.",
                 paste(names(dirs), collapse = "', '"), bp$ncoef), call. = FALSE)
  }
  w <- Reduce(`*`, dirs)
  .scs_block(bp, bp$W, bp$G,
             .scs_B(bp, coef, length(dirs) + 1L) * matrix(w[bp$ix], bp$nc, bp$m))
}

S7::method(term_block_deriv3, ScsTerm) <- function(term, coef = NULL, v, u, s, ...) {
  .scs_block_deriv_k(term, coef, list(v = v, u = u, s = s))
}

S7::method(term_block_deriv4, ScsTerm) <- function(term, coef = NULL, v, u, s, r, ...) {
  .scs_block_deriv_k(term, coef, list(v = v, u = u, s = s, r = r))
}

# the start from a target: each level of a factor `by` from its own rows; a
# numeric `by` from the constant start, with the level of f at the
# least-squares coefficient of the target on the variable
S7::method(term_coef_start, ScsTerm) <- function(term, target = NULL, ...) {
  .assert_built(term)
  bp <- term@build_info
  if (is.null(target)) return(bp$coef0)
  target <- as.numeric(target)
  if (length(target) != bp$n) return(bp$coef0)
  out <- matrix(bp$coef0[bp$ix], bp$nc, bp$m)
  if (identical(bp$by_kind, "numeric")) {
    if (bp$keep_level) {
      z <- bp$G[, 1L]
      cf <- tryCatch(stats::coef(stats::lm.fit(cbind(1, z), target))[[2L]],
                     error = function(e) NA_real_)
      if (is.finite(cf)) out[1L, 1L] <- cf
    }
    return(as.vector(out))
  }
  for (l in seq_len(bp$m)) {
    pc <- bp$pieces[[bp$piece_of[l]]]
    rows <- bp$G[, l] != 0
    if (sum(rows) < 2L) next
    b0 <- tryCatch(shape_start(pc$shape, bp$x[rows], target[rows], bp$xi,
                               rule = term@spec$start_rule),
                   error = function(e) NULL)
    if (is.null(b0)) next
    g <- .shape_gamma(bp$link, b0[setdiff(seq_along(b0), pc$drop_b)], pc$cons, pc$lower)
    if (all(is.finite(g))) out[, l] <- g
  }
  .scs_collapse(bp, out)
}

# THE INDEX OF A LOCATED SHAPE IS SETTLED AT THE BUILD, where the Greville
# abscissae exist: the index the shape carries, else the abscissa nearest to
# its location (mode or inflection), else the middle one, which a fitting
# layer replaces through term_variants()
.scs_resolve_shape <- function(shape, xi) {
  if (!S7::S7_inherits(shape, ShapeLocated)) return(shape)
  p <- shape@shape_params
  K <- length(xi)
  rng <- .shape_index_range(shape, K)
  if (!is.null(p$index)) {
    if (!p$index %in% rng) {
      stop(sprintf("the index of the location (%d) is outside %d to %d.",
                   p$index, min(rng), K), call. = FALSE)
    }
    return(shape)
  }
  M <- if (!is.null(p$at)) rng[which.min(abs(xi[rng] - p$at))] else rng[(length(rng) + 1L) %/% 2L]
  shape@shape_params$index <- as.integer(M)
  shape
}

.scs_auto_location <- function(term) {
  S7::S7_inherits(term@shape, ShapeLocated) &&
    is.null(term@shape@shape_params$index) && is.null(term@shape@shape_params$at)
}

# where the pilot puts the location: its extremum for a unimodal shape, its
# steepest point for a sigmoid, read at the build's rows
.scs_pilot_center <- function(shape, x, pilot) {
  sg <- shape@shape_params$sign
  if (S7::S7_inherits(shape, ShapeSigmoid)) {
    o <- order(x); xs <- x[o]; ps <- pilot[o]
    keep <- diff(xs) > 0
    if (!any(keep)) return(NA_real_)
    sl <- (diff(ps) / diff(xs))[keep]
    mid <- ((xs[-1L] + xs[-length(xs)]) / 2)[keep]
    return(mid[which.max(sg * sl)])
  }
  x[which.max(sg * pilot)]
}

# A located shape whose location is GIVEN is searched too, over its two
# starts: the stand-in is the same, and term_variants() returns the two
.scs_two_starts <- function(term) {
  S7::S7_inherits(term@shape, ShapeLocated) && !.scs_auto_location(term) &&
    is.null(term@spec$start_rule) &&
    !(term_is_built(term) && identical(term@build_info$by_kind, "numeric"))
}

S7::method(term_pilot, ScsTerm) <- function(term, ...) {
  if (!.scs_auto_location(term) && !.scs_two_starts(term)) return(NULL)
  # an unconstrained smooth on the same basis and with the same `by`, under
  # the term's own label
  args <- list(s, term@var, smoother = term@spec$smoother, label = term@label)
  if (!is.null(term@by)) args$by <- term@by
  eval(as.call(args))
}

# THE INDEX NEAREST TO THE PILOT'S LOCATION, level by level. A numeric `by`
# reads f = pilot / z on the rows where z is not zero, the pilot being z f(x).
.scs_pilot_index <- function(term, bp, pilot) {
  K <- length(bp$xi)
  rng <- .shape_index_range(term@shape, K)
  mid <- rng[(length(rng) + 1L) %/% 2L]
  ok <- !is.null(pilot) && length(pilot) == bp$n && any(is.finite(pilot))
  near <- function(rows, f, shape = term@shape) {
    if (!ok || !any(rows)) return(mid)
    x0 <- .scs_pilot_center(shape, bp$x[rows], f[rows])
    if (is.na(x0)) mid else rng[which.min(abs(bp$xi[rng] - x0))]
  }
  if (identical(bp$by_kind, "factor")) {
    return(stats::setNames(vapply(seq_len(bp$m), function(l) {
      near(bp$G[, l] != 0, pilot, .shape_level(term@shape, bp$by_levels[l]))
    }, integer(1)), bp$by_levels))
  }
  if (identical(bp$by_kind, "numeric")) {
    z <- bp$G[, 1L]
    rows <- abs(z) > 1e-8 * max(abs(z))
    return(near(rows, if (ok) pilot / ifelse(rows, z, 1) else pilot))
  }
  near(rep(TRUE, bp$n), pilot)
}

S7::method(term_variants, ScsTerm) <- function(term, pilot = NULL, ...) {
  # THE TWO STARTS OF A GIVEN LOCATION. Away from the mode of the data the
  # first start is nearly flat and the fit stays on the flat curve; the
  # spline's slopes leave it, and neither start dominates, so both are
  # fitted and the criterion keeps the better (misure_scs M15: 82 fits of
  # 960 below the better of the two starts with the first alone, 0 with both)
  # ... and from two starting smoothing parameters, 1 and 0.01: far from the
  # mode both starts stay on the flat curve when the search starts at 1
  # (misure_scs M17: 0.01 better in 69 fits of 960, worse in 11)
  if (.scs_two_starts(term)) {
    grid <- expand.grid(rule = c("sides", "spline"), lambda = c(NA, 0.01),
                        stringsAsFactors = FALSE)
    out <- lapply(seq_len(nrow(grid)), function(i) {
      v <- term
      v@spec$start_rule <- grid$rule[i]
      if (!is.na(grid$lambda[i])) v@spec$lambda_start <- grid$lambda[i]
      v@X <- NULL; v@coef_names <- character(0); v@build_info <- list(); v@penalty <- NULL
      v
    })
    names(out) <- ifelse(is.na(grid$lambda), sprintf("start: %s", grid$rule),
                         sprintf("start: %s, lambda %g", grid$rule, grid$lambda))
    return(out)
  }
  if (!.scs_auto_location(term)) return(list(term))
  .assert_built(term)
  bp <- term@build_info
  K <- length(bp$xi)
  rng <- .shape_index_range(term@shape, K)
  M0 <- .scs_pilot_index(term, bp, pilot)
  reset <- function(v) {
    v@X <- NULL; v@coef_names <- character(0); v@build_info <- list(); v@penalty <- NULL
    v
  }
  # ONE LOCATION FOR EVERY LEVEL: the versions around the median of the
  # levels' pilot indices
  shared <- identical(term@shape@shape_params$location_by, "shared")
  if (identical(bp$by_kind, "factor") && shared) {
    M0 <- as.integer(round(stats::median(M0)))
  }
  if (!identical(bp$by_kind, "factor") || shared) {
    Ms <- intersect(M0 + (-2L):2L, rng)
    out <- lapply(Ms, function(M) {
      v <- reset(term)
      sh <- term@shape
      sh@shape_params$index <- as.integer(M)
      v@shape <- sh
      v
    })
    names(out) <- sprintf("M = %d", Ms)
    return(out)
  }
  # ONE LEVEL AT A TIME UNDER A FACTOR `by`. Stage k moves the index of
  # level k around its pilot index, the other levels at the indices of
  # `current` (the version chosen at the stage before) or, at the first
  # stage, at their pilot indices. The index that `current` already holds
  # for level k was fitted at the stage before and is left out.
  dots <- list(...)
  stage <- if (is.null(dots$stage)) 1L else as.integer(dots$stage)
  current <- dots$current
  idx <- if (!is.null(current) && !is.null(current@spec$index_by))
    current@spec$index_by else M0
  Ms <- intersect(M0[stage] + (-2L):2L, rng)
  if (!is.null(current)) Ms <- setdiff(Ms, idx[stage])
  out <- lapply(Ms, function(M) {
    v <- reset(term)
    v@spec$index_by <- replace(idx, stage, as.integer(M))
    v
  })
  names(out) <- sprintf("%s: M = %d", bp$by_levels[stage], Ms)
  out
}

S7::method(term_variant_stages, ScsTerm) <- function(term, ...) {
  if (!.scs_auto_location(term) || !term_is_built(term) ||
      !identical(term@build_info$by_kind, "factor") ||
      identical(term@shape@shape_params$location_by, "shared")) return(1L)
  term@build_info$m
}

# the covariate at the build's rows, recovered from the stored design: the
# build keeps the basis and the centered design, not the data
.scs_x_built <- function(term) term@build_info$x

S7::method(print, ScsTerm) <- function(x, ...) {
  if (term_is_built(x)) {
    bp <- x@build_info
    by <- if (identical(bp$by_kind, "factor")) sprintf(", by a factor of %d levels", bp$m)
          else if (identical(bp$by_kind, "numeric")) ", by a numeric variable"
          else ""
    cat(sprintf("<ScsTerm> '%s' built: %d coefficient%s; shape %s%s\n",
                x@label, ncol(x@X), if (ncol(x@X) == 1L) "" else "s",
                x@shape@shape_name, by))
  } else {
    cat(sprintf("<ScsTerm> '%s' (specification; call term_build() with data)\n",
                x@label))
  }
  invisible(x)
}

S7::method(term_hyper_start, ScsTerm) <- function(term, ...) {
  if (is.null(term@spec$lambda_start)) NULL else c(lambda = term@spec$lambda_start)
}
