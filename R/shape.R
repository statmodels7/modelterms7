#' @include term_classes.R generics.R
NULL

#' @title S7 Class for Shape Constraints
#' @name shape_constraint
#'
#' @description
#' `shape_constraint` is the parent class of the shapes that [scs()] imposes
#' on a smooth: [shape_increasing()], [shape_decreasing()], [shape_convex()]
#' and [shape_concave()]. A shape writes the B-spline coefficients of the
#' smooth as \eqn{\beta = A\,b(\gamma)}, where \eqn{A} is a fixed matrix of
#' cumulative sums on the Greville abscissae of the basis and \eqn{b} applies
#' the inverse of a link to the constrained coordinates one at a time, so
#' every free vector \eqn{\gamma} gives a curve of that shape.
#'
#' @details
#' # The map
#'
#' With \eqn{\xi_1 < \dots < \xi_K} the Greville abscissae of a B-spline
#' basis and \eqn{h_j = \xi_{j+1} - \xi_j}, a monotone shape parametrizes the
#' slopes \eqn{s_j = (\beta_j - \beta_{j-1}) / h_{j-1}} and a convex or
#' concave shape the second divided differences of the coefficients. A
#' B-spline whose coefficients have nonnegative slopes is nondecreasing, and
#' one whose slopes on the Greville abscissae are nondecreasing is convex,
#' which holds for any knot placement, including the knots repeated at the
#' ends that \pkg{basis7} uses.
#'
#' On a constrained coordinate \eqn{b_j = \ell_j + g^{-1}(\gamma_j)}, with
#' \eqn{g} the link of the shape (a link onto \eqn{(0, \infty)}) and
#' \eqn{\ell_j} a lower bound (zero, or the minimum slope). On a free
#' coordinate \eqn{b_j = \gamma_j}. The derivatives of the map are therefore
#' diagonal: \eqn{\partial^k b_j / \partial \gamma_j^k = (g^{-1})^{(k)}(\gamma_j)}
#' and every mixed derivative is zero.
#'
#' # A new shape
#'
#' A user-defined shape is a subclass of `shape_constraint` with methods for
#' [shape_map()] and [shape_start()].
#'
#' @param shape_name A single string naming the shape.
#' @param link A \pkg{linkfunctions7} link onto \eqn{(0, \infty)}.
#' @param shape_params A named list of the shape's own settings.
#'
#' @return An S7 object of class `shape_constraint`.
#'
#' @seealso [scs()], [shape_map()], [shape_start()].
#'
#' @examples
#' sh <- shape_increasing()
#' S7::S7_inherits(sh, shape_constraint)
#' sh@shape_name
#'
#' @export
shape_constraint <- S7::new_class(
  name = "shape_constraint",
  properties = list(
    shape_name = S7::class_character,
    link = S7::class_any,
    shape_params = S7::class_list
  )
)

#' @rdname shape_constraint
#' @export
ShapeMonotone <- S7::new_class("ShapeMonotone", parent = shape_constraint)

#' @rdname shape_constraint
#' @export
ShapeCurvature <- S7::new_class("ShapeCurvature", parent = shape_constraint)

#' @rdname shape_constraint
#' @export
ShapeLocated <- S7::new_class("ShapeLocated", parent = shape_constraint)

#' @rdname shape_constraint
#' @export
ShapeUnimodal <- S7::new_class("ShapeUnimodal", parent = ShapeLocated)

#' @rdname shape_constraint
#' @export
ShapeSigmoid <- S7::new_class("ShapeSigmoid", parent = ShapeLocated)

#' Unimodal and Sigmoid Shapes
#' @name shape_located
#'
#' @description
#' Construct the shapes of [scs()] that have a location: `shape_unimodal()`
#' rises to a mode and falls after it (or, with `valley = TRUE`, falls to a
#' minimum and rises after it); `shape_sigmoid()` increases (or, with
#' `decreasing = TRUE`, decreases) with a slope that grows up to an
#' inflection and shrinks after it.
#'
#' @details
#' # Where the location sits
#'
#' The location is an index \eqn{M} of the Greville abscissae
#' \eqn{\xi_1 < \dots < \xi_K} of the basis. For `shape_unimodal()` the slopes
#' of the coefficients are \eqn{\pm g^{-1}(\gamma_j)}, of one sign up to
#' \eqn{M} and of the other after it, so the coefficients, and with them the
#' curve, have a single extremum; \eqn{M = 1} and \eqn{M = K} are monotone
#' curves. For `shape_sigmoid()` every slope has the sign of the shape: the
#' slope at each end is constrained, and the second divided differences are
#' summed from each end toward the inflection, up to index \eqn{M} from the
#' left and down to \eqn{M + 1} from the right, so the slopes grow in absolute
#' value up to the inflection and shrink after it; \eqn{M} runs from 2 to
#' \eqn{K}.
#'
#' The penalty of [scs()] is one random walk on \eqn{\gamma} across the
#' location, whose null space is the same for every \eqn{M}, so fits at
#' different \eqn{M} are compared by their criterion.
#'
#' # How the location is chosen
#'
#' With `mode` or `inflection` given, \eqn{M} is the index of the Greville
#' abscissa nearest to it and the model is fitted once. Left at `NULL`, the
#' fitting layer chooses \eqn{M}: a pilot fit with [s()] in place of the term
#' gives the abscissa nearest to its extremum (unimodal) or to its steepest
#' point (sigmoid), \eqn{M_0}, and the model is fitted at the indices
#' \eqn{M_0 - 2, \dots, M_0 + 2} that exist; the best criterion decides (see
#' [term_pilot()] and [term_variants()]).
#'
#' @param link A \pkg{linkfunctions7} link onto \eqn{(0, \infty)}.
#' @param mode,inflection `NULL`, or a single number: the value of the
#'   covariate at which the mode or the inflection is placed.
#' @param index `NULL`, or the index \eqn{M} itself, which a fitting layer
#'   sets on the versions it tries. It takes precedence over `mode` and
#'   `inflection`.
#' @param valley `FALSE` for a maximum, `TRUE` for a minimum.
#' @param decreasing `FALSE` for an increasing sigmoid, `TRUE` for a
#'   decreasing one.
#'
#' @return An S7 object of class `ShapeUnimodal` or `ShapeSigmoid`, both
#'   inheriting from `ShapeLocated` and [shape_constraint].
#'
#' @seealso [scs()], [shape_increasing()].
#'
#' @examples
#' shape_unimodal()
#' shape_unimodal(mode = 0.4)
#' shape_unimodal(valley = TRUE)
#' shape_sigmoid(inflection = 0.5)
#' shape_sigmoid(decreasing = TRUE)
NULL

.check_location <- function(at, nm) {
  if (!is.null(at) && (!is.numeric(at) || length(at) != 1L || !is.finite(at))) {
    stop(sprintf("'%s' must be NULL or a single finite number.", nm), call. = FALSE)
  }
  at
}

.check_index <- function(index) {
  if (!is.null(index) && (!is.numeric(index) || length(index) != 1L ||
                          !is.finite(index) || index < 1 || index != round(index))) {
    stop("'index' must be NULL or a single positive whole number.", call. = FALSE)
  }
  if (is.null(index)) NULL else as.integer(index)
}

.check_flag <- function(v, nm) {
  if (!is.logical(v) || length(v) != 1L || is.na(v)) {
    stop(sprintf("'%s' must be TRUE or FALSE.", nm), call. = FALSE)
  }
  v
}

#' @rdname shape_located
#' @export
shape_unimodal <- function(link = linkfunctions7::log_link(), mode = NULL,
                           index = NULL, valley = FALSE) {
  valley <- .check_flag(valley, "valley")
  ShapeUnimodal(shape_name = if (valley) "unimodal (valley)" else "unimodal",
                link = .check_shape_link(link),
                shape_params = list(at = .check_location(mode, "mode"),
                                    index = .check_index(index),
                                    sign = if (valley) -1 else 1))
}

#' @rdname shape_located
#' @export
shape_sigmoid <- function(link = linkfunctions7::log_link(), inflection = NULL,
                          index = NULL, decreasing = FALSE) {
  decreasing <- .check_flag(decreasing, "decreasing")
  ShapeSigmoid(shape_name = if (decreasing) "sigmoid (decreasing)" else "sigmoid",
               link = .check_shape_link(link),
               shape_params = list(at = .check_location(inflection, "inflection"),
                                   index = .check_index(index),
                                   sign = if (decreasing) -1 else 1))
}

# the indices a location may take: a mode anywhere, an inflection from 2
.shape_index_range <- function(shape, K) {
  if (S7::S7_inherits(shape, ShapeSigmoid)) 2:K else seq_len(K)
}

#' Shape Constraints for a Smooth
#' @name shape_constructors
#'
#' @description
#' Construct the shapes that [scs()] imposes on a smooth.
#' `shape_increasing()` and `shape_decreasing()` give a monotone curve, with
#' a slope of at least `min_slope` in absolute value; `shape_convex()` and
#' `shape_concave()` give a curve of constant sign of curvature; the four
#' combinations `shape_increasing_convex()`, `shape_increasing_concave()`,
#' `shape_decreasing_convex()` and `shape_decreasing_concave()` give a monotone
#' curve of constant sign of curvature.
#'
#' @details
#' A monotone shape parametrizes the slopes of the B-spline coefficients on
#' the Greville abscissae, \eqn{s_j = \pm(\delta + g^{-1}(\gamma_j))} with
#' \eqn{\delta} = `min_slope`. A curvature shape keeps the slope at the left
#' end free and parametrizes the second divided differences,
#' \eqn{\pm g^{-1}(\gamma_j)}. A combination constrains the slope as well,
#' at the end where the slope is smallest in absolute value: the left end for
#' an increasing convex or a decreasing concave curve, the right end for an
#' increasing concave or a decreasing convex one, so that every slope keeps
#' its sign. See [shape_constraint] for the map.
#'
#' @param link A \pkg{linkfunctions7} link onto \eqn{(0, \infty)}, by default
#'   [linkfunctions7::log_link()]. [linkfunctions7::softplus_link()] is the
#'   other shipped choice.
#' @param min_slope A single nonnegative number, the smallest slope in
#'   absolute value.
#'
#' @return An S7 object of class `ShapeMonotone` (`shape_increasing()`,
#'   `shape_decreasing()`) or `ShapeCurvature` (the other six), both
#'   inheriting from [shape_constraint].
#'
#' @seealso [scs()], [shape_constraint].
#'
#' @examples
#' shape_increasing()
#' shape_decreasing(min_slope = 0.1)
#' shape_convex(link = linkfunctions7::softplus_link())
#' shape_increasing_concave()
#'
#' # A link onto the whole real line is rejected.
#' try(shape_concave(link = linkfunctions7::identity_link()))
NULL

.check_shape_link <- function(link) {
  if (!S7::S7_inherits(link, linkfunctions7::link)) {
    stop("'link' must be a linkfunctions7 link.", call. = FALSE)
  }
  bd <- link@link_bounds
  if (!isTRUE(bd[[1L]] == 0) || !isTRUE(is.infinite(bd[[2L]]) && bd[[2L]] > 0)) {
    stop(sprintf(paste0(
      "the link '%s' maps onto (%s, %s); a shape needs a link onto (0, Inf),",
      "\n  such as log_link() or softplus_link()."
    ), link@link_name, format(bd[[1L]]), format(bd[[2L]])), call. = FALSE)
  }
  link
}

.check_min_slope <- function(min_slope) {
  if (!is.numeric(min_slope) || length(min_slope) != 1L || !is.finite(min_slope) ||
      min_slope < 0) {
    stop("'min_slope' must be a single nonnegative finite number.", call. = FALSE)
  }
  as.numeric(min_slope)
}

#' @rdname shape_constructors
#' @export
shape_increasing <- function(link = linkfunctions7::log_link(), min_slope = 0) {
  ShapeMonotone(shape_name = "increasing", link = .check_shape_link(link),
                shape_params = list(sign = 1, min_slope = .check_min_slope(min_slope)))
}

#' @rdname shape_constructors
#' @export
shape_decreasing <- function(link = linkfunctions7::log_link(), min_slope = 0) {
  ShapeMonotone(shape_name = "decreasing", link = .check_shape_link(link),
                shape_params = list(sign = -1, min_slope = .check_min_slope(min_slope)))
}

#' @rdname shape_constructors
#' @export
shape_convex <- function(link = linkfunctions7::log_link()) {
  ShapeCurvature(shape_name = "convex", link = .check_shape_link(link),
                 shape_params = list(sign = 1))
}

#' @rdname shape_constructors
#' @export
shape_concave <- function(link = linkfunctions7::log_link()) {
  ShapeCurvature(shape_name = "concave", link = .check_shape_link(link),
                 shape_params = list(sign = -1))
}

#' @rdname shape_constructors
#' @export
shape_increasing_convex <- function(link = linkfunctions7::log_link()) {
  ShapeCurvature(shape_name = "increasing convex", link = .check_shape_link(link),
                 shape_params = list(sign = 1, mono = 1))
}

#' @rdname shape_constructors
#' @export
shape_increasing_concave <- function(link = linkfunctions7::log_link()) {
  ShapeCurvature(shape_name = "increasing concave", link = .check_shape_link(link),
                 shape_params = list(sign = -1, mono = 1))
}

#' @rdname shape_constructors
#' @export
shape_decreasing_convex <- function(link = linkfunctions7::log_link()) {
  ShapeCurvature(shape_name = "decreasing convex", link = .check_shape_link(link),
                 shape_params = list(sign = 1, mono = -1))
}

#' @rdname shape_constructors
#' @export
shape_decreasing_concave <- function(link = linkfunctions7::log_link()) {
  ShapeCurvature(shape_name = "decreasing concave", link = .check_shape_link(link),
                 shape_params = list(sign = -1, mono = -1))
}

S7::method(print, shape_constraint) <- function(x, ...) {
  p <- x@shape_params
  extra <- if (!is.null(p$min_slope) && p$min_slope > 0)
    sprintf(", min_slope = %s", format(p$min_slope)) else ""
  if (!is.null(p$index)) extra <- sprintf("%s, located at index %d", extra, p$index)
  else if (!is.null(p$at)) extra <- sprintf("%s, located at %s", extra, format(p$at))
  cat(sprintf("<shape: %s, link = %s%s>\n", x@shape_name, x@link@link_name, extra))
  invisible(x)
}


#' The Map of a Shape Constraint
#'
#' @description
#' Returns the linear part of the map \eqn{\beta = A\,b(\gamma)} of a shape on
#' a B-spline basis with Greville abscissae `xi`: the matrix \eqn{A}, the
#' positions of the constrained coordinates and their lower bounds.
#'
#' @details
#' The first coordinate is always the level, \eqn{\beta_1 = b_1}, and the
#' first column of \eqn{A} is a column of ones, because a B-spline basis sums
#' to one. [scs()] removes that column against the intercept of the model.
#'
#' @param shape A [shape_constraint].
#' @param xi The Greville abscissae of the basis, a strictly increasing
#'   numeric vector of length \eqn{K \ge 3}.
#' @param ... Passed to methods.
#'
#' @return A list with `A`, a \eqn{K \times K} matrix; `cons`, the integer
#'   positions of the constrained coordinates; `lower`, a numeric vector of
#'   length \eqn{K}, the bound \eqn{\ell_j} added to \eqn{g^{-1}(\gamma_j)} on a
#'   constrained coordinate (zero elsewhere); and `pen`, the positions whose
#'   first differences [scs()] penalizes (the constrained slopes of a monotone
#'   shape; the second divided differences of a curvature shape, preceded by
#'   the anchored slope in a combination). A method that returns no `pen` has
#'   its `cons` penalized.
#'
#' @seealso [shape_constraint], [shape_start()].
#'
#' @examples
#' xi <- c(0, 0.2, 0.5, 0.8, 1)
#' m <- shape_map(shape_increasing(), xi)
#' m$cons
#'
#' # Every positive b gives increasing coefficients.
#' b <- c(0.3, exp(rnorm(4)))
#' all(diff(drop(m$A %*% b)) > 0)
#'
#' @export
shape_map <- S7::new_generic("shape_map", "shape",
  function(shape, xi, ...) S7::S7_dispatch())

S7::method(shape_map, ShapeMonotone) <- function(shape, xi, ...) {
  xi <- .check_greville(xi)
  K <- length(xi); h <- diff(xi); sg <- shape@shape_params$sign
  L <- lower.tri(diag(K), diag = TRUE) * 1
  # beta_j = beta_{j-1} + h_{j-1} s_j, with s_j = sign * b_j
  A <- L %*% diag(c(1, sg * h), K)
  list(A = A, cons = 2:K, lower = c(0, rep(shape@shape_params$min_slope, K - 1L)),
       pen = 2:K)
}

S7::method(shape_map, ShapeCurvature) <- function(shape, xi, ...) {
  xi <- .check_greville(xi)
  K <- length(xi); h <- diff(xi); sg <- shape@shape_params$sign
  mono <- shape@shape_params$mono
  if (is.null(mono)) mono <- 0
  # slopes s_2..s_K: s_2 = b_2 (free), s_j = s_{j-1} + sign * b_j (xi_j - xi_{j-2})/2,
  # so that a constant b gives constant second divided differences of the
  # coefficients: a parabola inside, a smaller curvature in the end panels
  # where the knots repeat (0.83 of it, measured); then
  # beta_j = beta_{j-1} + h_{j-1} s_j
  # A COMBINATION ANCHORS THE SLOPE b_2, constrained with the sign of the
  # monotonicity, at the end where the slope is smallest in absolute value:
  # the left end when the two signs agree (increasing convex, decreasing
  # concave), the right end when they differ. The second divided differences
  # b_3..b_K keep their meaning; only the end from which they are summed
  # moves, so every slope keeps the sign of the anchor.
  w <- c(0, 0, (xi[3:K] - xi[1:(K - 2L)]) / 2)
  right <- mono != 0 && mono != sg
  Ms <- matrix(0, K, K)
  for (j in 2:K) {
    Ms[j, 2L] <- if (mono == 0) 1 else mono
    if (right) {
      if (j < K) for (i in (j + 1L):K) Ms[j, i] <- -sg * w[i]
    } else if (j >= 3L) {
      for (i in 3:j) Ms[j, i] <- sg * w[i]
    }
  }
  A <- matrix(0, K, K); A[1L, 1L] <- 1
  for (j in 2:K) A[j, ] <- A[j - 1L, ] + h[j - 1L] * Ms[j, ]
  # THE ANCHOR IS THE FIRST STEP OF THE PENALIZED WALK. Left out, it is the
  # one coordinate nothing holds: where the slope at its end is near zero it
  # runs to -Inf on a flat likelihood (measured on five curves, 40 fits each:
  # a flat direction in 114 fits of 200, 39 not converged, 5 errors; inside
  # the walk none of the three, at 5-20 per cent of RMSE where the end slope
  # is exactly zero; misure_scs/fase2/anchor2.R). The free slope of a plain
  # curvature shape stays outside, as the linear part of s() does.
  list(A = A, cons = if (mono == 0) 3:K else 2:K, lower = rep(0, K),
       pen = if (mono == 0) 3:K else 2:K)
}

.shape_index <- function(shape, K) {
  M <- shape@shape_params$index
  if (is.null(M)) {
    stop(paste0("the index of the location is not set: scs() sets it from",
                " 'mode' or\n  'inflection', or a fitting layer from its pilot."),
         call. = FALSE)
  }
  if (!M %in% .shape_index_range(shape, K)) {
    stop(sprintf("the index of the location (%d) is outside %d to %d.", M,
                 min(.shape_index_range(shape, K)), K), call. = FALSE)
  }
  M
}

S7::method(shape_map, ShapeUnimodal) <- function(shape, xi, ...) {
  xi <- .check_greville(xi)
  K <- length(xi); h <- diff(xi)
  M <- .shape_index(shape, K)
  L <- lower.tri(diag(K), diag = TRUE) * 1
  # slopes s_2..s_M of one sign, s_{M+1}..s_K of the other (a valley turns
  # both): beta_j = beta_{j-1} + h s_j
  sg <- shape@shape_params$sign * c(1, ifelse(2:K <= M, 1, -1))
  sg[1L] <- 1
  list(A = L %*% diag(c(1, h) * sg, K), cons = 2:K, lower = rep(0, K), pen = 2:K)
}

S7::method(shape_map, ShapeSigmoid) <- function(shape, xi, ...) {
  xi <- .check_greville(xi)
  K <- length(xi); h <- diff(xi)
  M <- .shape_index(shape, K)
  sg <- shape@shape_params$sign
  w <- c(0, 0, (xi[3:K] - xi[1:(K - 2L)]) / 2)
  # TWO ANCHORS AND THE INCREMENTS SUMMED TOWARD THE INFLECTION. From the
  # left, s_2 = sign b_2 and s_j = s_{j-1} + sign b_j w_j up to j = M; from
  # the right, s_K = sign b_K and s_j = s_{j+1} + sign b_j w_{j+1} down to
  # j = M + 1. Every slope has the sign of the shape and grows in absolute
  # value toward the inflection, between xi_M and xi_{M+1}; M = K has no
  # right part and is a monotone convex (or concave) curve.
  Ms <- matrix(0, K, K)
  for (j in 2:K) {
    if (j <= M) {
      Ms[j, 2L] <- sg
      if (j >= 3L) for (i in 3:j) Ms[j, i] <- sg * w[i]
    } else {
      Ms[j, K] <- sg
      if (j < K) for (i in j:(K - 1L)) Ms[j, i] <- Ms[j, i] + sg * w[i + 1L]
    }
  }
  A <- matrix(0, K, K); A[1L, 1L] <- 1
  for (j in 2:K) A[j, ] <- A[j - 1L, ] + h[j - 1L] * Ms[j, ]
  list(A = A, cons = 2:K, lower = rep(0, K), pen = 2:K)
}

.check_greville <- function(xi) {
  if (!is.numeric(xi) || length(xi) < 3L || anyNA(xi) || any(diff(xi) <= 0)) {
    stop("'xi' must be a strictly increasing numeric vector of length 3 or more.",
         call. = FALSE)
  }
  as.numeric(xi)
}


#' The Starting Value of a Shape Constraint
#'
#' @description
#' Returns the starting value of \eqn{b} for a shape: constant increments,
#' read from a straight line (monotone shapes) or a parabola (curvature
#' shapes) fitted by least squares to a working response.
#'
#' @details
#' A monotone shape starts from equal slopes, the slope of the least-squares
#' line when its sign agrees with the shape. A curvature shape starts from
#' the least-squares parabola: its slope at the left end (for a combination,
#' at the anchored end) and its second derivative. A slope or a curvature of the wrong sign, or below the lower
#' bound, is replaced by a floor of one hundredth of the scale of the
#' response over the range of the covariate, so that \eqn{g(b_j - \ell_j)} is
#' finite.
#'
#' @param shape A [shape_constraint].
#' @param x The covariate.
#' @param target A working response on the scale of the linear predictor, of
#'   the length of `x`, or `NULL`, in which case the start is a flat curve at
#'   the floor.
#' @param xi The Greville abscissae of the basis.
#' @param ... Passed to methods.
#'
#' @return A numeric vector of length `length(xi)`, the starting \eqn{b} on
#'   the natural scale (the level first).
#'
#' @seealso [shape_map()], [scs()].
#'
#' @examples
#' set.seed(1)
#' x <- runif(50); y <- 2 * x + rnorm(50, sd = 0.1)
#' shape_start(shape_increasing(), x, y, xi = seq(0, 1, length.out = 6))
#'
#' @export
shape_start <- S7::new_generic("shape_start", "shape",
  function(shape, x, target, xi, ...) S7::S7_dispatch())

.shape_floor <- function(x, target) {
  sc <- if (is.null(target) || !is.finite(stats::sd(target)) || stats::sd(target) == 0) 1
        else stats::sd(target)
  0.01 * sc / max(diff(range(x)), .Machine$double.eps)
}

S7::method(shape_start, ShapeMonotone) <- function(shape, x, target, xi, ...) {
  xi <- .check_greville(xi); K <- length(xi)
  sg <- shape@shape_params$sign; lo <- shape@shape_params$min_slope
  fl <- .shape_floor(x, target)
  if (is.null(target)) return(c(0, rep(lo + fl, K - 1L)))
  cf <- stats::coef(stats::lm.fit(cbind(1, x), target))
  sl <- max(sg * cf[[2L]], lo + fl)
  c(cf[[1L]] + cf[[2L]] * xi[1L], rep(sl, K - 1L))
}

S7::method(shape_start, ShapeCurvature) <- function(shape, x, target, xi, ...) {
  xi <- .check_greville(xi); K <- length(xi)
  sg <- shape@shape_params$sign
  mono <- shape@shape_params$mono
  if (is.null(mono)) mono <- 0
  fl <- .shape_floor(x, target)
  if (is.null(target)) return(c(0, if (mono == 0) 0 else fl, rep(fl, K - 2L)))
  cf <- stats::coef(stats::lm.fit(cbind(1, x, x^2), target))
  c2 <- max(sg * 2 * cf[[3L]], fl)
  # the slope of the parabola at the anchor: the left end, or the right end
  # for a combination whose two signs differ
  xa <- if (mono != 0 && mono != sg) xi[K] else xi[1L]
  sl <- cf[[2L]] + 2 * cf[[3L]] * xa
  c(cf[[1L]] + cf[[2L]] * xi[1L] + cf[[3L]] * xi[1L]^2,
    if (mono == 0) sl else max(mono * sl, fl),
    rep(c2, K - 2L))
}


S7::method(shape_start, ShapeUnimodal) <- function(shape, x, target, xi, ...) {
  xi <- .check_greville(xi); K <- length(xi)
  M <- shape@shape_params$index
  if (is.null(M)) M <- (K + 1L) %/% 2L
  sg <- shape@shape_params$sign
  fl <- .shape_floor(x, target)
  if (is.null(target)) return(c(0, rep(fl, K - 1L)))
  # one least-squares line on each side of the extremum, a floor where a side
  # has too few points or the wrong sign
  xm <- xi[M]
  side <- function(sel) {
    if (sum(sel) < 3L) return(0)
    stats::coef(stats::lm.fit(cbind(1, x[sel]), target[sel]))[[2L]]
  }
  up <- max(sg * side(x <= xm), fl)
  dn <- max(-sg * side(x >= xm), fl)
  c(0, rep(up, max(M - 1L, 0L)), rep(dn, K - M))
}

S7::method(shape_start, ShapeSigmoid) <- function(shape, x, target, xi, ...) {
  xi <- .check_greville(xi); K <- length(xi)
  sg <- shape@shape_params$sign
  fl <- .shape_floor(x, target)
  if (is.null(target)) return(c(0, rep(fl, K - 1L)))
  # half the least-squares slope at each anchor and the floor elsewhere: a
  # straight line split between the two ends
  cf <- stats::coef(stats::lm.fit(cbind(1, x), target))
  a <- max(sg * cf[[2L]] / 2, fl)
  out <- c(0, rep(fl, K - 1L))
  out[2L] <- a
  out[K] <- a
  out
}

# b(gamma) and its derivatives of order k, componentwise: on the constrained
# coordinates lower + g^{-1}(gamma) and (g^{-1})^{(k)}(gamma); on the free
# ones gamma, 1 and 0. ONE CALL PER ORDER, each order its own generic of
# linkfunctions7, so a caller asks only for the order it uses.
.shape_b <- function(link, gamma, cons, lower, order = 0L) {
  K <- length(gamma)
  out <- switch(as.character(order),
    "0" = gamma,
    "1" = rep(1, K),
    rep(0, K))
  if (length(cons)) {
    g <- gamma[cons]
    out[cons] <- switch(as.character(order),
      "0" = lower[cons] + linkfunctions7::linkinv(link, g),
      "1" = linkfunctions7::dlinkinv(link, g),
      "2" = linkfunctions7::d2linkinv(link, g),
      "3" = linkfunctions7::d3linkinv(link, g),
      "4" = linkfunctions7::d4linkinv(link, g),
      "5" = linkfunctions7::d5linkinv(link, g),
      stop("order must be 0 to 5.", call. = FALSE))
  }
  out
}

# gamma from b on the constrained coordinates, the inverse of .shape_b(order = 0)
.shape_gamma <- function(link, b, cons, lower) {
  g <- b
  if (length(cons)) g[cons] <- linkfunctions7::linkfun(link, b[cons] - lower[cons])
  g
}
