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
#' abscissa nearest to it, and the model is fitted from two starts (the
#' least-squares lines of [shape_start()] and the slopes of a smoothing
#' spline) whose better criterion is kept: away from the mode of the data
#' the first start is nearly flat and can leave the fit on the flat curve. Left at `NULL`, the
#' fitting layer chooses \eqn{M}: a pilot fit with [s()] in place of the term
#' gives the abscissa nearest to its extremum (unimodal) or to its steepest
#' point (sigmoid), \eqn{M_0}, and the model is fitted at the indices
#' \eqn{M_0 - 2, \dots, M_0 + 2} that exist; the best criterion decides (see
#' [term_pilot()] and [term_variants()]).
#'
#' # Under a factor `by`
#'
#' With a factor `by` in [scs()], `mode`, `inflection`, `index`, `valley`
#' and `decreasing` may be given by level, as a vector named by the levels
#' (`valley = c(a = FALSE, b = TRUE)` gives a maximum in level `a` and a
#' minimum in level `b`), and the point and the value of `value` may be given
#' by level in the list form. A location that is not given is chosen by
#' level (`mode_by` or `inflection_by` at `"level"`, the default): the
#' fitting layer moves the index of one level at a time, from the pilot
#' indices of the levels. With `"shared"` one index holds for every level,
#' and the versions are taken around the median of the levels' pilot
#' indices.
#'
#' @param link A \pkg{linkfunctions7} link onto \eqn{(0, \infty)}.
#' @param mode,inflection `NULL`, a single number (the value of the covariate
#'   at which the mode or the inflection is placed), or a vector of such
#'   numbers named by the levels of a factor `by`.
#' @param index `NULL`, the index \eqn{M} itself, or a vector of indices
#'   named by level. A fitting layer sets it on the versions it tries. It
#'   takes precedence over `mode` and `inflection`.
#' @param valley `FALSE` for a maximum, `TRUE` for a minimum, or a logical
#'   vector named by level.
#' @param decreasing `FALSE` for an increasing sigmoid, `TRUE` for a
#'   decreasing one, or a logical vector named by level.
#' @param value `NULL`, two finite numbers `c(at = a, value = v)`, or a list
#'   `list(at = a, value = v)` whose elements are one number or a vector
#'   named by level: the smooth is \eqn{v} at \eqn{a}, in place of the
#'   centering (see [shape_increasing()] for when \eqn{v} changes the fit).
#' @param mode_by,inflection_by Under a factor `by`, whether a location that
#'   is not given is chosen for each level (`"level"`, the default) or once
#'   for every level (`"shared"`). `"shared"` without a factor `by` is an
#'   error.
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
#'
#' # Under a factor by with levels a and b: a peak in a, a valley in b, and
#' # one location for both.
#' shape_unimodal(valley = c(a = FALSE, b = TRUE), mode_by = "shared")
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

# A SETTING BY LEVEL: one value for every level, or a vector named by the
# levels of a factor `by`, which scs() resolves at the build (.shape_level())
.check_by_level <- function(v, nm, check) {
  if (is.null(v)) return(NULL)
  nms <- names(v)
  if (length(v) == 1L && is.null(nms)) return(check(v))
  if (is.null(nms) || anyNA(nms) || any(!nzchar(nms)) || anyDuplicated(nms)) {
    stop(sprintf(paste0("'%s' must be one value, or a vector named by the levels",
                        " of the factor 'by' of scs()."), nm), call. = FALSE)
  }
  stats::setNames(unlist(lapply(seq_along(v), function(i) check(unname(v[[i]])))), nms)
}

.check_by_choice <- function(x, nm) {
  if (identical(x, c("level", "shared"))) return("level")
  if (!is.character(x) || length(x) != 1L || !x %in% c("level", "shared")) {
    stop(sprintf("'%s' must be \"level\" or \"shared\".", nm), call. = FALSE)
  }
  x
}

# a point constraint: NULL, or two finite numbers, the point and the value,
# named "at" and "value" in either order or unnamed in that order; or a list
# with `at` and `value`, each one number or a vector named by level
.check_point <- function(p, nm) {
  if (is.null(p)) return(NULL)
  if (is.list(p)) {
    if (length(p) != 2L || !setequal(names(p), c("at", "value"))) {
      stop(sprintf(paste0("'%s' as a list must have the elements 'at' and",
                          " 'value'."), nm), call. = FALSE)
    }
    num <- function(w) function(u) {
      if (!is.numeric(u) || length(u) != 1L || !is.finite(u)) {
        stop(sprintf("'%s$%s' must hold finite numbers.", nm, w), call. = FALSE)
      }
      as.numeric(u)
    }
    return(list(at = .check_by_level(p$at, sprintf("%s$at", nm), num("at")),
                value = .check_by_level(p$value, sprintf("%s$value", nm), num("value"))))
  }
  nms <- names(p)
  ok <- is.numeric(p) && length(p) == 2L && all(is.finite(p)) &&
    (is.null(nms) || setequal(nms, c("at", "value")))
  if (!ok) {
    stop(sprintf(paste0("'%s' must be NULL or two finite numbers, such as",
                        " c(at = 0, value = 0)."), nm), call. = FALSE)
  }
  p <- as.numeric(if (is.null(nms)) p else p[c("at", "value")])
  list(at = p[1L], value = p[2L])
}

#' @rdname shape_located
#' @export
shape_unimodal <- function(link = linkfunctions7::log_link(), mode = NULL,
                           index = NULL, valley = FALSE, value = NULL,
                           mode_by = c("level", "shared")) {
  valley <- .check_by_level(valley, "valley", function(u) .check_flag(u, "valley"))
  ShapeUnimodal(shape_name = if (all(valley)) "unimodal (valley)"
                             else if (any(valley)) "unimodal (peak or valley by level)"
                             else "unimodal",
                link = .check_shape_link(link),
                shape_params = list(at = .check_by_level(mode, "mode", function(u)
                                      .check_location(u, "mode")),
                                    index = .check_by_level(index, "index", .check_index),
                                    sign = ifelse(valley, -1, 1),
                                    value = .check_point(value, "value"),
                                    location_by = .check_by_choice(mode_by, "mode_by")))
}

#' @rdname shape_located
#' @export
shape_sigmoid <- function(link = linkfunctions7::log_link(), inflection = NULL,
                          index = NULL, decreasing = FALSE, value = NULL,
                          inflection_by = c("level", "shared")) {
  decreasing <- .check_by_level(decreasing, "decreasing", function(u)
    .check_flag(u, "decreasing"))
  ShapeSigmoid(shape_name = if (all(decreasing)) "sigmoid (decreasing)"
                            else if (any(decreasing)) "sigmoid (increasing or decreasing by level)"
                            else "sigmoid",
               link = .check_shape_link(link),
               shape_params = list(at = .check_by_level(inflection, "inflection", function(u)
                                     .check_location(u, "inflection")),
                                   index = .check_by_level(index, "index", .check_index),
                                   sign = ifelse(decreasing, -1, 1),
                                   value = .check_point(value, "value"),
                                   location_by = .check_by_choice(inflection_by,
                                                                  "inflection_by")))
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
#' # A value or a slope at a point
#'
#' With `value = c(at = a, value = v)` the smooth satisfies \eqn{f(a) = v}
#' exactly, and this constraint replaces the centering of [scs()]. With
#' `slope = c(at = a, value = v)`, which only `shape_convex()` and
#' `shape_concave()` accept, it satisfies \eqn{f'(a) = v} exactly: the slope
#' at the left end, the one free slope of these two shapes, is then written in
#' terms of the second divided differences. The other shapes constrain every
#' slope, so a slope at a point cannot be imposed on them by this
#' elimination. The point \eqn{a} must lie in the interval of the basis,
#' which the `lower` and `upper` arguments of the smoother can widen.
#'
#' A value other than zero changes the fit only when no other term of the
#' equation absorbs a constant (an intercept, the effects of a factor, an
#' uncentered term), as in a formula that starts with `0 +`. Otherwise the
#' fit is the same for every \eqn{v}, and \eqn{v} only moves the level
#' between the intercept and the smooth; with \eqn{v = 0} the intercept is
#' the linear predictor at \eqn{a}.
#'
#' # Under a factor `by`
#'
#' With a factor `by` in [scs()], `min_slope` may be a vector named by the
#' levels, and the point and the value of `value` or `slope` may be given by
#' level in the list form, `list(at = 0, value = c(a = 0, b = 1))`. A
#' curvature shape has one left slope per level by default; with
#' `left_slope_by = "shared"` that slope (the free one of `shape_convex()`
#' and `shape_concave()`, the anchored one of a combination) is one
#' coefficient for every level. The anchored slope of a combination is
#' penalized, so tying it is not available with `by_hyper = "level"` in
#' [scs()], which would penalize it once per level; a slope at a point
#' removes the left slope and cannot be combined with a tie.
#'
#' @param link A \pkg{linkfunctions7} link onto \eqn{(0, \infty)}, by default
#'   [linkfunctions7::log_link()]. [linkfunctions7::softplus_link()] is the
#'   other shipped choice.
#' @param min_slope A single nonnegative number, the smallest slope in
#'   absolute value, or a vector of such numbers named by the levels of a
#'   factor `by`.
#' @param value `NULL`, two finite numbers `c(at = a, value = v)`, or a list
#'   `list(at = a, value = v)` whose elements are one number or a vector
#'   named by level: the smooth is \eqn{v} at \eqn{a}. Unnamed, the first
#'   number is the point.
#' @param slope `NULL`, or the same forms as `value`: the derivative of the
#'   smooth is \eqn{v} at \eqn{a}.
#' @param left_slope_by Under a factor `by`, whether the left slope of a
#'   curvature shape is a coefficient of each level (`"level"`, the default)
#'   or one coefficient for every level (`"shared"`). `"shared"` without a
#'   factor `by` is an error.
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
#' shape_increasing(value = c(at = 0, value = 0))
#' shape_convex(slope = c(at = 0.5, value = 0))
#'
#' # Settings by level, for a factor by with levels a and b.
#' shape_increasing(min_slope = c(a = 0, b = 0.2))
#' shape_increasing(value = list(at = 0, value = c(a = 0, b = 1)))
#' shape_convex(left_slope_by = "shared")
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
shape_increasing <- function(link = linkfunctions7::log_link(), min_slope = 0,
                             value = NULL) {
  ShapeMonotone(shape_name = "increasing", link = .check_shape_link(link),
                shape_params = list(sign = 1,
                                    min_slope = .check_by_level(min_slope, "min_slope",
                                                                .check_min_slope),
                                    value = .check_point(value, "value")))
}

#' @rdname shape_constructors
#' @export
shape_decreasing <- function(link = linkfunctions7::log_link(), min_slope = 0,
                             value = NULL) {
  ShapeMonotone(shape_name = "decreasing", link = .check_shape_link(link),
                shape_params = list(sign = -1,
                                    min_slope = .check_by_level(min_slope, "min_slope",
                                                                .check_min_slope),
                                    value = .check_point(value, "value")))
}

#' @rdname shape_constructors
#' @export
shape_convex <- function(link = linkfunctions7::log_link(), value = NULL,
                         slope = NULL, left_slope_by = c("level", "shared")) {
  ShapeCurvature(shape_name = "convex", link = .check_shape_link(link),
                 shape_params = list(sign = 1, value = .check_point(value, "value"),
                                     slope = .check_point(slope, "slope"),
                                     left_slope_by = .check_by_choice(left_slope_by,
                                                                      "left_slope_by")))
}

#' @rdname shape_constructors
#' @export
shape_concave <- function(link = linkfunctions7::log_link(), value = NULL,
                          slope = NULL, left_slope_by = c("level", "shared")) {
  ShapeCurvature(shape_name = "concave", link = .check_shape_link(link),
                 shape_params = list(sign = -1, value = .check_point(value, "value"),
                                     slope = .check_point(slope, "slope"),
                                     left_slope_by = .check_by_choice(left_slope_by,
                                                                      "left_slope_by")))
}

#' @rdname shape_constructors
#' @export
shape_increasing_convex <- function(link = linkfunctions7::log_link(), value = NULL,
                     left_slope_by = c("level", "shared")) {
  ShapeCurvature(shape_name = "increasing convex", link = .check_shape_link(link),
                 shape_params = list(sign = 1, mono = 1,
                                     value = .check_point(value, "value"),
                                     left_slope_by = .check_by_choice(left_slope_by,
                                                                      "left_slope_by")))
}

#' @rdname shape_constructors
#' @export
shape_increasing_concave <- function(link = linkfunctions7::log_link(), value = NULL,
                     left_slope_by = c("level", "shared")) {
  ShapeCurvature(shape_name = "increasing concave", link = .check_shape_link(link),
                 shape_params = list(sign = -1, mono = 1,
                                     value = .check_point(value, "value"),
                                     left_slope_by = .check_by_choice(left_slope_by,
                                                                      "left_slope_by")))
}

#' @rdname shape_constructors
#' @export
shape_decreasing_convex <- function(link = linkfunctions7::log_link(), value = NULL,
                     left_slope_by = c("level", "shared")) {
  ShapeCurvature(shape_name = "decreasing convex", link = .check_shape_link(link),
                 shape_params = list(sign = 1, mono = -1,
                                     value = .check_point(value, "value"),
                                     left_slope_by = .check_by_choice(left_slope_by,
                                                                      "left_slope_by")))
}

#' @rdname shape_constructors
#' @export
shape_decreasing_concave <- function(link = linkfunctions7::log_link(), value = NULL,
                     left_slope_by = c("level", "shared")) {
  ShapeCurvature(shape_name = "decreasing concave", link = .check_shape_link(link),
                 shape_params = list(sign = -1, mono = -1,
                                     value = .check_point(value, "value"),
                                     left_slope_by = .check_by_choice(left_slope_by,
                                                                      "left_slope_by")))
}

S7::method(print, shape_constraint) <- function(x, ...) {
  p <- x@shape_params
  # one value, or the values by level in parentheses
  fv <- function(v) if (is.null(names(v))) format(v)
    else sprintf("(%s)", paste(names(v), format(v), sep = " = ", collapse = ", "))
  extra <- if (!is.null(p$min_slope) && any(p$min_slope > 0))
    sprintf(", min_slope = %s", fv(p$min_slope)) else ""
  if (!is.null(p$index)) extra <- sprintf("%s, located at index %s", extra, fv(p$index))
  else if (!is.null(p$at)) extra <- sprintf("%s, located at %s", extra, fv(p$at))
  if (identical(p$location_by, "shared")) {
    extra <- sprintf("%s, one location for every level", extra)
  }
  if (identical(p$left_slope_by, "shared")) {
    extra <- sprintf("%s, one left slope for every level", extra)
  }
  if (!is.null(p$value)) {
    extra <- sprintf("%s, f(%s) = %s", extra, fv(p$value$at), fv(p$value$value))
  }
  if (!is.null(p[["slope"]])) {
    extra <- sprintf("%s, f'(%s) = %s", extra, fv(p[["slope"]]$at), fv(p[["slope"]]$value))
  }
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
  function(shape, xi, ...) {
    .shape_require_one(shape)
    S7::S7_dispatch()
  })

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
#' A unimodal shape starts from the least-squares line on each side of its
#' location, and a sigmoid from half the least-squares slope at each anchor.
#' With `rule = "spline"` passed through `...`, both start instead from the
#' slopes of a smoothing spline ([stats::smooth.spline()]) on the segments
#' between the Greville abscissae, each given the sign that the shape
#' requires at its position and floored; for a sigmoid the slopes are made
#' to grow up to the inflection and shrink after it before they are written
#' as anchors and increments. [scs()] fits a term whose location is given
#' from both starts and keeps the better (see [term_variants()]).
#'
#' @param shape A [shape_constraint].
#' @param x The covariate.
#' @param target A working response on the scale of the linear predictor, of
#'   the length of `x`, or `NULL`, in which case the start is a flat curve at
#'   the floor.
#' @param xi The Greville abscissae of the basis.
#' @param ... Passed to methods; `rule = "spline"` for the second start of
#'   a located shape.
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
  function(shape, x, target, xi, ...) {
    .shape_require_one(shape)
    S7::S7_dispatch()
  })

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


# THE SLOPES OF A SMOOTHING SPLINE on the Greville segments, the second
# start of a located shape (rule = "spline"): where the location is far from
# the mode of the data, the lines of the first start have the wrong sign on
# a side, the start is nearly flat, and the flat curve, stationary for every
# smoothing parameter, holds the fit (misure_scs M15)
.shape_spline_slopes <- function(x, target, xi) {
  sp <- tryCatch(stats::smooth.spline(x, target), error = function(e) NULL)
  if (is.null(sp)) return(NULL)
  sl <- diff(stats::predict(sp, xi)$y) / diff(xi)
  if (all(is.finite(sl))) sl else NULL
}

S7::method(shape_start, ShapeUnimodal) <- function(shape, x, target, xi, ...) {
  xi <- .check_greville(xi); K <- length(xi)
  M <- shape@shape_params$index
  if (is.null(M)) M <- (K + 1L) %/% 2L
  sg <- shape@shape_params$sign
  fl <- .shape_floor(x, target)
  if (is.null(target)) return(c(0, rep(fl, K - 1L)))
  if (identical(list(...)$rule, "spline")) {
    sl <- .shape_spline_slopes(x, target, xi)
    if (!is.null(sl)) return(c(0, pmax(sg * ifelse(2:K <= M, 1, -1) * sl, fl)))
  }
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
  if (identical(list(...)$rule, "spline")) {
    sl <- .shape_spline_slopes(x, target, xi)
    M <- shape@shape_params$index
    if (is.null(M)) M <- (K + 2L) %/% 2L
    if (!is.null(sl)) {
      # the slopes of the shape's sign, growing up to the inflection and
      # shrinking after it, then the two anchors and the increments
      a <- c(NA, pmax(sg * sl, fl))
      if (M >= 3L) for (j in 3:M) a[j] <- max(a[j], a[j - 1L])
      if (M < K - 1L) for (j in (K - 1L):(M + 1L)) a[j] <- max(a[j], a[j + 1L])
      w <- c(0, 0, (xi[3:K] - xi[1:(K - 2L)]) / 2)
      b <- numeric(K)
      b[2L] <- a[2L]
      if (M >= 3L) for (j in 3:M) b[j] <- max((a[j] - a[j - 1L]) / w[j], fl)
      if (M < K) {
        b[K] <- a[K]
        if (M + 1L <= K - 1L) for (j in (M + 1L):(K - 1L)) b[j] <- max((a[j] - a[j + 1L]) / w[j + 1L], fl)
      }
      return(b)
    }
  }
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


# THE LEVELS THAT A SHAPE NAMES in its settings by level, and the shape of
# one level: every setting given by level is replaced by that level's value
.shape_level_names <- function(shape) {
  p <- shape@shape_params
  vals <- list(p$min_slope, p$at, p$index, p$sign, p$value$at, p$value$value,
               p[["slope"]]$at, p[["slope"]]$value)
  unique(unlist(lapply(vals, names)))
}

.shape_require_one <- function(shape) {
  if (length(.shape_level_names(shape))) {
    stop(paste0("the shape carries values by level, which scs() resolves for",
                " each level of\n  its factor 'by'."), call. = FALSE)
  }
  invisible(NULL)
}

.shape_level <- function(shape, lev) {
  pick <- function(v, nm) {
    if (is.null(v) || is.null(names(v))) return(v)
    if (!lev %in% names(v)) {
      stop(sprintf("'%s' has no value for the level '%s' of the factor 'by'.", nm, lev),
           call. = FALSE)
    }
    unname(v[[lev]])
  }
  sig <- S7::S7_inherits(shape, ShapeSigmoid)
  p <- shape@shape_params
  p$min_slope <- pick(p$min_slope, "min_slope")
  p$at <- pick(p$at, if (sig) "inflection" else "mode")
  p$index <- pick(p$index, "index")
  p$sign <- pick(p$sign, if (sig) "decreasing" else "valley")
  if (!is.null(p$value)) {
    p$value <- list(at = pick(p$value$at, "value$at"),
                    value = pick(p$value$value, "value$value"))
  }
  if (!is.null(p[["slope"]])) {
    p[["slope"]] <- list(at = pick(p[["slope"]]$at, "slope$at"),
                    value = pick(p[["slope"]]$value, "slope$value"))
  }
  shape@shape_params <- p
  if (S7::S7_inherits(shape, ShapeUnimodal)) {
    shape@shape_name <- if (p$sign < 0) "unimodal (valley)" else "unimodal"
  } else if (sig) {
    shape@shape_name <- if (p$sign < 0) "sigmoid (decreasing)" else "sigmoid"
  }
  shape
}
