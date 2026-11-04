#' @include scs.R
NULL

#' @title S7 Class for Shape-Constrained Tensor Product Smooths
#' @name ScteTerm
#'
#' @description
#' `ScteTerm` is the subclass of [additive_term()] for a smooth of two
#' covariates with a shape (monotone, convex or concave, or a combination) in
#' one or both of them: [scte()] constructs it. The tensor-product coefficients are written through the maps of the
#' margins' shapes, so the contribution is nonlinear in the free coefficients
#' \eqn{\gamma} and the block is its Jacobian at the current coefficients, as
#' for [scs()].
#'
#' @inheritParams additive_term
#' @inheritParams model_term
#' @param vars A list of the two covariate expressions, unevaluated.
#' @param by The `by` expression, unevaluated, or `NULL`.
#' @param shapes A list of two elements, each a monotone or curvature
#'   [shape_constraint] or `NULL` for a free margin.
#' @param spec A named list of the settings that the build uses: the
#'   `smoothers`, `default_label` and `by_hyper`.
#' @param sparse `FALSE`.
#'
#' @return An S7 object of class `ScteTerm`, inheriting from
#'   [additive_term()] and [model_term()].
#'
#' @seealso [scte()], the constructor.
#'
#' @examples
#' S7::S7_inherits(scte(x, z), ScteTerm)
#'
#' @export
ScteTerm <- S7::new_class(
  name = "ScteTerm",
  parent = additive_term,
  properties = list(
    vars = S7::class_list,
    by = S7::class_any,
    shapes = S7::class_list,
    spec = S7::class_list,
    sparse = S7::class_any
  )
)

#' Shape-Constrained Tensor Product Smooth of Two Covariates
#'
#' @description
#' Constructs a smooth function of two covariates that has a shape in one of
#' them for every value of the other, or a shape in each: monotone, convex or
#' concave, or monotone and convex or concave together. The term enters the
#' formula of any distribution parameter, like [te()].
#'
#' @details
#' # The construction
#'
#' Each covariate is expanded in the B-spline basis of its smoother, and the
#' surface is \eqn{f(x, z) = B_x(x)^\top \Theta B_z(z)} with \eqn{\Theta} the
#' \eqn{K_x \times K_z} matrix of the tensor-product coefficients. A margin
#' with a shape writes its direction of \eqn{\Theta} through the map
#' \eqn{A} of that shape on the Greville abscissae of its basis, as
#' [scs()] does for one covariate: the coordinates that the shape constrains
#' (the increments of a monotone shape, the second differences of a convex or
#' concave one) are nonnegative through the inverse of the log link, and the
#' others (the level, and the slope of a plain curvature) are free.
#'
#' With the shape on \eqn{x} alone, \eqn{\Theta = A_x B}. Every column of
#' \eqn{\Theta} has the shape, and since the basis in \eqn{z} is nonnegative
#' and sums to one, \eqn{f(\cdot, z)} is a convex combination of the columns
#' and has the shape for every \eqn{z}. The rows of \eqn{B} on the free
#' coordinates of the shape are free functions of \eqn{z}; the coefficients of
#' the level row sum to zero, the constant being left to the intercept of the
#' model.
#'
#' With shapes on both margins, \eqn{\Theta = A_x B A_z^\top}. Every column
#' of \eqn{\Theta} has the shape in \eqn{x} when an element \eqn{b_{jk}} on a
#' coordinate \eqn{j} constrained in \eqn{x} has the sign of column \eqn{k} of
#' \eqn{A_z}, and every row has the shape in \eqn{z} when an element on a
#' coordinate constrained in \eqn{z} has the sign of its column of
#' \eqn{A_x}; the free coordinates of a plain curvature are written as the two
#' hat functions of the ends of its interval, which are nonnegative. The
#' shapes hold for every coefficient vector. An element constrained in both
#' directions must then have both signs; the covariates are reflected (which
#' reverses a monotonicity and keeps a curvature) and the surface changes
#' sign where this makes the two signs agree, which it does for every pair
#' except a convex margin with a concave one. For that pair the elements
#' constrained in both directions are zero, and the curvature in one
#' covariate varies with the other only linearly. The elements free in both
#' directions sum to zero, the constant being left to the intercept.
#'
#' # The interaction under two shapes
#'
#' With shapes on both margins the cross elements \eqn{b_{jk}} are
#' constrained as well, so the interaction has a fixed sign. Under two
#' monotone shapes the mixed differences of \eqn{\Theta} are nonnegative when
#' both margins are increasing or both decreasing, and nonpositive when one is
#' increasing and the other decreasing. A surface that is increasing in both
#' covariates with a negative interaction, such as \eqn{x + z - xz/2}, is not
#' of this form, and the fit is biased for it; with the shape on one margin
#' the interaction is free. The construction of two monotone margins is that
#' of Pya and Wood (2015).
#'
#' # The penalty
#'
#' Two smoothing parameters, `lambda1` and `lambda2`, are selected with the
#' rest of the model. `lambda1` penalizes the first differences of the
#' constrained coordinates of \eqn{\gamma} along \eqn{x}, and `lambda2` those
#' along \eqn{z}; the free row or column has the second differences of its
#' coefficients along its own direction under the smoothing parameter of that
#' direction.
#'
#' # Varying the surface by another variable
#'
#' A factor `by` gives one surface per level, each with the shapes of
#' `shapes`, and the level means come from the factor written in the formula;
#' the coefficients are named `<label>.<level>.g1`, and so on. Under
#' `by_hyper = "level"` each level has its own pair of smoothing parameters,
#' named `lambda1.<level>` and `lambda2.<level>`. A numeric `by` \eqn{w}
#' gives the varying-coefficient term \eqn{w\,f(x, z)}: the level of
#' \eqn{f} is then a free and unpenalized coefficient, `<label>.const`, and
#' the default label is `scte(x, z):w`.
#'
#' # Prediction outside the data
#'
#' Beyond the interval of a basis the surface is held at its value at the
#' boundary, which preserves the shapes.
#'
#' @param x,z The covariates, bare names or expressions in the data.
#' @param shapes A list of two elements, for `x` and for `z`: each is a
#'   monotone shape ([shape_increasing()], [shape_decreasing()]), a curvature
#'   shape ([shape_convex()], [shape_concave()], [shape_increasing_convex()]
#'   and the other three combinations) or `NULL` for a free margin, and at
#'   least one is a shape. A shape here takes neither `min_slope` nor a value
#'   or a slope at a point.
#' @param smooths A \pkg{basis7} smoother whose basis is a B-spline, or a list
#'   of two, one per covariate: [basis7::bspline_smooth()] (the default, with
#'   six basis functions) or [basis7::pspline_smooth()]. Only the bases are
#'   used.
#' @param by An optional factor or numeric variable, given as a bare name or
#'   an expression in the data, as in [s()].
#' @param by_hyper With a factor `by`, whether the levels share the two
#'   smoothing parameters (`"shared"`, the default) or each level has its own
#'   pair (`"level"`).
#' @param hyper `NULL` to estimate the smoothing parameters, or values at
#'   which to hold them, named `lambda1` and `lambda2` (with the level
#'   appended under `by_hyper = "level"`, as `lambda1.a`).
#' @param id `NULL`, or strings shared with other terms whose smoothing
#'   parameters are the same, named as `hyper`.
#' @param label The term's label; `NULL` builds `scte(<x>, <z>)` from the
#'   call.
#'
#' @return An unbuilt [ScteTerm].
#'
#' @references
#' Pya, N. and Wood, S. N. (2015). Shape constrained additive models.
#' *Statistics and Computing*, 25, 543--559.
#'
#' @seealso [scs()] for one covariate, [te()] for the unconstrained tensor
#'   product.
#'
#' @examples
#' set.seed(1)
#' dd <- data.frame(x = runif(150), z = runif(150))
#' tm <- term_build(scte(x, z, shapes = list(shape_increasing(), NULL)), dd)
#' f <- term_value(tm, coef = rnorm(length(term_coef_names(tm))),
#'                 newdata = expand.grid(x = seq(0, 1, 0.1), z = 0.3))
#' all(diff(f) >= -1e-12)
#'
#' # Increasing in both covariates.
#' scte(x, z, shapes = list(shape_increasing(), shape_increasing()))
#'
#' # Convex in x, increasing and concave in z.
#' scte(x, z, shapes = list(shape_convex(), shape_increasing_concave()))
#'
#' @export
scte <- function(x, z, shapes = list(shape_increasing(), NULL),
                 smooths = basis7::bspline_smooth(k = 6), by = NULL,
                 by_hyper = c("shared", "level"), hyper = NULL, id = NULL,
                 label = NULL) {
  xe <- substitute(x)
  ze <- substitute(z)
  be <- substitute(by)
  by_hyper <- match.arg(by_hyper)
  if (!is.list(shapes) || length(shapes) != 2L) {
    stop("'shapes' must be a list of two elements, one per covariate.", call. = FALSE)
  }
  for (k in 1:2) {
    sh <- shapes[[k]]
    if (is.null(sh)) next
    if (!S7::S7_inherits(sh, ShapeMonotone) && !S7::S7_inherits(sh, ShapeCurvature)) {
      stop(paste0("each element of 'shapes' must be a monotone or a curvature shape ",
                  "(shape_increasing(), shape_convex(), shape_increasing_concave(), ...) ",
                  "or NULL."), call. = FALSE)
    }
    p <- sh@shape_params
    if (length(.shape_level_names(sh)) || any(p$min_slope > 0) || !is.null(p$value) ||
        !is.null(p[["slope"]])) {
      stop("a shape of scte() takes neither 'min_slope' nor a value or a slope at a point.",
           call. = FALSE)
    }
  }
  if (is.null(shapes[[1L]]) && is.null(shapes[[2L]])) {
    stop("at least one element of 'shapes' must be a shape; te() is the free surface.",
         call. = FALSE)
  }
  if (S7::S7_inherits(smooths, basis7::smoother)) smooths <- list(smooths, smooths)
  if (length(smooths) != 2L || !all(vapply(smooths, function(sm)
    S7::S7_inherits(sm, basis7::BsplineSmoother) || S7::S7_inherits(sm, basis7::PsplineSmoother),
    logical(1)))) {
    stop("'smooths' must be bspline_smooth() or pspline_smooth(), or a list of two.",
         call. = FALSE)
  }
  default <- is.null(label)
  if (default) {
    label <- sprintf("scte(%s, %s)", paste(deparse(xe), collapse = ""),
                     paste(deparse(ze), collapse = ""))
  }
  if (!is.character(label) || length(label) != 1L || is.na(label) || !nzchar(label)) {
    stop("'label' must be a single non-empty character string.", call. = FALSE)
  }
  per_level <- identical(by_hyper, "level")
  if (per_level && is.null(be)) {
    stop(sprintf(paste0("'%s' sets by_hyper = \"level\" and has no 'by'. One",
                        " pair of smoothing parameters per level\n  needs a factor",
                        " 'by'."), label), call. = FALSE)
  }
  nms <- c("lambda1", "lambda2")
  ScteTerm(label = label, vars = list(xe, ze), by = be, shapes = shapes, sparse = FALSE,
           spec = list(smoothers = smooths, default_label = default, by_hyper = by_hyper),
           hyper = if (per_level) as_hyper(hyper, label) else smooth_hyper(hyper, nms, label),
           ids = check_ids(id, if (per_level) NULL else nms, label),
           X = NULL, coef_names = character(0), build_info = list(), penalty = NULL)
}

# the shape that a margin carries once the surface is multiplied by `s` and
# its covariate reflected under `rev`: a reflection reverses a monotonicity
# and keeps a curvature, a change of sign reverses both
.scte_effective <- function(shape, s, rev) {
  p <- shape@shape_params
  r <- if (rev) -1 else 1
  if (S7::S7_inherits(shape, ShapeMonotone)) {
    shape@shape_params$sign <- p$sign * s * r
  } else {
    shape@shape_params$sign <- p$sign * s
    if (!is.null(p$mono)) shape@shape_params$mono <- p$mono * s * r
  }
  shape
}

# one margin: its basis and the map of its shape on the Greville abscissae,
# for the surface multiplied by `s` and the covariate reflected under `rev`
# (the basis then read in reverse); a free margin has the identity. `cons`
# are the coordinates that the shape constrains, `sigma` the sign of each
# column of the map
.scte_margin <- function(sm, v, shape, s = 1, rev = FALSE) {
  b <- basis7::smoother_basis(sm, v)
  xi <- .scs_greville(b)
  K <- length(xi)
  if (is.null(shape)) {
    return(list(basis = b, K = K, rev = FALSE, xi = xi, A = diag(K), cons = integer(0),
                free = seq_len(K), sigma = rep(1, K), shaped = FALSE, shape = NULL))
  }
  xr <- if (rev) -rev(xi) else xi
  eff <- .scte_effective(shape, s, rev)
  mp <- shape_map(eff, xr)
  sigma <- apply(mp$A, 2L, function(a) if (all(a >= 0)) 1 else if (all(a <= 0)) -1 else NA_real_)
  if (anyNA(sigma)) stop("internal: a column of the shape's map has no sign.", call. = FALSE)
  list(basis = b, K = K, rev = rev, xi = xr, A = mp$A, cons = mp$cons,
       free = setdiff(seq_len(K), mp$cons), sigma = sigma, shaped = TRUE, shape = eff)
}

# THE ORIENTATION: the sign of the surface and the reflections of the shaped
# covariates under which the fewest cross elements of B are dropped (an
# element constrained in both directions must have the sign of both columns),
# then the fewest constrained columns of negative sign, then no change of
# sign, then the fewest reflections. Monotone shapes keep the construction
# that reflects a decreasing margin.
.scte_orient <- function(term, xs) {
  sh <- term@shapes
  sm <- term@spec$smoothers
  revs <- function(k) if (is.null(sh[[k]])) FALSE else c(FALSE, TRUE)
  best <- NULL
  for (s in c(1, -1)) for (rx in revs(1L)) for (rz in revs(2L)) {
    mx <- .scte_margin(sm[[1L]], xs[[1L]], sh[[1L]], s, rx)
    mz <- .scte_margin(sm[[2L]], xs[[2L]], sh[[2L]], s, rz)
    drop <- sum(outer(mx$sigma[mx$cons], mz$sigma[mz$cons], "!="))
    neg <- sum(mx$sigma[mx$cons] < 0) + sum(mz$sigma[mz$cons] < 0)
    key <- c(drop, neg, s < 0, rx + rz)
    better <- is.null(best) || {
      d <- key - best$key
      k <- which(d != 0)
      length(k) > 0L && d[k[1L]] < 0
    }
    if (better) best <- list(key = key, mx = mx, mz = mz, s = s)
  }
  best
}

# the free coordinates of a curvature margin under a second shaped margin,
# as the two hat functions of its ends: both nonnegative, so that the
# constraint that the other margin puts on them holds at both ends and in
# between
.scte_hats <- function(mg) {
  if (!mg$shaped || length(mg$free) != 2L) return(mg)
  t <- (mg$xi - mg$xi[1L]) / (mg$xi[mg$K] - mg$xi[1L])
  mg$A[, mg$free] <- cbind(1 - t, t)
  mg$sigma[mg$free] <- 1
  mg
}

.scte_eval <- function(mg, v) {
  E <- basis7::basis_eval(mg$basis, pmin(pmax(v, mg$basis@lower), mg$basis@upper))
  if (mg$rev) E[, rev(seq_len(mg$K)), drop = FALSE] else E
}

.scte_vals <- function(term, data) {
  out <- lapply(term@vars, function(e) as.numeric(eval(e, data, baseenv())))
  if (any(lengths(out) != nrow(data))) {
    stop(sprintf("the covariates of '%s' must give one value per row.", term@label),
         call. = FALSE)
  }
  if (any(vapply(out, function(v) anyNA(v) || !all(is.finite(v)), logical(1)))) {
    stop(sprintf("the covariates of '%s' have missing or infinite values.", term@label),
         call. = FALSE)
  }
  out
}

# the row tensor product, the index of x fastest, as vec(Theta)
.scte_W <- function(bp, xs) {
  Ex <- .scte_eval(bp$mx, xs[[1L]]); Ez <- .scte_eval(bp$mz, xs[[2L]])
  Kx <- ncol(Ex); Kz <- ncol(Ez)
  (Ex[, rep(seq_len(Kx), Kz), drop = FALSE] *
     Ez[, rep(seq_len(Kz), each = Kx), drop = FALSE]) %*% bp$M
}

S7::method(term_build, ScteTerm) <- function(term, data, ...) {
  xs <- .scte_vals(term, data)
  orient <- .scte_orient(term, xs)
  mx <- orient$mx; mz <- orient$mz
  both <- mx$shaped && mz$shaped
  if (both) { mx <- .scte_hats(mx); mz <- .scte_hats(mz) }
  Kx <- mx$K; Kz <- mz$K
  idx <- function(j, k) (k - 1L) * Kx + j
  pos <- matrix(0L, Kx, Kz)
  # THE COORDINATES. vec(Theta) = (A_z (x) A_x) vec(B). Every column of
  # Theta is in the cone of the shape in x when, for each coordinate j that
  # the shape constrains, row j of B A_z' is nonnegative, which holds where
  # B_jk has the sign of column k of A_z; every row is in the cone of the
  # shape in z where B_jk, k constrained, has the sign of column j of A_x.
  # An element constrained in both directions whose two signs differ is
  # dropped. The elements free in both directions carry the constant, which
  # is left to the intercept: under two shapes they sum to zero (the free
  # columns of a margin sum to one), under one shape the level line does.
  sgn <- matrix(0, Kx, Kz)
  for (j in seq_len(Kx)) for (k in seq_len(Kz)) {
    cx <- j %in% mx$cons
    cz <- k %in% mz$cons
    sgn[j, k] <- if (cx && cz) {
      if (mx$sigma[j] == mz$sigma[k]) mx$sigma[j] else NA_real_
    } else if (cx) mz$sigma[k] else if (cz) mx$sigma[j] else 0
  }
  L <- kronecker(mz$A, mx$A)
  if (both) {
    lines <- list(which(sgn == 0))
  } else if (mx$shaped) {
    lines <- lapply(mx$free, function(j) idx(j, seq_len(Kz)))
  } else {
    lines <- lapply(mz$free, function(k) idx(seq_len(Kx), k))
  }
  # the level line alone (b_11 under two monotone shapes) is the constant
  lead <- length(lines[[1L]]) > 1L
  Ds <- lapply(seq_along(lines), function(i) {
    if (i == 1L) stats::contr.sum(max(length(lines[[i]]), 2L)) else diag(length(lines[[i]]))
  })
  if (!lead) {
    lines <- lines[-1L]
    Ds <- Ds[-1L]
  }
  free_b <- unlist(lines, use.names = FALSE)
  D <- if (length(Ds)) as.matrix(Matrix::bdiag(Ds)) else NULL
  cons_b <- which(!is.na(sgn) & sgn != 0)
  M <- cbind(if (length(free_b)) L[, free_b, drop = FALSE] %*% D,
             L[, cons_b, drop = FALSE] %*% diag(sgn[cons_b], length(cons_b)))
  M <- orient$s * M
  nf <- if (is.null(D)) 0L else ncol(D)
  pos[cons_b] <- nf + seq_along(cons_b)
  p <- ncol(M)
  # first differences of the constrained gamma along each direction
  diffs <- function(lines) {
    rows <- list()
    for (v in lines) {
      v <- v[v > 0L]
      if (length(v) > 1L) for (a in seq_len(length(v) - 1L)) {
        rows[[length(rows) + 1L]] <- replace(numeric(p), v[a:(a + 1L)], c(-1, 1))
      }
    }
    if (length(rows)) crossprod(do.call(rbind, rows)) else matrix(0, p, p)
  }
  Px <- diffs(lapply(seq_len(Kz), function(k) pos[, k]))
  Pz <- diffs(lapply(seq_len(Kx), function(j) pos[j, ]))
  # a free line under one shape has the second differences of its
  # coefficients along its own direction; the few elements free in both
  # directions under two shapes are not penalized
  if (nf > 0L && !both) {
    at <- 0L
    for (i in seq_along(lines)) {
      Di <- Ds[[i]]
      ii <- at + seq_len(ncol(Di))
      P2 <- crossprod(diff(diag(nrow(Di)), differences = 2L) %*% Di)
      if (mx$shaped) Pz[ii, ii] <- Pz[ii, ii] + P2 else Px[ii, ii] <- Px[ii, ii] + P2
      at <- at + ncol(Di)
    }
  }
  # THE BY VARIABLE, as for scs(): the multipliers G of the levels, and under
  # a numeric `by` the level of f put back in front as a free, unpenalized
  # coordinate, the coefficient of w that no other term carries
  by <- .smooth_by(term@by, data)
  kind <- if (is.null(by)) "none" else by$kind
  lev <- if (identical(kind, "factor")) by$levels else NULL
  if (identical(term@spec$by_hyper, "level") && !identical(kind, "factor")) {
    stop(sprintf(paste0("'%s' sets by_hyper = \"level\" and has no factor 'by'."),
                 term@label), call. = FALSE)
  }
  G <- .scs_G(term, kind, lev, data)
  m <- ncol(G)
  if (identical(kind, "numeric") && isTRUE(term@spec$default_label)) {
    term@label <- paste0(term@label, ":", paste(deparse(term@by), collapse = ""))
  }
  keep_level <- identical(kind, "numeric")
  if (keep_level) {
    M <- cbind(1, M)
    pad <- function(P) { Q <- matrix(0, p + 1L, p + 1L); Q[-1L, -1L] <- P; Q }
    Px <- pad(Px); Pz <- pad(Pz)
    pos[pos > 0L] <- pos[pos > 0L] + 1L
    nf <- nf + 1L
    p <- p + 1L
  }
  cons <- sort(unique(pos[pos > 0L]))
  bp <- list(mx = mx, mz = mz, M = M, pos = pos, nf = nf, s = orient$s, lines = lines,
             pieces = list(list(cons = cons, lower = numeric(p), drop_b = integer(0))),
             piece_of = rep(1L, m), link = linkfunctions7::log_link(), nc = p, m = m,
             ncoef = p * m, ix = matrix(seq_len(p * m), p, m), C = NULL,
             n = length(xs[[1L]]), xs = xs, by = term@by, by_kind = kind,
             by_levels = lev, keep_level = keep_level, G = G)
  bp$W <- list(.scte_W(bp, xs))
  bp$fixed <- list(numeric(bp$n))
  g0 <- numeric(p)
  g0[cons] <- log(.shape_floor(xs[[1L]], NULL))
  bp$coef0 <- rep(g0, m)
  gn <- if (keep_level) c("const", sprintf("g%d", seq_len(p - 1L))) else sprintf("g%d", seq_len(p))
  term@coef_names <- if (identical(kind, "factor"))
    as.character(t(outer(lev, gn, function(a, g) sprintf("%s.%s.%s", term@label, a, g))))
  else sprintf("%s.%s", term@label, gn)
  bp$penalties <- NULL
  if (identical(term@spec$by_hyper, "level")) {
    parts <- lapply(seq_len(m), function(l) {
      .penalty_entry(lev[l], (l - 1L) * p + seq_len(p),
                     penalties7::additive_penalty(list(Px, Pz)))
    })
    parts <- .entry_hyper(parts, term@hyper, term@ids, term@label,
                          "this smooth's penalty")
    term@hyper <- do.call(c, c(list(list()), lapply(parts, function(en) {
      if (!length(en$fixed) || !nzchar(en$name)) return(en$fixed)
      stats::setNames(en$fixed, paste0(names(en$fixed), ".", en$name))
    })))
    bp$penalties <- parts
    term@penalty <- NULL
  } else {
    term@penalty <- if (m == 1L) penalties7::additive_penalty(list(Px, Pz))
      else penalties7::additive_penalty(list(kronecker(diag(m), Px), kronecker(diag(m), Pz)))
  }
  term@build_info <- bp
  term@X <- matrix(0, bp$n, p * m)
  term_refresh(term, bp$coef0)
}

S7::method(term_penalties, ScteTerm) <- function(term, ...) {
  ent <- term@build_info$penalties
  if (is.null(ent)) return(S7::method(term_penalties, model_term)(term, ...))
  ent
}

# THE REFRESH AND THE DERIVATIVES OF THE BLOCK are those of scs(): the
# contribution is W b(gamma) with b acting on one coordinate at a time
S7::method(term_refresh, ScteTerm) <- function(term, coef, ...) {
  S7::method(term_refresh, ScsTerm)(term, coef, ...)
}
S7::method(term_block_deriv, ScteTerm) <- function(term, coef = NULL, v, ...) {
  S7::method(term_block_deriv, ScsTerm)(term, coef = coef, v = v, ...)
}
S7::method(term_block_contract, ScteTerm) <- function(term, coef = NULL, A, ...) {
  S7::method(term_block_contract, ScsTerm)(term, coef = coef, A = A, ...)
}
S7::method(term_block_deriv2, ScteTerm) <- function(term, coef = NULL, v, u, ...) {
  S7::method(term_block_deriv2, ScsTerm)(term, coef = coef, v = v, u = u, ...)
}
S7::method(term_block_deriv3, ScteTerm) <- function(term, coef = NULL, v, u, s, ...) {
  S7::method(term_block_deriv3, ScsTerm)(term, coef = coef, v = v, u = u, s = s, ...)
}
S7::method(term_block_deriv4, ScteTerm) <- function(term, coef = NULL, v, u, s, r, ...) {
  S7::method(term_block_deriv4, ScsTerm)(term, coef = coef, v = v, u = u, s = s,
                                         r = r, ...)
}

S7::method(term_value, ScteTerm) <- function(term, coef = NULL, newdata = NULL, ...) {
  .assert_built(term)
  bp <- term@build_info
  if (is.null(coef) && is.null(newdata)) return(bp$value)
  coef <- .scs_coef(bp, coef)
  B0 <- .scs_B(bp, coef, 0L)
  if (is.null(newdata)) return(.scs_value(bp, bp$W, bp$fixed, bp$G, B0))
  W <- .scte_W(bp, .scte_vals(term, newdata))
  .scs_value(bp, list(W), list(numeric(nrow(W))),
             .scs_G(term, bp$by_kind, bp$by_levels, newdata), B0)
}

S7::method(term_predict, ScteTerm) <- function(term, newdata, ...) {
  .assert_built(term)
  bp <- term@build_info
  W <- .scte_W(bp, .scte_vals(term, newdata))
  X <- .scs_block(bp, list(W), .scs_G(term, bp$by_kind, bp$by_levels, newdata),
                  .scs_B(bp, bp$coef, 1L))
  colnames(X) <- term@coef_names
  X
}

# the start: the least-squares plane, its slopes (in the shaped directions,
# with the sign of the shape) on the edge increments, the floor elsewhere
S7::method(term_coef_start, ScteTerm) <- function(term, target = NULL, ...) {
  .assert_built(term)
  bp <- term@build_info
  if (is.null(target) || length(target) != bp$n) return(bp$coef0)
  target <- as.numeric(target)
  out <- matrix(bp$coef0, bp$nc, bp$m)
  if (identical(bp$by_kind, "numeric")) {
    cf <- tryCatch(stats::coef(stats::lm.fit(cbind(1, bp$G[, 1L]), target))[[2L]],
                   error = function(e) NA_real_)
    if (is.finite(cf)) out[1L, 1L] <- cf
    return(as.vector(out))
  }
  pos <- bp$pos
  for (l in seq_len(bp$m)) {
    rows <- bp$G[, l] != 0
    if (sum(rows) < 4L) next
    cf <- tryCatch(stats::coef(stats::lm.fit(cbind(1, bp$xs[[1L]][rows], bp$xs[[2L]][rows]),
                                             target[rows])), error = function(e) NULL)
    if (is.null(cf) || !all(is.finite(cf))) next
    fl <- .shape_floor(bp$xs[[1L]], target)
    g <- out[, l]
    if (.scte_curved(bp)) {
      out[, l] <- .scte_curved_start(bp, g, cf, fl)
      next
    }
    sx <- max(if (bp$mx$rev) -cf[[2L]] else cf[[2L]], fl)
    sz <- max(if (bp$mz$rev) -cf[[3L]] else cf[[3L]], fl)
    if (bp$mx$shaped && bp$mz$shaped) {
      g[pos[-1L, 1L]] <- log(sx)
      g[pos[1L, -1L]] <- log(sz)
      g[pos[-1L, -1L]] <- log(fl)
    } else {
      g[pos[pos > 0L]] <- log(if (bp$mx$shaped) sx else sz)
    }
    out[, l] <- g
  }
  as.vector(out)
}

.scte_curved <- function(bp) {
  any(vapply(list(bp$mx, bp$mz), function(mg)
    mg$shaped && !S7::S7_inherits(mg$shape, ShapeMonotone), logical(1)))
}

# the start under a curvature: the plane's slope on the coordinates that
# carry a slope (every constrained coordinate of a monotone margin, the
# anchor of a monotone curvature, the free slope line of a plain curvature
# under one shape), each read with the direction in which it moves the
# surface; the floor on the curvature coordinates
.scte_curved_start <- function(bp, g, cf, fl) {
  pos <- bp$pos
  carries <- function(mg) {
    if (!mg$shaped) return(integer(0))
    if (S7::S7_inherits(mg$shape, ShapeMonotone)) return(mg$cons)
    if (!is.null(mg$shape@shape_params$mono)) return(2L)
    integer(0)
  }
  slope <- function(mg, c0) if (mg$rev) -c0 else c0
  Kx <- bp$mx$K; Kz <- bp$mz$K
  for (j in carries(bp$mx)) for (k in setdiff(seq_len(Kz), bp$mz$cons)) {
    if (pos[j, k] > 0L) {
      d <- bp$s * bp$mx$sigma[j]
      g[pos[j, k]] <- log(max(d * slope(bp$mx, cf[[2L]]), fl))
    }
  }
  for (k in carries(bp$mz)) for (j in setdiff(seq_len(Kx), bp$mx$cons)) {
    if (pos[j, k] > 0L) {
      d <- bp$s * bp$mz$sigma[k]
      g[pos[j, k]] <- log(max(d * slope(bp$mz, cf[[3L]]), fl))
    }
  }
  # the free slope line of a plain curvature under one shape
  one <- xor(bp$mx$shaped, bp$mz$shaped)
  mg <- if (bp$mx$shaped) bp$mx else bp$mz
  if (one && length(mg$free) == 2L && length(bp$lines) == 2L) {
    n1 <- length(bp$lines[[1L]]) - 1L
    ii <- n1 + seq_along(bp$lines[[2L]])
    g[ii] <- bp$s * slope(mg, cf[[if (bp$mx$shaped) 2L else 3L]])
  }
  g
}

S7::method(print, ScteTerm) <- function(x, ...) {
  desc <- vapply(1:2, function(k) {
    sh <- x@shapes[[k]]
    v <- paste(deparse(x@vars[[k]]), collapse = "")
    if (is.null(sh)) sprintf("free in %s", v) else sprintf("%s in %s", sh@shape_name, v)
  }, character(1))
  if (term_is_built(x)) {
    bp <- x@build_info
    by <- if (identical(bp$by_kind, "factor")) sprintf(", by a factor of %d levels", bp$m)
          else if (identical(bp$by_kind, "numeric")) ", by a numeric variable" else ""
    cat(sprintf("<ScteTerm> '%s' built: %d coefficients; %s%s\n", x@label, ncol(x@X),
                paste(desc, collapse = ", "), by))
  } else {
    cat(sprintf("<ScteTerm> '%s' (specification; call term_build() with data); %s\n",
                x@label, paste(desc, collapse = ", ")))
  }
  invisible(x)
}
