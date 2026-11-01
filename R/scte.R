#' @include scs.R
NULL

#' @title S7 Class for Shape-Constrained Tensor Product Smooths
#' @name ScteTerm
#'
#' @description
#' `ScteTerm` is the subclass of [additive_term()] for a smooth of two
#' covariates that is monotone in one or both of them: [scte()] constructs
#' it. The tensor-product coefficients are written through the maps of the
#' margins' shapes, so the contribution is nonlinear in the free coefficients
#' \eqn{\gamma} and the block is its Jacobian at the current coefficients, as
#' for [scs()].
#'
#' @inheritParams additive_term
#' @inheritParams model_term
#' @param vars A list of the two covariate expressions, unevaluated.
#' @param by The `by` expression, unevaluated, or `NULL`.
#' @param shapes A list of two elements, each a monotone
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
#' Constructs a smooth function of two covariates that is monotone in one of
#' them for every value of the other, or monotone in both. The term enters the
#' formula of any distribution parameter, like [te()].
#'
#' @details
#' # The construction
#'
#' Each covariate is expanded in the B-spline basis of its smoother, and the
#' surface is \eqn{f(x, z) = B_x(x)^\top \Theta B_z(z)} with \eqn{\Theta} the
#' \eqn{K_x \times K_z} matrix of the tensor-product coefficients. A margin
#' with a shape writes its direction of \eqn{\Theta} as cumulative sums of
#' increments on the Greville abscissae of its basis, as [shape_increasing()]
#' does for one covariate, each increment constrained through the inverse of
#' the log link. With the shape on \eqn{x} alone,
#' \eqn{\Theta = A_x B}: every column of \eqn{\Theta}, a function of \eqn{x}
#' at one coefficient of \eqn{z}, is monotone, and the first row of \eqn{B} is
#' a free function of \eqn{z} whose coefficients sum to zero, the constant
#' being left to the intercept of the model. With shapes on both margins,
#' \eqn{\Theta = A_x B A_z^\top} and every element of \eqn{B} except the level
#' \eqn{b_{11}} is constrained. A decreasing margin is the increasing one in
#' the reflected covariate, which reverses the order of its basis.
#'
#' # The interaction under two shapes
#'
#' With shapes on both margins every increment is constrained, the cross
#' increments \eqn{b_{jk}}, \eqn{j, k \ge 2}, included. The surface is then
#' monotone in both covariates, and its interaction has a fixed sign: the
#' mixed differences of \eqn{\Theta} are nonnegative when both margins are
#' increasing or both decreasing, and nonpositive when one is increasing and
#' the other decreasing. A surface that is increasing in both covariates with
#' a negative interaction, such as \eqn{x + z - xz/2}, is not of this form,
#' and the fit is biased for it; with the shape on one margin the interaction
#' is free. The construction is that of Pya and Wood (2015).
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
#' @param shapes A list of two elements, for `x` and for `z`: each is
#'   [shape_increasing()], [shape_decreasing()] or `NULL` for a free margin,
#'   and at least one is a shape. A shape here takes neither `min_slope` nor
#'   a value or a slope at a point.
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
    if (!S7::S7_inherits(sh, ShapeMonotone)) {
      stop("each element of 'shapes' must be shape_increasing(), shape_decreasing() or NULL.",
           call. = FALSE)
    }
    p <- sh@shape_params
    if (length(.shape_level_names(sh)) || any(p$min_slope > 0) || !is.null(p$value)) {
      stop("a shape of scte() takes neither 'min_slope' nor a value at a point.",
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

# one margin: its basis, and the Greville abscissae on which its map is
# written, reflected for a decreasing shape (the basis then read in reverse)
.scte_margin <- function(sm, v, shape) {
  b <- basis7::smoother_basis(sm, v)
  xi <- .scs_greville(b)
  dec <- !is.null(shape) && shape@shape_params$sign < 0
  xr <- if (dec) -rev(xi) else xi
  list(basis = b, K = length(xi), rev = dec, xi = xr,
       A = shape_map(shape_increasing(), xr)$A, shaped = !is.null(shape))
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
  mx <- .scte_margin(term@spec$smoothers[[1L]], xs[[1L]], term@shapes[[1L]])
  mz <- .scte_margin(term@spec$smoothers[[2L]], xs[[2L]], term@shapes[[2L]])
  Kx <- mx$K; Kz <- mz$K
  idx <- function(j, k) (k - 1L) * Kx + j
  pos <- matrix(0L, Kx, Kz)
  # THE COORDINATES: b_11 is the level, left to the intercept. With one shape
  # the first line of B across the free direction is free, its coefficients
  # summing to zero (contr.sum) so that the constant is not among them.
  if (mx$shaped && mz$shaped) {
    L <- kronecker(mz$A, mx$A)
    cons_b <- seq_len(Kx * Kz)[-1L]
    M <- L[, cons_b, drop = FALSE]
    nf <- 0L
    D <- NULL
  } else {
    L <- if (mx$shaped) kronecker(diag(Kz), mx$A) else kronecker(mz$A, diag(Kx))
    free_b <- if (mx$shaped) idx(1L, seq_len(Kz)) else idx(seq_len(Kx), 1L)
    D <- stats::contr.sum(length(free_b))
    cons_b <- setdiff(seq_len(Kx * Kz), free_b)
    M <- cbind(L[, free_b, drop = FALSE] %*% D, L[, cons_b, drop = FALSE])
    nf <- ncol(D)
  }
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
  if (nf > 0L) {
    D2 <- crossprod(diff(diag(nf + 1L), differences = 2L) %*% D)
    i <- seq_len(nf)
    if (mx$shaped) Pz[i, i] <- Pz[i, i] + D2 else Px[i, i] <- Px[i, i] + D2
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
  bp <- list(mx = mx, mz = mz, M = M, pos = pos, nf = nf,
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
    sx <- max(if (bp$mx$rev) -cf[[2L]] else cf[[2L]], fl)
    sz <- max(if (bp$mz$rev) -cf[[3L]] else cf[[3L]], fl)
    g <- out[, l]
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
