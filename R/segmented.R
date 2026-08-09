#' @include term_classes.R generics.R nonlinear.R
NULL

#' @title S7 Class for Segmented and Stepmented Terms
#' @name SegTerm
#'
#' @description
#' A subclass of \code{\link{additive_term}} for a covariate whose effect
#' changes at estimated break-points: a change of slope
#' (\code{\link{seg}}), a jump in level (\code{\link{jump}}), or both at
#' the same points (\code{\link{jseg}}). The design block is the working
#' one of the iteration that estimates the break-points, and is
#' recomputed by \code{\link{term_refresh}} as they move.
#'
#' @inheritParams additive_term
#' @param kind Which of the three constructions.
#' @param var The covariate expression.
#' @param npsi The number of break-points.
#' @param by An optional grouping expression.
#' @param linear Whether the block carries the linear effect.
#' @param penalty_kind The penalty on the changes, if any.
#' @param spec The resolved construction settings.
#'
#' @return An object of class \code{SegTerm}.
#'
#' @seealso \code{\link{seg}}
#' @examples
#' S7::S7_inherits(seg(x), SegTerm)
#' @export
SegTerm <- S7::new_class(
  name = "SegTerm",
  parent = additive_term,
  properties = list(
    kind = S7::class_character,
    var = S7::class_any,
    npsi = S7::class_integer,
    by = S7::class_any,
    linear = S7::class_logical,
    penalty_kind = S7::class_character,
    spec = S7::class_list
  )
)

#' Segmented, Stepmented and Segmented-with-Jump Terms
#'
#' @description
#' A covariate whose effect changes at break-points estimated with
#' everything else. \code{seg} changes slope at each break-point and
#' stays continuous (\cite{muggeo2003}); \code{jump} steps to a new level
#' and is discontinuous (\cite{fasola2018}); \code{jseg} does both at the
#' same points.
#'
#' @details
#' Both constructions rest on the same device: a quantity that depends on
#' the break-point non-linearly is written as a linear function of it
#' once something is frozen at the previous iterate, so that a linear fit
#' returns the new break-point. Iterating that is the estimation
#' algorithm, and here it is the \code{\link{term_refresh}} contract of
#' \code{\link{nl}} with a different block: refreshing the term at the
#' current coefficients and fitting is one step of it.
#'
#' \subsection{The continuous case}{
#' With \eqn{f(x) = \delta (x-\psi)_+} the contribution is differentiable
#' in \eqn{\psi} away from the break-point, and
#' \eqn{\partial f/\partial\psi = -\delta\,\mathbb{1}(x>\psi)}, so the
#' design block is the ordinary Jacobian and the break-point is an
#' ordinary coefficient, updated by the increment a Gauss-Newton step
#' solves for. That is the algorithm of \cite{muggeo2003} written in the
#' coordinates the rest of this package uses: his working variables
#' \eqn{U} and \eqn{V} are its columns, and his update
#' \eqn{\psi \leftarrow \psi + \gamma/\delta} is that step.
#' }
#'
#' \subsection{The discontinuous case}{
#' A jump is not differentiable in \eqn{\psi} at all, and yet needs no
#' search over candidate positions. The identity
#' \deqn{\mathbb{1}(x>\psi) = \frac{1}{2}
#'   + \frac{x-\psi}{2\lvert x-\psi \rvert}}
#' holds exactly for \eqn{x \neq \psi}, and is linear in \eqn{\psi} once
#' the weight \eqn{1/(2\lvert x - \psi\rvert)} is held at the previous
#' iterate. Writing \eqn{W = 1/(2\lvert x-\psi^{0}\rvert)} and
#' \eqn{Z = xW + 1/2}, a jump of size \eqn{\kappa} at \eqn{\psi} is
#' \eqn{\kappa Z + gW} with \eqn{g = -\kappa\psi}, so a linear fit on
#' \eqn{(Z, W)} returns the break-point as \eqn{\psi = -g/\kappa}
#' (\cite{fasola2018}). The break-point is therefore not a coefficient
#' here but a quantity read off two of them, which is why refreshing the
#' term recovers it before rebuilding the weights.
#'
#' The weight is unbounded as \eqn{x} approaches \eqn{\psi}, and that
#' matters more than it looks. Since \eqn{Z - \psi W} is the step
#' itself, \eqn{Z} is \eqn{\psi W} plus a quantity of order one: let
#' \eqn{W} grow without bound and the two columns become numerically
#' collinear, drowning the very signal the fit is meant to read. The
#' denominator is therefore held at or above \code{band} times the
#' covariate's range, which caps \eqn{W}.
#'
#' That is a bandwidth and not a guard. Within the band the step is
#' replaced by a ramp, so the fixed point of the iteration is that of a
#' slightly smoothed problem; a narrower band is more faithful and worse
#' conditioned. The segmented literature makes the same trade by
#' displacing the observations nearest a break-point instead of capping
#' the weight.
#' }
#'
#' \subsection{What the term carries}{
#' With \code{linear = TRUE} the block carries the linear effect too, so
#' the term is the whole relationship rather than the change in it.
#' \code{by} gives an independent set of break-points and changes per
#' level of a factor. \code{penalty} puts a penalty on the changes
#' themselves -- the slope changes for \code{seg}, the jump sizes for
#' \code{jump} -- through a map that selects those coefficients and
#' leaves the linear effect and the break-points alone; with
#' \code{"lasso"} that is a selection of how many break-points are
#' really there.
#'
#' A break-point is confined to the interval between the 5th and the 95th
#' percentile of the covariate. Outside it the block is singular rather
#' than merely ill-conditioned: with \eqn{\psi} below the smallest
#' observation the indicator is constant, so the truncated line and that
#' constant are linearly dependent, and with the linear effect present so
#' are all three columns. A run that ends against the limit has not
#' located a break-point, and \code{\link{seg_psi}} then returns the
#' limit itself.
#'
#' The objective has local optima in the break-points, and the iteration
#' converges from within a basin around the true position rather than
#' from anywhere. On a joint jump and change of slope at \eqn{x = 5} in
#' 500 observations, swept over eight samples and damping factors from
#' 0.05 to 1, every start at 4 or above recovers the break-point at
#' every damping below 1, every start at 2 or below fails at all of
#' them, and a start at 3 succeeds for some samples and not others. The
#' step also has to be damped for its own sake: taken whole it
#' overshoots even from a good start. A run should therefore be started
#' from several positions, which is what
#' \code{\link[optimizers7]{multistart}} does and what the bootstrap
#' restarting of the segmented literature is for.
#' }
#'
#' @param x The covariate, an expression evaluated in the data.
#' @param npsi The number of break-points. Defaults to 1.
#' @param psi Optional starting positions; defaults to evenly spaced
#'   quantiles of the covariate.
#' @param by An optional factor giving an independent set of
#'   break-points per level.
#' @param linear Whether the block carries the linear effect. Defaults to
#'   \code{TRUE}.
#' @param penalty One of \code{"none"} (default), \code{"lasso"} or
#'   \code{"ridge"}, applied to the changes.
#' @param band For a discontinuous term, the half-width of the band
#'   around a break-point over which the step is replaced by a ramp, as a
#'   fraction of the covariate's range. Defaults to \code{0.02}; see
#'   Details.
#' @param label A single non-empty string prefixed to the coefficient
#'   names.
#'
#' @return An object of class \code{\link{SegTerm}} (a specification; see
#'   \code{\link{term_build}}).
#'
#' @references
#' Muggeo, V. M. R. (2003). Estimating regression models with unknown
#' break-points. \emph{Statistics in Medicine}, 22(19), 3055--3071.
#'
#' Fasola, S., Muggeo, V. M. R. and Kuchenhoff, H. (2018). A heuristic,
#' iterative algorithm for change-point detection in abrupt change
#' models. \emph{Computational Statistics}, 33, 997--1015.
#'
#' @examples
#' set.seed(1)
#' dd <- data.frame(x = sort(runif(200, 0, 10)))
#' dd$y <- 1 + 0.5 * dd$x + 2 * pmax(dd$x - 6, 0) + rnorm(200, sd = 0.3)
#'
#' built <- term_build(seg(x), dd)
#' term_coef_names(built)
#' seg_psi(built)
#'
#' @export
seg <- function(x, npsi = 1, psi = NULL, by = NULL, linear = TRUE,
                penalty = c("none", "lasso", "ridge"), band = 0.02,
                label = "seg") {
  .seg_spec("seg", substitute(x), npsi, psi, substitute(by), linear,
            match.arg(penalty), band, label)
}

#' @rdname seg
#' @export
jump <- function(x, npsi = 1, psi = NULL, by = NULL, linear = TRUE,
                 penalty = c("none", "lasso", "ridge"), band = 0.02,
                 label = "jump") {
  .seg_spec("jump", substitute(x), npsi, psi, substitute(by), linear,
            match.arg(penalty), band, label)
}

#' @rdname seg
#' @export
jseg <- function(x, npsi = 1, psi = NULL, by = NULL, linear = TRUE,
                 penalty = c("none", "lasso", "ridge"), band = 0.02,
                 label = "jseg") {
  .seg_spec("jseg", substitute(x), npsi, psi, substitute(by), linear,
            match.arg(penalty), band, label)
}

.seg_spec <- function(kind, var, npsi, psi, by, linear, penalty, band,
                      label) {
  if (!is.numeric(band) || length(band) != 1L || is.na(band) ||
      band <= 0 || band >= 0.5) {
    stop("'band' must be a single number strictly between 0 and 0.5.",
         call. = FALSE)
  }
  if (!is.numeric(npsi) || length(npsi) != 1L || is.na(npsi) || npsi < 1 ||
      npsi != round(npsi)) {
    stop("'npsi' must be a single integer of at least 1.", call. = FALSE)
  }
  if (!is.null(psi) && (!is.numeric(psi) || length(psi) != npsi ||
                        anyNA(psi))) {
    stop("'psi' must give one starting position per break-point.",
         call. = FALSE)
  }
  if (!is.logical(linear) || length(linear) != 1L || is.na(linear)) {
    stop("'linear' must be TRUE or FALSE.", call. = FALSE)
  }
  if (!is.character(label) || length(label) != 1L || is.na(label) ||
      !nzchar(label)) {
    stop("'label' must be a single non-empty character string.",
         call. = FALSE)
  }
  SegTerm(label = label, kind = kind, var = var, npsi = as.integer(npsi),
          by = by, linear = linear, penalty_kind = penalty,
          spec = list(psi = psi, band = band),
          X = NULL, coef_names = character(0),
          blueprint = list(), penalty = NULL)
}

# the coefficient names of one level, and which of them are the changes
# a penalty may reach
.seg_names <- function(kind, npsi, linear) {
  k <- seq_len(npsi)
  nm <- if (linear) "lin" else character(0)
  chg <- character(0)
  if (kind %in% c("seg", "jseg")) {
    nm <- c(nm, paste0("delta", k))
    chg <- c(chg, paste0("delta", k))
  }
  if (kind == "seg") {
    nm <- c(nm, paste0("psi", k))
  } else {
    nm <- c(nm, paste0("kappa", k), paste0("g", k))
    chg <- c(chg, paste0("kappa", k))
  }
  list(names = nm, changes = chg)
}

# the break-points implied by one level's coefficients: a coefficient in
# the continuous case, and -g/kappa in the discontinuous one.
#
# They are held inside `lim`, and that is not cosmetic. A break-point
# outside the data leaves the indicator constant, and then the linear
# effect, the truncated line and that constant are linearly dependent:
# the working block goes exactly singular. Confining the break-point to
# the interval where data lie on both sides is what the segmented
# literature does for the same reason.
.seg_psi_of <- function(kind, cf, npsi, linear, lim) {
  off <- if (linear) 1L else 0L
  psi <- if (kind == "seg") {
    cf[off + npsi + seq_len(npsi)]
  } else {
    d <- if (kind == "jseg") npsi else 0L
    kap <- cf[off + d + seq_len(npsi)]
    g <- cf[off + d + npsi + seq_len(npsi)]
    kap[abs(kap) < 1e-12] <- 1e-12
    -g / kap
  }
  pmin(pmax(psi, lim[1L]), lim[2L])
}

# The working block and the true contribution of one level. The
# arithmetic is elementwise and is compiled (src/seg_block.cpp), which
# measures 1.2 to 3.2 times the R form below over n from 1e3 to 1e6 and
# one to five break-points. The operations are the same in the same
# order, so the two agree to a rounding; .seg_block_r is kept as the
# twin the tests compare against.
.seg_block <- function(kind, xv, cf, npsi, linear, floor_w, lim) {
  psi <- .seg_psi_of(kind, cf, npsi, linear, lim)
  off <- if (linear) 1L else 0L
  del <- if (kind %in% c("seg", "jseg")) cf[off + seq_len(npsi)] else numeric(0)
  d <- if (kind == "jseg") npsi else 0L
  kap <- if (kind == "seg") numeric(0) else cf[off + d + seq_len(npsi)]
  out <- seg_block_cpp(switch(kind, seg = 0L, jump = 1L, jseg = 2L),
                       xv, psi, del, kap,
                       if (linear) cf[1L] else 0, linear, floor_w)
  list(X = out$X, value = out$value, psi = psi)
}

.seg_block_r <- function(kind, xv, cf, npsi, linear, floor_w, lim) {
  n <- length(xv)
  off <- if (linear) 1L else 0L
  psi <- .seg_psi_of(kind, cf, npsi, linear, lim)
  cols <- list()
  if (linear) cols[[1L]] <- xv
  value <- if (linear) cf[1L] * xv else numeric(n)

  if (kind %in% c("seg", "jseg")) {
    del <- cf[off + seq_len(npsi)]
    for (j in seq_len(npsi)) {
      cols[[length(cols) + 1L]] <- pmax(xv - psi[j], 0)
      value <- value + del[j] * pmax(xv - psi[j], 0)
    }
  }

  if (kind == "seg") {
    # the Jacobian in the break-point: an ordinary column, so the
    # break-point is an ordinary coefficient
    for (j in seq_len(npsi)) {
      cols[[length(cols) + 1L]] <- -del[j] * (xv > psi[j])
    }
  } else {
    d <- if (kind == "jseg") npsi else 0L
    kap <- cf[off + d + seq_len(npsi)]
    for (j in seq_len(npsi)) {
      # the identity of Fasola et al.: exact at x != psi, and linear in
      # the break-point once the weight is frozen
      den <- pmax(2 * abs(xv - psi[j]), floor_w)
      W <- 1 / den
      cols[[length(cols) + 1L]] <- xv * W + 0.5
      cols[[length(cols) + 1L]] <- W
      value <- value + kap[j] * (xv > psi[j])
    }
    # the two columns of a jump are interleaved above; reorder so that
    # every kappa precedes every g, matching .seg_names()
    base <- length(cols) - 2L * npsi
    zw <- cols[base + seq_len(2L * npsi)]
    cols[base + seq_len(2L * npsi)] <-
      c(zw[seq(1L, 2L * npsi, by = 2L)], zw[seq(2L, 2L * npsi, by = 2L)])
  }

  list(X = do.call(cbind, cols), value = value, psi = psi)
}

.seg_x <- function(expr, data) {
  v <- as.numeric(eval(expr, data, baseenv()))
  if (length(v) != nrow(data)) {
    stop("the covariate of a segmented term must give one value per row.",
         call. = FALSE)
  }
  if (anyNA(v)) {
    stop("the covariate of a segmented term must not contain missing values.",
         call. = FALSE)
  }
  v
}

.seg_by <- function(by, data, levels = NULL) {
  if (is.null(by)) return(NULL)
  v <- eval(by, data, baseenv())
  f <- if (is.null(levels)) factor(v) else factor(v, levels = levels)
  if (any(is.na(f) & !is.na(v))) {
    stop("a 'by' level absent at build time cannot be predicted.",
         call. = FALSE)
  }
  f
}

# the block over every level of by, and the contribution
.seg_assemble <- function(bp, xv, grp, coef) {
  per <- bp$per_level
  m <- length(bp$levels)
  n <- length(xv)
  X <- matrix(0, n, m * per)
  value <- numeric(n)
  psi <- numeric(0)
  for (l in seq_len(m)) {
    idx <- (l - 1L) * per + seq_len(per)
    rows <- if (is.null(grp)) rep(TRUE, n) else grp == bp$levels[l]
    b <- .seg_block(bp$kind, xv[rows], coef[idx], bp$npsi, bp$linear,
                    bp$floor_w, bp$lim)
    X[rows, idx] <- b$X
    value[rows] <- b$value
    psi <- c(psi, b$psi)
  }
  list(X = X, value = value, psi = psi)
}

S7::method(term_build, SegTerm) <- function(term, data, ...) {
  xv <- .seg_x(term@var, data)
  grp <- .seg_by(term@by, data)
  levs <- if (is.null(grp)) "" else levels(grp)
  nmi <- .seg_names(term@kind, term@npsi, term@linear)
  per <- length(nmi$names)

  rng <- diff(range(xv))
  if (rng <= 0) {
    stop("the covariate of a segmented term must vary.", call. = FALSE)
  }
  # the cap on the weight of a jump. Z is psi * W plus a quantity of
  # order one, so an unbounded W makes the two columns collinear and the
  # fit unreadable; the band is what keeps them apart, at the price of a
  # ramp instead of a step within it
  floor_w <- 2 * term@spec$band * rng

  # the interval a break-point is held in: far enough inside the data
  # that both sides carry observations
  lim <- as.numeric(stats::quantile(xv, c(0.05, 0.95), names = FALSE))

  start_psi <- if (!is.null(term@spec$psi)) term@spec$psi else {
    as.numeric(stats::quantile(xv, seq_len(term@npsi) / (term@npsi + 1),
                               names = FALSE))
  }

  # the starting coefficients: unit changes, the break-points at their
  # starting positions, and g = -kappa * psi where a jump reads them off
  coef0 <- numeric(length(levs) * per)
  for (l in seq_along(levs)) {
    idx <- (l - 1L) * per
    off <- if (term@linear) 1L else 0L
    if (term@linear) coef0[idx + 1L] <- 0
    if (term@kind %in% c("seg", "jseg")) {
      coef0[idx + off + seq_len(term@npsi)] <- 1
    }
    if (term@kind == "seg") {
      coef0[idx + off + term@npsi + seq_len(term@npsi)] <- start_psi
    } else {
      d <- if (term@kind == "jseg") term@npsi else 0L
      coef0[idx + off + d + seq_len(term@npsi)] <- 1
      coef0[idx + off + d + term@npsi + seq_len(term@npsi)] <- -start_psi
    }
  }

  bp <- list(kind = term@kind, npsi = term@npsi, linear = term@linear,
             per_level = per, levels = levs, floor_w = floor_w,
             lim = lim,
             var = term@var, by = term@by, coef = coef0,
             xv = xv, grp = grp,
             names_one = nmi$names, changes_one = nmi$changes)

  asm <- .seg_assemble(bp, xv, grp, coef0)
  cn <- if (is.null(grp)) {
    paste(term@label, nmi$names, sep = ".")
  } else {
    as.character(t(outer(levs, nmi$names,
                         function(a, b) paste(term@label, a, b, sep = "."))))
  }
  X <- asm$X
  colnames(X) <- cn

  pen <- NULL
  if (term@penalty_kind != "none") {
    # a map selecting the changes: the linear effect and the break-points
    # are not shrunk, only the sizes of the changes
    keep <- which(rep(nmi$names %in% nmi$changes, times = length(levs)))
    D <- matrix(0, length(keep), ncol(X))
    D[cbind(seq_along(keep), keep)] <- 1
    pen <- switch(term@penalty_kind,
      lasso = penalties7::lasso_penalty(map = D),
      ridge = penalties7::ridge_penalty(map = D))
  }

  bp$value <- asm$value
  bp$psi <- asm$psi
  term@X <- X
  term@coef_names <- cn
  term@blueprint <- bp
  term@penalty <- pen
  term
}

S7::method(term_refresh, SegTerm) <- function(term, coef, ...) {
  .assert_built(term)
  bp <- term@blueprint
  coef <- as.numeric(coef)
  if (length(coef) != ncol(term@X)) {
    stop(sprintf("'coef' must have length %d.", ncol(term@X)), call. = FALSE)
  }
  asm <- .seg_assemble(bp, bp$xv, bp$grp, coef)
  X <- asm$X
  colnames(X) <- term@coef_names
  bp$coef <- coef
  bp$value <- asm$value
  bp$psi <- asm$psi
  term@X <- X
  term@blueprint <- bp
  term
}

S7::method(term_value, SegTerm) <- function(term, coef = NULL, ...) {
  .assert_built(term)
  if (is.null(coef)) return(term@blueprint$value)
  term_refresh(term, coef)@blueprint$value
}

S7::method(term_predict, SegTerm) <- function(term, newdata, ...) {
  .assert_built(term)
  bp <- term@blueprint
  xv <- .seg_x(bp$var, newdata)
  grp <- if (is.null(bp$by)) NULL else .seg_by(bp$by, newdata, bp$levels)
  X <- .seg_assemble(bp, xv, grp, bp$coef)$X
  colnames(X) <- term@coef_names
  X
}

#' The Break-Points of a Segmented Term
#'
#' @description
#' The estimated positions of the break-points, one per break-point and
#' per level of \code{by}. For a continuous term these are coefficients;
#' for a discontinuous one they are read off two coefficients as
#' \eqn{-g/\kappa}, so this is the function that reports them either way.
#'
#' @param term A built \code{\link{SegTerm}}.
#' @param coef Optional coefficients; defaults to the ones the term was
#'   last refreshed at.
#'
#' @return A numeric vector of break-point positions.
#'
#' @examples
#' set.seed(1)
#' dd <- data.frame(x = sort(runif(100, 0, 10)))
#' seg_psi(term_build(seg(x, npsi = 2), dd))
#'
#' @seealso \code{\link{seg}}
#' @export
seg_psi <- function(term, coef = NULL) {
  if (!S7::S7_inherits(term, SegTerm)) {
    stop("'term' must be a segmented term.", call. = FALSE)
  }
  .assert_built(term)
  if (is.null(coef)) return(term@blueprint$psi)
  term_refresh(term, coef)@blueprint$psi
}

S7::method(print, SegTerm) <- function(x, ...) {
  if (term_is_built(x)) {
    cat(sprintf("<SegTerm> '%s': %s, %d break-point%s%s\n", x@label, x@kind,
                x@npsi, if (x@npsi == 1L) "" else "s",
                if (x@penalty_kind != "none")
                  sprintf("; %s on the changes", x@penalty_kind) else ""))
    cat("  at: ", paste(format(x@blueprint$psi, digits = 4),
                        collapse = ", "), "\n", sep = "")
  } else {
    cat(sprintf("<SegTerm> '%s': %s (specification)\n", x@label, x@kind))
  }
  invisible(x)
}
