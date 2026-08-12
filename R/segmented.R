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
    penalty_kind = S7::class_any,
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
#' The weight is unbounded as \eqn{x} approaches \eqn{\psi}, and since
#' \eqn{Z - \psi W} is the step itself, \eqn{Z} is \eqn{\psi W} plus a
#' quantity of order one: an unbounded \eqn{W} makes the two columns
#' numerically collinear and drowns the signal the fit reads. The remedy
#' of \cite{fasola2018} is to move the observations rather than to cap
#' the weight. With a scaling factor \eqn{c} the two intervals
#' \eqn{[x_{(1)}, \psi]} and \eqn{(\psi, x_{(n)}]} are mapped onto
#' \deqn{[x_{(1)},\, \psi - c(\psi - x_{(1)})]
#'   \quad\text{and}\quad
#'   (\psi + c(x_{(n)} - \psi),\, x_{(n)}],}
#' which leaves a gap of relative width \eqn{c} around \eqn{\psi} and
#' bounds \eqn{W} without altering the model: the working covariates are
#' computed on the rescaled covariate, while the truncated line, the
#' linear column and the reported contribution stay on the original one.
#'
#' The factor is not a constant. It governs how far the break-point may
#' travel in one step, so a large value lets the estimate leave a
#' spurious optimum and a small one is faithful to the step function.
#' \code{c0} is its starting value, and \code{\link{term_refresh}}
#' halves it whenever the break-point reverses direction, which is the
#' signal that the iteration has begun to circle an optimum rather than
#' travel towards one. The run has converged when the change in every
#' break-point falls below a hundredth of the smallest distance between
#' distinct observations, which \code{\link{seg_converged}} reports.
#' }
#'
#' \subsection{What the term carries}{
#' With \code{linear = TRUE} the block carries the linear effect too, so
#' the term is the whole relationship rather than the change in it.
#' \code{by} gives an independent set of break-points and changes per
#' level of a factor. \code{penalty} puts a penalty on the changes
#' themselves -- the slope changes for \code{seg}, the jump sizes for
#' \code{jump}, both for \code{jseg} -- and leaves the linear effect and
#' the break-points alone; with \code{"lasso"} that is a selection of how
#' many break-points are really there.
#'
#' The penalty is declared through \code{\link{term_penalties}}, which
#' names the coefficients it covers, rather than attached to the whole
#' block through a map that selects them. The two describe the same
#' function of the same coefficients and are not interchangeable to a
#' fitting layer: a separable penalty under a selection map is the
#' generalized-lasso problem, whose proximal operator does not split by
#' coordinate, so \code{\link[penalties7]{has_prox}} is \code{FALSE} for
#' it and neither a proximal step nor a coordinate descent can be taken.
#' Named as coordinates the map is the identity and both are available
#' unchanged. \code{jseg} declares two penalties, one over the slope
#' changes and one over the jump sizes, since a change of slope and a
#' change of level are not comparable quantities and cannot share a
#' hyperparameter.
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
#' The objective has local optima in the break-points, and the scaling
#' schedule widens the basin the iteration converges from rather than
#' removing the problem. Where the run begins therefore decides what it
#' finds, and \code{\link{seg_start}} is the answer: it scores an
#' equally spaced grid on the least-squares profile and returns the
#' specification with \code{psi} set to the best of it, which is what
#' \cite{fasola2018} recommend and what measurement supports over both
#' a conventional single start and bootstrap restarting.
#' A continuous term has no scaling factor, its working block being
#' bounded already; where its iteration alternates between two values
#' the remedy is to shrink the increment, as \code{segmented}'s
#' \code{h} does.
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
#' @param penalty The penalty on the changes: \code{"none"} (default),
#'   \code{"lasso"}, \code{"ridge"}, a \pkg{penalties7} penalty over as many
#'   coefficients as there are changes, or a function of that count
#'   returning one. A joint term declares two penalties, one on the slope
#'   changes and one on the jumps, and a penalty given as an object is used
#'   for both.
#' @param c0 For a discontinuous term, the starting value of the scaling
#'   factor that separates the observations from the break-point, as a
#'   fraction of the distance to the ends of the range. Defaults to
#'   \code{0.05}, the value \cite{fasola2018} recommend; see Details.
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
#' @seealso \code{\link{seg_psi}}, \code{\link{seg_start}}, \code{\link{seg_step}}, \code{\link{nl}}
#' @export
seg <- function(x, npsi = 1, psi = NULL, by = NULL, linear = TRUE,
                penalty = "none", c0 = 0.05,
                label = "seg") {
  .seg_spec("seg", substitute(x), npsi, psi, substitute(by), linear,
            .penalty_arg(penalty), c0, label)
}

#' @rdname seg
#' @export
jump <- function(x, npsi = 1, psi = NULL, by = NULL, linear = TRUE,
                 penalty = "none", c0 = 0.05,
                 label = "jump") {
  .seg_spec("jump", substitute(x), npsi, psi, substitute(by), linear,
            .penalty_arg(penalty), c0, label)
}

#' @rdname seg
#' @export
jseg <- function(x, npsi = 1, psi = NULL, by = NULL, linear = TRUE,
                 penalty = "none", c0 = 0.05,
                 label = "jseg") {
  .seg_spec("jseg", substitute(x), npsi, psi, substitute(by), linear,
            .penalty_arg(penalty), c0, label)
}

.seg_spec <- function(kind, var, npsi, psi, by, linear, penalty, c0,
                      label) {
  if (!is.numeric(c0) || length(c0) != 1L || is.na(c0) ||
      c0 <= 0 || c0 >= 1) {
    stop("'c0' must be a single number strictly between 0 and 1.",
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
          spec = list(psi = psi, c0 = c0),
          X = NULL, coef_names = character(0),
          blueprint = list(), penalty = NULL)
}

# the coefficient names of one level, and which of them are the changes a
# penalty may reach, grouped by what the change is: a slope change and a
# jump are different quantities, so a term carrying both carries two
# penalties rather than one over their union
.seg_names <- function(kind, npsi, linear) {
  k <- seq_len(npsi)
  nm <- if (linear) "lin" else character(0)
  groups <- list()
  if (kind %in% c("seg", "jseg")) {
    nm <- c(nm, paste0("delta", k))
    groups$delta <- paste0("delta", k)
  }
  if (kind == "seg") {
    nm <- c(nm, paste0("psi", k))
  } else {
    nm <- c(nm, paste0("kappa", k), paste0("g", k))
    groups$kappa <- paste0("kappa", k)
  }
  list(names = nm, changes = unlist(groups, use.names = FALSE),
       groups = groups)
}

# The break-points implied by one level's coefficients.
#
# A continuous term carries them as coefficients. A pure jump reads them
# off two, psi = -g/kappa, after Fasola et al. A JOINT term needs more
# than that, because the truncated line depends on the break-point as
# well and reading only the jump pair discards what the slope change
# says about it. Linearizing both parts about the previous position, and
# writing h for the increment, U for the truncated line and
# I = Z - psi0 W for the indicator,
#
#   delta (x-psi)_+  ~  delta U - delta h I
#   kappa 1(x>psi)    =  kappa Z - kappa psi W
#
# so the fitted coefficients of Z and W are a = kappa - delta h and
# b = delta h psi0 - kappa psi0 - kappa h, and the increment solves
#
#   delta h^2 + a h + (b + a psi0) = 0,
#
# the root of smaller modulus. At delta = 0 the quadratic degenerates to
# h = -(b + a psi0)/a, that is psi = -b/a, so the pure jump is the case
# the general form contains rather than an exception to it.
#
# The result is held inside `lim`, and that is not cosmetic. A
# break-point outside the data leaves the indicator constant, and then
# the linear effect, the truncated line and that constant are linearly
# dependent: the working block goes exactly singular. Confining the
# break-point to the interval where data lie on both sides is what the
# segmented literature does for the same reason.
.seg_psi_of <- function(kind, cf, npsi, linear, lim, psi_prev = NULL) {
  off <- if (linear) 1L else 0L
  psi <- if (kind == "seg") {
    cf[off + npsi + seq_len(npsi)]
  } else {
    d <- if (kind == "jseg") npsi else 0L
    a <- cf[off + d + seq_len(npsi)]
    b <- cf[off + d + npsi + seq_len(npsi)]
    a[abs(a) < 1e-12] <- 1e-12
    if (kind == "jseg" && !is.null(psi_prev)) {
      del <- cf[off + seq_len(npsi)]
      cc <- b + a * psi_prev
      disc <- a^2 - 4 * del * cc
      # the smaller root where the quadratic is one, the linear solution
      # where delta vanishes or the roots are complex
      lin_h <- -cc / a
      r <- sqrt(pmax(disc, 0))
      h1 <- (-a + r) / (2 * del)
      h2 <- (-a - r) / (2 * del)
      h <- ifelse(abs(h1) <= abs(h2), h1, h2)
      ok <- disc >= 0 & abs(del) > 1e-10 & is.finite(h)
      psi_prev + ifelse(ok, h, lin_h)
    } else {
      -b / a
    }
  }
  pmin(pmax(psi, lim[1L]), lim[2L])
}

# The rescaled covariate of Fasola, Muggeo and Kuchenhoff: the two
# intervals on either side of the break-point are shrunk towards the
# ends of the range, leaving a gap of width c around psi.
.seg_rescale <- function(xv, psi, cs, lo, hi) {
  keep <- 1 - cs
  ifelse(xv <= psi, lo + (xv - lo) * keep,
         (psi + cs * (hi - psi)) + (xv - psi) * keep)
}

# The working block and the true contribution of one level. The
# arithmetic is elementwise and is compiled (src/seg_block.cpp), which
# measures 1.2 to 3.2 times the R form below over n from 1e3 to 1e6 and
# one to five break-points. The operations are the same in the same
# order, so the two agree to a rounding; .seg_block_r is kept as the
# twin the tests compare against.
.seg_block <- function(kind, xv, cf, npsi, linear, cs, lo, hi, lim,
                       psi_prev = NULL) {
  psi <- .seg_psi_of(kind, cf, npsi, linear, lim, psi_prev)
  off <- if (linear) 1L else 0L
  del <- if (kind %in% c("seg", "jseg")) cf[off + seq_len(npsi)] else numeric(0)
  d <- if (kind == "jseg") npsi else 0L
  kap <- if (kind == "seg") numeric(0) else cf[off + d + seq_len(npsi)]
  out <- seg_block_cpp(switch(kind, seg = 0L, jump = 1L, jseg = 2L),
                       xv, psi, del, kap, cs,
                       if (linear) cf[1L] else 0, linear, lo, hi)
  list(X = out$X, value = out$value, psi = psi)
}

.seg_block_r <- function(kind, xv, cf, npsi, linear, cs, lo, hi, lim,
                         psi_prev = NULL) {
  n <- length(xv)
  off <- if (linear) 1L else 0L
  psi <- .seg_psi_of(kind, cf, npsi, linear, lim, psi_prev)
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
      # the break-point once the weight is frozen. It is applied to the
      # rescaled covariate, which is what keeps the weight bounded.
      xs <- .seg_rescale(xv, psi[j], cs[j], lo, hi)
      W <- 1 / (2 * abs(xs - psi[j]))
      cols[[length(cols) + 1L]] <- xs * W + 0.5
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
.seg_assemble <- function(bp, xv, grp, coef, cscale = bp$cscale) {
  per <- bp$per_level
  m <- length(bp$levels)
  n <- length(xv)
  X <- matrix(0, n, m * per)
  value <- numeric(n)
  psi <- numeric(0)
  for (l in seq_len(m)) {
    idx <- (l - 1L) * per + seq_len(per)
    jj <- (l - 1L) * bp$npsi + seq_len(bp$npsi)
    cs <- cscale[jj]
    rows <- if (is.null(grp)) rep(TRUE, n) else grp == bp$levels[l]
    b <- .seg_block(bp$kind, xv[rows], coef[idx], bp$npsi, bp$linear,
                    cs, bp$lo, bp$hi, bp$lim, bp$psi[jj])
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

  rr <- range(xv)
  if (diff(rr) <= 0) {
    stop("the covariate of a segmented term must vary.", call. = FALSE)
  }

  # the interval a break-point is held in: far enough inside the data
  # that both sides carry observations
  lim <- as.numeric(stats::quantile(xv, c(0.05, 0.95), names = FALSE))

  # The convergence tolerance of Fasola et al., a hundredth of the
  # distance between consecutive distinct observations: below it the
  # objective, a step function, cannot change. They write the SMALLEST
  # such distance, which is the same number on the evenly spaced
  # covariates of their examples and is of order n^-2 on a random one,
  # so the median is taken instead and the two rules agree wherever
  # theirs is usable.
  ux <- sort(unique(xv))
  delta <- 0.01 * if (length(ux) > 1L) stats::median(diff(ux)) else diff(rr)

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

  npt <- length(levs) * term@npsi
  bp <- list(kind = term@kind, npsi = term@npsi, linear = term@linear,
             per_level = per, levels = levs,
             lo = rr[1L], hi = rr[2L], lim = lim, delta = delta,
             cscale = rep(term@spec$c0, npt), sgn = rep(0, npt),
             step = rep(NA_real_, npt), nref = 0L,
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

  # The penalty names the coefficients it covers instead of selecting them
  # with a map: over its own coordinates it is separable, and a fitting
  # layer can take a proximal step or a coordinate descent on it, which a
  # selection map would deny it. One penalty per kind of change, shared
  # across the levels of `by`.
  bp$penalties <- list()
  if (!.penalty_is_none(term@penalty_kind)) {
    factory <- .penalty_factory(term@penalty_kind)
    for (g in names(nmi$groups)) {
      keep <- which(rep(nmi$names %in% nmi$groups[[g]], times = length(levs)))
      bp$penalties[[length(bp$penalties) + 1L]] <- list(
        name = g, index = keep, penalty = factory(length(keep)))
    }
  }

  bp$value <- asm$value
  bp$psi <- asm$psi
  term@X <- X
  term@coef_names <- cn
  term@blueprint <- bp
  term@penalty <- NULL
  term
}

#' @title Penalties of a Segmented Term
#' @name term_penalties.SegTerm
#' @description
#' One entry per kind of change the term carries a penalty on, naming the
#' coefficients it covers: \code{"delta"} for the slope changes,
#' \code{"kappa"} for the jump sizes, shared across the levels of
#' \code{by}. The list is empty when \code{penalty = "none"}, and for a
#' specification, whose parameters there is nothing yet to index: a penalty
#' is attached at build, as it is for every penalized term here.
#' @param term A built \code{\link{SegTerm}}.
#' @param ... Unused.
#' @return A list of entries, as \code{\link{term_penalties}} documents.
#' @keywords internal
S7::method(term_penalties, SegTerm) <- function(term, ...) {
  pens <- term@blueprint$penalties
  if (is.null(pens)) list() else pens
}

# the break-points of every level, read off the coefficients without
# building the block
.seg_psi_all <- function(bp, coef) {
  per <- bp$per_level
  unlist(lapply(seq_along(bp$levels), function(l) {
    .seg_psi_of(bp$kind, coef[(l - 1L) * per + seq_len(per)], bp$npsi,
                bp$linear, bp$lim,
                bp$psi[(l - 1L) * bp$npsi + seq_len(bp$npsi)])
  }), use.names = FALSE)
}

S7::method(term_refresh, SegTerm) <- function(term, coef, ...) {
  .assert_built(term)
  bp <- term@blueprint
  coef <- as.numeric(coef)
  if (length(coef) != ncol(term@X)) {
    stop(sprintf("'coef' must have length %d.", ncol(term@X)), call. = FALSE)
  }

  # The scaling schedule of Fasola et al.: the factor is halved whenever
  # the break-point reverses direction, which is the signal that the
  # iteration has begun to circle an optimum rather than travel towards
  # one. A large factor moves the estimate far and lets it leave a
  # spurious optimum; a small one is faithful to the step function.
  psi_new <- .seg_psi_all(bp, coef)
  d <- psi_new - bp$psi
  s <- sign(d)
  flip <- s != 0 & bp$sgn != 0 & s != bp$sgn
  bp$cscale[flip] <- bp$cscale[flip] / 2
  # The factor is not allowed to collapse. Writing D and d for the
  # distances from the break-point to the further and the nearer end of
  # the range, the weight spans D/(c d) and Z = psi W + 1(x > psi) shears
  # the pair by a further |psi|/D, so the condition number of the block
  # is of order max(D, |psi|)/(c d). Holding that below eps^-1/2, which
  # leaves half the digits of a double to a QR of the design, is
  #     c >= sqrt(eps) * max(D, |psi|) / d.
  # The bound is on the DESIGN: a caller forming the normal equations
  # squares it and needs a factor a thousand times larger, which is why
  # the working model is fitted by a QR of X, as `segmented` does.
  dfar <- pmax(psi_new - bp$lo, bp$hi - psi_new)
  dnear <- pmin(psi_new - bp$lo, bp$hi - psi_new)
  cmin <- sqrt(.Machine$double.eps) * pmax(dfar, abs(psi_new)) / dnear
  bp$cscale <- pmax(bp$cscale, cmin)
  bp$sgn[s != 0] <- s[s != 0]
  # The first refresh evaluates the block at the starting coefficients,
  # so the break-point has not moved and the difference is zero by
  # construction rather than because the iteration has finished.
  bp$nref <- bp$nref + 1L
  bp$step <- if (bp$nref > 1L) abs(d) else rep(NA_real_, length(d))

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

S7::method(term_value, SegTerm) <- function(term, coef = NULL, newdata = NULL,
                                            ...) {
  .assert_built(term)
  bp <- term@blueprint
  if (!is.null(newdata)) {
    # the break-points are the ones the term carries, as term_predict()
    # reapplies them: refreshing here would read them off the new rows
    xv <- .seg_x(bp$var, newdata)
    grp <- if (is.null(bp$by)) NULL else .seg_by(bp$by, newdata, bp$levels)
    cf <- if (is.null(coef)) bp$coef else as.numeric(coef)
    return(.seg_assemble(bp, xv, grp, cf)$value)
  }
  if (is.null(coef)) return(bp$value)
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

#' The Progress of a Break-Point Iteration
#'
#' @description
#' \code{seg_step} returns how far each break-point moved at the last
#' call to \code{\link{term_refresh}}, and \code{seg_converged} compares
#' the largest of those with the tolerance of \cite{fasola2018}, a
#' hundredth of the distance between consecutive distinct observations of
#' the covariate. A term that has been built but not yet refreshed has
#' taken no step, so \code{seg_step} returns \code{NA} and
#' \code{seg_converged} returns \code{FALSE}.
#'
#' @details
#' With \eqn{x_{(1)} < \cdots < x_{(m)}} the distinct covariate values, the
#' run stops at
#'
#' \deqn{\max_k \lvert \psi_k^{(t)} - \psi_k^{(t-1)} \rvert < \Delta,
#'   \qquad \Delta = 0.01 \cdot
#'     \operatorname{median}_{i}\, (x_{(i+1)} - x_{(i)}).}
#'
#' \cite{fasola2018} take the smallest of those gaps rather than their
#' median, which agrees with this on the evenly spaced covariates of their
#' examples and is of order \eqn{m^{-2}} on a random one, hence unreachable.
#'
#' The rule is one of resolution: below that distance the objective of a
#' discontinuous term, a step function of the break-point, cannot change.
#' It therefore tightens as the sample grows while the precision the
#' fixed point is reached at does not, so on a large sample the last
#' iterations move the break-point by a little more than the rule allows
#' and the run continues past the point where the estimate has settled.
#' A caller that can evaluate the objective should stop on its relative
#' change instead, which is what \code{segmented} does and what costs a
#' continuous term nothing: its iteration can settle into a cycle of
#' period two in the break-point, in which case this rule is never met
#' while the objective has long since stopped moving.
#'
#' @param term A built \code{\link{SegTerm}}.
#'
#' @return \code{seg_step} returns a numeric vector with one entry per
#'   break-point and per level of \code{by}; \code{seg_converged} returns
#'   a single logical.
#'
#' @examples
#' set.seed(1)
#' dd <- data.frame(x = sort(runif(100, 0, 10)))
#' dd$y <- 2 * (dd$x > 6) + rnorm(100, sd = 0.3)
#' b <- term_build(jump(x, psi = 4, linear = FALSE), dd)
#' cf <- b@blueprint$coef
#' for (it in 1:30) {
#'   b <- term_refresh(b, cf)
#'   X <- term_matrix(b)
#'   cf <- as.numeric(qr.solve(crossprod(X), crossprod(X, dd$y)))
#'   if (seg_converged(b)) break
#' }
#' c(psi = seg_psi(b, cf), step = seg_step(b))
#'
#' @seealso \code{\link{seg}}, \code{\link{seg_psi}}, \code{\link{seg_start}}
#'   \code{\link{seg_start}}
#' @export
seg_step <- function(term) {
  if (!S7::S7_inherits(term, SegTerm)) {
    stop("'term' must be a segmented term.", call. = FALSE)
  }
  .assert_built(term)
  term@blueprint$step
}

#' @rdname seg_step
#' @export
seg_converged <- function(term) {
  st <- seg_step(term)
  !anyNA(st) && max(st) < term@blueprint$delta
}

#' @title Whether a Segmented Term's Break-Points Have Settled
#' @name term_converged.SegTerm
#' @description
#' \code{\link{seg_converged}}, so that a fitting layer reads the rule the
#' construction is stopped on without knowing it is holding a break-point
#' term. A term that has not been refreshed has not moved and reports
#' \code{FALSE}, the first refresh being the one that measures nothing.
#' @param term A built \code{\link{SegTerm}}.
#' @param ... Unused.
#' @return A single logical.
#' @keywords internal
S7::method(term_converged, SegTerm) <- function(term, ...) {
  isTRUE(seg_converged(term))
}

#' Starting Positions for a Break-Point Term
#'
#' @description
#' Chooses the starting positions of a \code{\link{seg}}, \code{\link{jump}}
#' or \code{\link{jseg}} term by scoring an equally spaced grid on the
#' least-squares profile of the term's own columns, and returns the
#' specification with \code{psi} set to the best combination found.
#'
#' @details
#' \cite{fasola2018} recommend fixing the starting value by evaluating
#' the objective on a small grid spanned over the range of the covariate
#' rather than at a single conventional point, and the recommendation
#' matters more than it sounds: the objective has local optima in the
#' break-point, and the iteration converges from within a basin around
#' the position it starts at. Measured on a joint jump and change of
#' slope in 500 observations, over eight samples, the fraction of runs
#' recovering the break-point is 0 to 0.5 depending on where a single
#' start is placed and 1 from the grid.
#'
#' Writing \eqn{X(\psi)} for the design the term produces at a candidate
#' position, the position chosen is
#'
#' \deqn{\hat\psi = \arg\min_{\psi \in \mathcal{G}}
#'   \bigl\lVert y - X(\psi)\,
#'     \widehat{\beta}(\psi) \bigr\rVert^{2},
#'   \qquad \widehat{\beta}(\psi)
#'     = \arg\min_{\beta} \lVert y - X(\psi)\beta \rVert^{2},}
#'
#' over an equally spaced grid \eqn{\mathcal{G}} of \code{k} points in the
#' range of the covariate. The inner minimization is a linear fit, so the
#' whole rule costs \code{k} of them.
#'
#' The grid is scored on the residual sum of squares of an intercept,
#' the term's columns at each candidate position and, where the term
#' carries one, the linear effect. That is the exact profile for a
#' gaussian response and an adequate starting rule for any other, the
#' quantity being used to place a starting value and not to fit. With
#' \code{by} each level is scored on its own rows. With several
#' break-points every increasing combination of grid points is scored,
#' so \code{k} should be kept small.
#'
#' @param spec An unbuilt \code{\link{SegTerm}}.
#' @param data A data frame in which the covariate is evaluated.
#' @param y The response, one value per row of \code{data}.
#' @param k The number of grid points. Defaults to 10.
#'
#' @return The specification, with \code{psi} set.
#'
#' @references
#' Fasola, S., Muggeo, V. M. R. and Kuchenhoff, H. (2018). A heuristic,
#' iterative algorithm for change-point detection in abrupt change
#' models. \emph{Computational Statistics}, 33, 997--1015.
#'
#' @examples
#' set.seed(1)
#' dd <- data.frame(x = sort(runif(200, 0, 10)))
#' dd$y <- 0.3 * dd$x + 2 * (dd$x > 6.5) + rnorm(200, sd = 0.3)
#' seg_start(jump(x), dd, dd$y)@spec$psi
#'
#' @seealso \code{\link{seg}}
#' @export
seg_start <- function(spec, data, y, k = 10) {
  if (!S7::S7_inherits(spec, SegTerm)) {
    stop("'spec' must be a segmented term.", call. = FALSE)
  }
  if (length(spec@X) > 0L) {
    stop("'spec' must be an unbuilt term; see term_build().", call. = FALSE)
  }
  if (!is.numeric(k) || length(k) != 1L || is.na(k) || k < spec@npsi + 1) {
    stop(sprintf("'k' must be a single number of at least %d.",
                 spec@npsi + 1L), call. = FALSE)
  }
  xv <- .seg_x(spec@var, data)
  y <- as.numeric(y)
  if (length(y) != length(xv)) {
    stop("'y' must have one value per row of 'data'.", call. = FALSE)
  }
  grp <- .seg_by(spec@by, data)
  levs <- if (is.null(grp)) "" else levels(grp)

  lim <- as.numeric(stats::quantile(xv, c(0.05, 0.95), names = FALSE))
  grid <- seq(lim[1L], lim[2L], length.out = as.integer(k))
  combos <- utils::combn(length(grid), spec@npsi, simplify = FALSE)

  cols <- function(x, psi) {
    Z <- if (spec@linear) cbind(1, x) else matrix(1, length(x), 1L)
    if (spec@kind %in% c("seg", "jseg")) {
      for (p in psi) Z <- cbind(Z, pmax(x - p, 0))
    }
    if (spec@kind %in% c("jump", "jseg")) {
      for (p in psi) Z <- cbind(Z, as.numeric(x > p))
    }
    Z
  }

  best <- unlist(lapply(levs, function(l) {
    rows <- if (is.null(grp)) rep(TRUE, length(xv)) else grp == l
    xs <- xv[rows]; ys <- y[rows]
    v <- vapply(combos, function(ii) {
      Z <- cols(xs, grid[ii])
      out <- tryCatch(sum(qr.resid(qr(Z), ys)^2), error = function(e) Inf)
      if (is.finite(out)) out else Inf
    }, numeric(1))
    grid[combos[[which.min(v)]]]
  }), use.names = FALSE)

  # one set of positions is carried, so the levels are averaged when
  # `by` splits them; the per-level fit moves them apart from there
  if (length(levs) > 1L) {
    best <- rowMeans(matrix(best, nrow = spec@npsi))
  }
  spec@spec$psi <- best
  spec
}

S7::method(print, SegTerm) <- function(x, ...) {
  if (term_is_built(x)) {
    cat(sprintf("<SegTerm> '%s': %s, %d break-point%s%s\n", x@label, x@kind,
                x@npsi, if (x@npsi == 1L) "" else "s",
                if (!.penalty_is_none(x@penalty_kind))
                  sprintf("; %s on the changes",
                          .penalty_label(x@penalty_kind)) else ""))
    cat("  at: ", paste(format(x@blueprint$psi, digits = 4),
                        collapse = ", "), "\n", sep = "")
  } else {
    cat(sprintf("<SegTerm> '%s': %s (specification)\n", x@label, x@kind))
  }
  invisible(x)
}
