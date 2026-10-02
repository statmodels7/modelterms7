#' @include term_classes.R generics.R segmented.R
NULL

#' @title The Coefficients a Term Holds
#'
#' @description
#' The positions, among the term's own coefficients, of those the term holds
#' at their values: the fitting layer keeps them out of the solve, out of the
#' information it inverts and out of the reading of the mode. A term whose
#' coefficients are all estimated by the derivative of the objective holds
#' none, which is the base method's answer.
#'
#' @details
#' The one shipped term that answers otherwise is a [jump()] or [jseg()] term
#' put in its held state by [seg_hold()]. There the break-point is fixed at a
#' minimum of the profile objective, which is constant between consecutive
#' observations and so has no derivative to estimate the position with, and
#' the coefficient slot that carried it in the working construction becomes a
#' column of zeros.
#'
#' @param term A term.
#' @param ... Passed to methods.
#'
#' @return An integer vector of positions in the term's coefficients, possibly
#'   empty.
#'
#' @seealso [seg_hold()], [term_kinks()].
#'
#' @examples
#' set.seed(1)
#' d <- data.frame(x = sort(runif(100, 0, 10)))
#' d$y <- 1 + 2 * (d$x > 6) + rnorm(100, sd = 0.3)
#' b <- term_build(jump(x, psi = 6), d)
#' term_held(b)
#' term_held(seg_hold(b))
#'
#' @export
#' @aliases term_held.model_term
term_held <- S7::new_generic("term_held", "term",
  function(term, ...) S7::S7_dispatch())

S7::method(term_held, model_term) <- function(term, ...) integer(0)

#' @title The Coefficients a Held Break-Point Term Holds
#' @name term_held.SegTerm
#' @description
#' The coefficient slots of the break-points where the term is held by
#' [seg_hold()], and none otherwise.
#' @param term A [SegTerm()].
#' @param ... Unused.
#' @return An integer vector of positions in the term's coefficients.
#' @keywords internal
S7::method(term_held, SegTerm) <- function(term, ...) {
  bp <- term@blueprint
  if (!length(bp) || !isTRUE(bp$held)) return(integer(0))
  as.integer(unlist(bp$index[paste0("psi", seq_len(bp$npsi))],
                    use.names = FALSE))
}

#' Hold a Break-Point Term at Its Positions
#'
#' @description
#' Puts a sharp [jump()] or [jseg()] term in its held state, or takes it out.
#' A held term keeps its break-points where they are and its block is the
#' exact design at those positions: the covariate for the linear part, the
#' truncated line \eqn{(x - \psi_k)_+} for each change of slope, the indicator
#' \eqn{1(x > \psi_k)} for each change of level, and a column of zeros in the
#' slot that carried the position.
#'
#' @details
#' The working construction of Fasola, Muggeo and Kuchenhoff (2018) reaches
#' its break-points through a block whose weight is frozen at the previous
#' iterate, so at its fixed point the coefficients are those of the working
#' model and not the least-squares coefficients at the positions reached. With
#' the positions held the term contributes a linear function of its
#' remaining coefficients, and a fit of those is an ordinary fit, which is how
#' `segmented` and `stepmented` report their coefficients. [term_held()]
#' names the slot that is held, [term_jacobian_block()] answers `TRUE`, and
#' [term_refresh()] no longer moves the positions.
#'
#' Taking the term out of the held state writes the position back into the
#' slot as \eqn{g_k = -\delta_k\psi_k}, so the read-off of the working
#' construction returns the held positions and the iteration can continue
#' from them.
#'
#' A continuous [seg()] and a smoothed term are rejected: their positions are
#' ordinary coefficients with a column of their own.
#'
#' @param term A built [jump()] or [jseg()] term.
#' @param hold `TRUE` to hold the term, `FALSE` to release it.
#'
#' @return The term in the requested state, its block rebuilt.
#'
#' @references
#' Fasola, S., Muggeo, V. M. R. and Kuchenhoff, H. (2018). A heuristic,
#' iterative algorithm for change-point detection in abrupt change models.
#' *Computational Statistics*, 33, 997--1015.
#'
#' @seealso [term_held()], [seg_polish_exact()], [seg_relocate()].
#'
#' @examples
#' set.seed(1)
#' d <- data.frame(x = sort(runif(100, 0, 10)))
#' d$y <- 1 + 2 * (d$x > 6) + rnorm(100, sd = 0.3)
#' b <- seg_hold(term_build(jump(x, psi = 6), d))
#' term_matrix(b)[c(1, 100), ]
#' seg_psi(b)
#'
#' @export
seg_hold <- function(term, hold = TRUE) {
  if (!S7::S7_inherits(term, SegTerm)) {
    stop("'term' must be a break-point term.", call. = FALSE)
  }
  .assert_built(term)
  bp <- term@blueprint
  if (identical(bp$kind, "seg") || !is.null(bp$smooth)) {
    stop(paste("only a sharp jump() or jseg() term can be held: the",
               "break-points of a seg() or a smoothed term are ordinary",
               "coefficients."), call. = FALSE)
  }
  cf <- bp$coef
  if (!isTRUE(hold)) {
    if (!isTRUE(bp$held)) return(term)
    # the slot carries the position again, g = -delta psi, which is what the
    # read-off of the working construction inverts
    for (k in seq_len(bp$npsi)) {
      dk <- cf[bp$index[[paste0("delta", k)]]]
      cf[bp$index[[paste0("psi", k)]]] <- -dk * bp$pk[[k]]
    }
    bp$held <- NULL
  } else {
    bp$held <- TRUE
  }
  bp$floor <- NULL
  bp$step <- rep(NA_real_, bp$npsi)
  asm <- .seg_assemble(bp, bp$xv, cf)
  X <- asm$X
  colnames(X) <- term@coef_names
  bp$coef <- cf
  bp$value <- asm$value
  bp$psi <- asm$psi
  bp$pk <- asm$pk
  term@X <- X
  term@blueprint <- bp
  term
}

# The block of a held term: the exact design at the held positions, with the
# position's own slot a column of zeros. The positions are the stored ones,
# read through .seg_pk(), which returns them whatever the coefficients are.
.seg_block_held <- function(bp, xv, coef, Z = bp$Z) {
  n <- length(xv)
  K <- bp$npsi
  pos <- .seg_positions(bp, coef, n, Z, bp$pk)
  psi <- pos$psi
  value <- numeric(n)
  cols <- list()
  put <- function(p, mult) {
    z <- Z[[p]]
    cols[[p]] <<- if (is.null(z)) matrix(mult, ncol = 1L) else mult * z
  }
  val <- function(p) .seg_pval(list(Z = Z, index = bp$index), coef, p, n)
  if (bp$linear) {
    put("beta", xv)
    value <- value + val("beta") * xv
  }
  for (k in seq_len(K)) {
    if (bp$kind == "jseg") {
      tr <- pmax(xv - psi[, k], 0)
      put(paste0("gamma", k), tr)
      value <- value + val(paste0("gamma", k)) * tr
    }
    ind <- as.numeric(xv > psi[, k])
    put(paste0("delta", k), ind)
    value <- value + val(paste0("delta", k)) * ind
    put(paste0("psi", k), rep(0, n))
  }
  list(X = .nl_bind(cols[bp$params]), value = value, psi = psi, pk = pos$pk)
}

#' Polish a Break-Point Term's Positions Over Every Interval
#'
#' @description
#' The exact minimum of the least-squares profile of a sharp [jump()] or
#' [jseg()] term, one break-point at a time with the others held. The
#' profile of such a term is constant between consecutive observations of the
#' covariate, so its minimum over one break-point is found by evaluating it
#' once in each interval between consecutive distinct values, and this
#' function does that for every interval inside the confinement limits.
#'
#' @details
#' Write \eqn{F} for the columns that do not move with the break-point being
#' polished (an intercept, the covariate where the term carries a linear
#' part, and the columns of the other break-points), and \eqn{a(q)} for the
#' columns that do: the indicator \eqn{1(x > q)}, and for a [jseg()] also the
#' truncated line \eqn{(x - q)1(x > q)}. The residual sum of squares at
#' \eqn{q} is that of the response on \eqn{F} less the part \eqn{a(q)}
#' explains of what \eqn{F} leaves, and every cross product \eqn{a(q)} enters
#' is a sum over the observations above \eqn{q} of a fixed quantity, or of
#' one fixed quantity minus \eqn{q} times another. Sorted once, those sums
#' are cumulative sums, so all the intervals cost \eqn{O(np^2)} together
#' rather than one linear fit each. The sweeps over the break-points repeat
#' until none moves. Each break-point is placed at the midpoint of its
#' interval.
#'
#' [seg_polish()] sweeps a grid of `k` points instead, which reaches a
#' neighbourhood of the minimum and not the interval: measured on 16 samples
#' of 200 observations, a fit polished by the grid reached the interval of
#' the global minimum in 11.
#'
#' The profile is least squares of `y` on the term's own columns plus an
#' intercept, exact for a gaussian response; a fitting layer for any other
#' family accepts the positions only where its own objective improves. A
#' term whose per-break-point coefficients carry a development is rejected,
#' its positions being one per observation.
#'
#' @inheritParams seg_polish
#'
#' @return The term at the polished positions (see [seg_relocate()]).
#'
#' @seealso [seg_polish()], [seg_hold()], [seg_profile_rss()].
#'
#' @examples
#' set.seed(1)
#' dd <- data.frame(x = sort(runif(300, 0, 10)))
#' dd$y <- 2 * (dd$x > 3) - 1.5 * (dd$x > 7) + rnorm(300, sd = 0.3)
#' b <- term_build(jump(x, npsi = 2, psi = c(2, 5)), dd)
#' seg_psi(seg_polish_exact(b, dd$y))
#'
#' @export
seg_polish_exact <- function(term, y, sweeps = 10, weights = NULL) {
  ev <- .seg_interval_eval(term, y, weights)
  pr <- ev$pr
  mid <- ev$mid
  if (!length(mid)) return(term)
  K <- term@npsi
  psi <- as.numeric(term@blueprint$psi[1L, ])
  best <- pr$rss(psi)
  for (s in seq_len(as.integer(sweeps))) {
    moved <- FALSE
    for (j in seq_len(K)) {
      v <- ev$rss_at(j, psi)
      i <- which.min(v)
      if (length(i) && is.finite(v[i]) && v[i] < best - 1e-10 * (best + 1)) {
        psi[j] <- mid[i]
        best <- v[i]
        moved <- TRUE
      }
    }
    if (!moved) break
  }
  seg_relocate(term, psi)
}

#' The Profile of One Break-Point Over Every Interval
#'
#' @description
#' The least-squares profile of a sharp [jump()] or [jseg()] term in one of its
#' break-points, the others held where they are: the residual sum of squares
#' at the midpoint of every interval between consecutive distinct values of
#' the covariate inside the confinement limits. It is the quantity
#' [seg_polish_exact()] minimizes, returned whole so that a caller whose
#' objective the profile only approximates can take several candidates from
#' it.
#'
#' @inheritParams seg_polish
#' @param k Which break-point, an integer.
#'
#' @return A data frame with columns `psi`, the midpoints, and `rss`, the
#'   profile there (`Inf` where the design is singular).
#'
#' @seealso [seg_polish_exact()].
#'
#' @examples
#' set.seed(1)
#' dd <- data.frame(x = sort(runif(200, 0, 10)))
#' dd$y <- 1 + 1.5 * (dd$x > 6) + rnorm(200, sd = 0.4)
#' b <- term_build(jump(x, psi = 4), dd)
#' pr <- seg_profile_intervals(b, dd$y)
#' pr$psi[which.min(pr$rss)]
#'
#' @export
seg_profile_intervals <- function(term, y, k = 1L, weights = NULL) {
  ev <- .seg_interval_eval(term, y, weights)
  K <- term@npsi
  k <- as.integer(k)
  if (length(k) != 1L || is.na(k) || k < 1L || k > K) {
    stop(sprintf("'k' must be one of 1 to %d.", K), call. = FALSE)
  }
  if (!length(ev$mid)) return(data.frame(psi = numeric(0), rss = numeric(0)))
  psi <- as.numeric(term@blueprint$psi[1L, ])
  data.frame(psi = ev$mid, rss = ev$rss_at(k, psi))
}

# The machinery of the two: the midpoints and an evaluator of the profile in
# break-point j over every midpoint, the other positions given.
.seg_interval_eval <- function(term, y, weights = NULL) {
  pr <- .seg_profile(term, y, weights)
  bp <- term@blueprint
  if (identical(bp$kind, "seg") || !is.null(bp$smooth)) {
    stop("only a sharp jump() or jseg() term is polished over intervals.",
         call. = FALSE)
  }
  xv <- bp$xv
  n <- length(xv)
  w <- if (is.null(weights)) rep(1, n) else as.numeric(weights)
  yv <- as.numeric(y)
  o <- order(xv)
  xs <- xv[o]
  ws <- w[o]
  ys <- yv[o]
  u <- unique(xs)
  mid <- (u[-1L] + u[-length(u)]) / 2
  mid <- mid[mid > pr$lim[1L] & mid < pr$lim[2L]]
  # the first sorted observation above each candidate
  first <- findInterval(mid, xs) + 1L
  jseg <- identical(bp$kind, "jseg")
  K <- bp$npsi
  cols_of <- function(q) {
    if (jseg) cbind(pmax(xs - q, 0), as.numeric(xs > q))
    else matrix(as.numeric(xs > q), ncol = 1L)
  }
  # suffix sums over the sorted observations, read at `first`
  suffix <- function(v) {
    v <- as.matrix(v)
    cs <- apply(v, 2L, function(c) rev(cumsum(rev(c))))
    cs <- rbind(as.matrix(cs), 0)
    cs[first, , drop = FALSE]
  }
  rss_at <- function(j, psi) {
    Fm <- if (bp$linear) cbind(1, xs) else matrix(1, n, 1L)
    for (k in seq_len(K)) if (k != j) Fm <- cbind(Fm, cols_of(psi[k]))
    sw <- sqrt(ws)
    qrF <- qr(Fm * sw)
    if (qrF$rank < ncol(Fm)) return(rep(Inf, length(mid)))
    r <- as.numeric(qr.resid(qrF, ys * sw)) / ifelse(sw > 0, sw, 1)
    r[sw == 0] <- 0
    rss_f <- sum(ws * r^2)
    Mi <- chol2inv(qr.R(qrF))
    Sw <- suffix(ws)[, 1L]
    SwF <- suffix(ws * Fm)
    SwR <- suffix(ws * r)[, 1L]
    # the indicator: a'Wa, a'WF, a'Wr
    qI <- rowSums((SwF %*% Mi) * SwF)
    sII <- Sw - qI
    cI <- SwR
    if (!jseg) {
      out <- rss_f - cI^2 / sII
      out[!(sII > 1e-12 * pmax(Sw, 1))] <- Inf
      return(out)
    }
    Swx <- suffix(ws * xs)[, 1L]
    Swxx <- suffix(ws * xs^2)[, 1L]
    SwxF <- suffix(ws * xs * Fm)
    SwxR <- suffix(ws * xs * r)[, 1L]
    q <- mid
    # the truncated line (x - q) 1(x > q)
    TF <- SwxF - q * SwF
    tt <- Swxx - 2 * q * Swx + q^2 * Sw
    ti <- Swx - q * Sw
    sTT <- tt - rowSums((TF %*% Mi) * TF)
    sTI <- ti - rowSums((TF %*% Mi) * SwF)
    cT <- SwxR - q * SwR
    det <- sTT * sII - sTI^2
    out <- rss_f - (sII * cT^2 - 2 * sTI * cT * cI + sTT * cI^2) / det
    out[!(det > 1e-12 * pmax(sTT * sII, 1e-300))] <- Inf
    out
  }
  list(pr = pr, mid = mid, rss_at = rss_at)
}
