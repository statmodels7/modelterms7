#' @include segmented.R block_deriv3.R block_deriv4.R
NULL

#' @title The Third and Fourth Derivatives of a Break-Point Term's Block
#' @name term_block_deriv3.SegTerm
#'
#' @description
#' Return the third ([term_block_deriv3()]) and fourth ([term_block_deriv4()])
#' derivatives of the block of a smoothed [seg()], [jump()] or [jseg()] along
#' three and four directions. A sharp construction returns zeros, for the
#' reasons that [term_block_deriv2()] records for its second derivative.
#'
#' @details
#' With \eqn{u = x - \psi}, the smoother \eqn{s} of the term and
#' \eqn{d_j = s^{(j)}/2}, the columns of a smoothed term are \eqn{H(u)} for a
#' change \eqn{\gamma}, \eqn{S(u)} for a step \eqn{\delta} and
#' \eqn{-(\gamma S(u) + \delta d_2(u))} for a break-point, where
#' \eqn{H' = S = 1/2 + d_1}. Write \eqn{p_i} for the movement of the
#' break-point along direction \eqn{i} (its development's design times the
#' direction, gated where the position sits against a confinement limit) and
#' \eqn{g_i}, \eqn{e_i} for the movements of the change and of the step. A
#' derivative in a break-point coefficient multiplies a function of \eqn{u} by
#' \eqn{-p_i} and raises its order by one, and the change and the step enter
#' the break-point column linearly. At order three the columns are
#' \deqn{-d_3\,p_ap_bp_c,\qquad -d_4\,p_ap_bp_c,}
#' \deqn{(\gamma d_4 + \delta d_5)\,p_ap_bp_c - (g_ad_3 + e_ad_4)\,p_bp_c
#'   - (g_bd_3 + e_bd_4)\,p_ap_c - (g_cd_3 + e_cd_4)\,p_ap_b,}
#' and at order four
#' \deqn{d_4\,p_ap_bp_cp_e,\qquad d_5\,p_ap_bp_cp_e,}
#' \deqn{-(\gamma d_5 + \delta d_6)\,p_ap_bp_cp_e + \sum_i (g_id_4 + e_id_5)
#'   \prod_{j\ne i}p_j.}
#' The fourth derivative of a term with a step reads the sixth derivative of
#' the smoother, which the smoothers of numericals7 do not carry; the method
#' signals an error there.
#'
#' @param term A built [SegTerm()].
#' @param coef The coefficients, or `NULL` for the stored ones.
#' @param v,u,s,r Numeric vectors as long as the term's coefficients.
#' @param ... Unused.
#'
#' @return A numeric matrix, one row per observation and one column per
#'   coefficient of the term.
#'
#' @examples
#' set.seed(1)
#' d <- data.frame(x = sort(runif(40, 0, 10)))
#' b <- term_build(seg(x, smoothed = numericals7::smooth_probit()), d)
#' m <- length(term_coef_names(b))
#' dim(term_block_deriv3(b, v = rep(1, m), u = rep(1, m), s = rep(1, m)))
#'
#' @aliases term_block_deriv4.SegTerm
#' @keywords internal
S7::method(term_block_deriv3, SegTerm) <- function(term, coef = NULL, v, u, s,
                                                   ...) {
  .seg_deriv_high(term, coef, list(v = v, u = u, s = s))
}

# documented with term_block_deriv3.SegTerm
S7::method(term_block_deriv4, SegTerm) <- function(term, coef = NULL, v, u, s,
                                                   r, ...) {
  .seg_deriv_high(term, coef, list(v = v, u = u, s = s, r = r))
}

# The shared body of the two methods: the checks, the zeros of a sharp
# construction and the closed forms of the smoothed one at order three or four.
.seg_deriv_high <- function(term, coef, dirs) {
  bp <- term@build_info
  if (!length(bp)) stop("the term is not built.", call. = FALSE)
  cf <- if (is.null(coef)) bp$coef else as.numeric(coef)
  ncoef <- length(term_coef_names(term))
  n <- length(bp$xv)
  for (nm in names(dirs)) {
    dirs[[nm]] <- as.numeric(dirs[[nm]])
    if (length(dirs[[nm]]) != ncoef) {
      stop(sprintf("'%s' must have length %d, the term's coefficients.", nm,
                   ncoef), call. = FALSE)
    }
  }
  out <- matrix(0, n, ncoef)
  if (is.null(bp$smooth)) return(out)
  order <- length(dirs)
  has_g <- bp$kind %in% c("seg", "jseg")
  has_d <- bp$kind %in% c("jump", "jseg")
  top <- order + if (has_d) 2L else 1L
  if (length(bp$smooth$sm@s) < top + 1L) {
    stop(sprintf(paste0("the derivative of order %d of this block reads the ",
                        "derivative of order %d\n  of the smoother, which it ",
                        "does not carry."), order, top), call. = FALSE)
  }
  along <- function(vec, p) {
    idx <- bp$index[[p]]
    z <- bp$Z[[p]]
    if (is.null(z)) rep(vec[idx], n) else as.numeric(as.matrix(z) %*% vec[idx])
  }
  put <- function(p, w) {
    idx <- bp$index[[p]]
    z <- bp$Z[[p]]
    out[, idx] <<- out[, idx] + if (is.null(z)) w else as.matrix(z) * w
  }
  val <- function(p) .seg_pval(list(Z = bp$Z, index = bp$index), cf, p, n)
  sm <- bp$smooth
  wv <- .seg_smooth_w(sm, n, bp$Z)
  pos <- .seg_positions(bp, cf, n)
  for (k in seq_len(bp$npsi)) {
    pk <- paste0("psi", k)
    gk <- paste0("gamma", k)
    dk <- paste0("delta", k)
    uu <- bp$xv - pos$psi[, k]
    zp <- bp$Z[[pk]]
    lv <- if (is.null(zp)) rep(pos$pk[[k]], n) else
      as.numeric(as.matrix(zp) %*% pos$pk[[k]])
    free <- as.numeric(lv > .seg_lim(bp, k)[1L] & lv < .seg_lim(bp, k)[2L])
    dj <- function(j) .seg_sw(sm$sm, uu, wv, j) / 2
    p <- lapply(dirs, function(x) free * along(x, pk))
    g <- if (has_g) lapply(dirs, function(x) along(x, gk))
    e <- if (has_d) lapply(dirs, function(x) along(x, dk))
    prod_p <- Reduce(`*`, p)
    # the product of the break-point movements without direction i
    rest <- function(i) Reduce(`*`, p[-i])
    if (order == 3L) {
      d3 <- dj(3L)
      d4 <- dj(4L)
      self <- numeric(n)
      side <- numeric(n)
      if (has_g) {
        put(gk, -d3 * prod_p)
        self <- self + val(gk) * d4
        for (i in 1:3) side <- side + g[[i]] * d3 * rest(i)
      }
      if (has_d) {
        put(dk, -d4 * prod_p)
        self <- self + val(dk) * dj(5L)
        for (i in 1:3) side <- side + e[[i]] * d4 * rest(i)
      }
      put(pk, self * prod_p - side)
    } else {
      d4 <- dj(4L)
      d5 <- dj(5L)
      self <- numeric(n)
      side <- numeric(n)
      if (has_g) {
        put(gk, d4 * prod_p)
        self <- self + val(gk) * d5
        for (i in 1:4) side <- side + g[[i]] * d4 * rest(i)
      }
      if (has_d) {
        put(dk, d5 * prod_p)
        self <- self + val(dk) * dj(6L)
        for (i in 1:4) side <- side + e[[i]] * d5 * rest(i)
      }
      put(pk, side - self * prod_p)
    }
  }
  out
}
