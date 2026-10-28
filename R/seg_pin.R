#' @include term_classes.R generics.R segmented.R
NULL

#' Confine the Free Break-Points of a Term Away From the Held Ones
#'
#' @description
#' Prepares a continuous [seg()] term for a fit in which some of its
#' coefficients are held at given values. Each break-point whose position is
#' not held stays in the interval between the held positions that surround it,
#' with at least `min_obs` observations between it and each of them, and the
#' coefficients of a break-point of which a coefficient is held are not
#' exchanged with those of another when the positions cross.
#'
#' @details
#' Two break-points of a continuous term at one position make two columns of
#' the block equal. As the free break-point approaches the held one, the
#' coefficients of the two changes of slope grow without bound in opposite
#' directions, because the pair then imitates a jump of the level. Without
#' this function the relabeling of crossed break-points in [term_refresh()]
#' exchanges the coefficients of the held break-point with those of the free
#' one, so that the value held for the first is written onto the second and
#' the fit alternates between two states without converging.
#'
#' A change of slope is estimated from the observations between its
#' break-point and the next one, so that two observations are the least that
#' identifies it, and the interval of a free break-point leaves that many
#' between it and a held one. The interval is the one, among those that the held
#' positions cut out of the confinement interval of the term, that contains
#' the position in `start` or, when that position leaves fewer observations
#' than that, the nearest interval with room. The side is therefore the one on
#' which the break-point starts, and the break-point cannot cross a held one
#' during the fit. The held positions are clamped to the confinement interval
#' of the term. Free break-points that share an interval exchange their
#' labels when their positions cross, as in an ordinary fit.
#'
#' The columns of a free break-point at an end of its interval are zero, as
#' they are at the confinement limits of the term, so that the coefficient
#' stays at the end.
#'
#' @param term A built [seg()] term whose break-point positions carry no
#'   development.
#' @param held A named numeric vector of the coefficients of the term that
#'   are held. The names are those of `term@coef_names`. The values of the
#'   coefficients that carry the positions of break-points are positions; the
#'   other values are not read, and the coefficient of the linear part does not
#'   affect the break-points.
#' @param start The positions of the break-points from which the fit starts,
#'   one for each break-point, ordered as the term orders them. The default is
#'   the positions of the term.
#' @param min_obs The smallest number of observations between a free and a
#'   held break-point, a positive integer.
#'
#' @return The term, with the intervals of its break-points and the groups
#'   within which they exchange labels recorded. [seg_limits()] returns the
#'   intervals.
#'
#' @seealso [seg_limits()], [term_refresh()], [seg_hold()].
#'
#' @examples
#' set.seed(1)
#' d <- data.frame(x = sort(runif(100, 0, 10)))
#' b <- term_build(seg(x, npsi = 2, psi = c(3, 7)), d)
#' p <- seg_pin(b, held = c(seg.psi2 = 5), start = c(3, 7))
#' seg_limits(p)
#'
#' @export
seg_pin <- function(term, held, start = NULL, min_obs = 2L) {
  if (!S7::S7_inherits(term, SegTerm)) {
    stop("'term' must be a break-point term.", call. = FALSE)
  }
  .assert_built(term)
  bp <- term@build_info
  if (!identical(bp$kind, "seg")) {
    stop("only a seg() term can be pinned: the positions of a jump() or a ",
         "jseg() are held by seg_hold().", call. = FALSE)
  }
  K <- bp$npsi
  pk <- paste0("psi", seq_len(K))
  if (any(vapply(pk, function(p) !is.null(bp$Z[[p]]), logical(1)))) {
    stop("the positions of a pinned term must carry no development.",
         call. = FALSE)
  }
  if (!is.numeric(held) || !length(held) || is.null(names(held)) ||
      !all(is.finite(held))) {
    stop("'held' must be a named numeric vector of finite positions.",
         call. = FALSE)
  }
  if (anyNA(match(names(held), term@coef_names)) || anyDuplicated(names(held))) {
    stop("'held' must name distinct coefficients of the term.", call. = FALSE)
  }
  nm_of <- function(prefix) {
    term@coef_names[vapply(paste0(prefix, seq_len(K)),
                           function(p) bp$index[[p]][1L], integer(1))]
  }
  at_pos <- match(names(held), nm_of("psi"))
  at_chg <- match(names(held), nm_of("gamma"))
  # the slots of which a coefficient is held, and the positions held
  slot <- which(!is.na(at_pos))
  hold_pos <- as.numeric(held)[slot]
  slot <- at_pos[slot]
  fixed <- sort(unique(c(slot, at_chg[!is.na(at_chg)])))
  if (!is.numeric(min_obs) || length(min_obs) != 1L || !is.finite(min_obs) ||
      min_obs < 1 || min_obs != round(min_obs)) {
    stop("'min_obs' must be a positive integer.", call. = FALSE)
  }
  min_obs <- as.integer(min_obs)
  if (is.null(start)) start <- as.numeric(bp$psi[1L, ])
  if (!is.numeric(start) || length(start) != K || !all(is.finite(start))) {
    stop(sprintf("'start' must give %d finite positions.", K), call. = FALSE)
  }

  lim <- bp$lim
  pos <- pmin(pmax(hold_pos, lim[1L]), lim[2L])
  hs <- sort(pos)
  # the observations that are active for a held break-point are those above
  # it by more than the tolerance that the block reads them with; the k-th one
  # above it is the nearest position of a free break-point to its right that
  # leaves k observations between the two, and the k-th one at or below it,
  # less the tolerance, the nearest to its left
  tie <- .seg_tie(bp)
  xs <- sort(bp$xv)
  right <- vapply(hs, function(h) {
    above <- xs[xs - h > tie]
    if (length(above) >= min_obs) above[min_obs] else Inf
  }, numeric(1))
  left <- vapply(hs, function(h) {
    below <- sort(xs[xs - h <= tie], decreasing = TRUE)
    if (length(below) >= min_obs) below[min_obs] - 2 * tie else -Inf
  }, numeric(1))
  # the intervals that the held positions cut out of the confinement interval
  lo <- pmax(c(lim[1L], right), lim[1L])
  hi <- pmin(c(left, lim[2L]), lim[2L])
  room <- lo < hi
  lims <- matrix(rep(lim, each = K), K, 2L)
  free <- setdiff(seq_len(K), slot)
  for (j in free) {
    if (!length(hs)) break
    if (!any(room)) {
      stop("the held positions leave no room for the free break-points.",
           call. = FALSE)
    }
    dist <- pmax(lo - start[j], start[j] - hi, 0)
    dist[!room] <- Inf
    i <- which.min(dist)
    lims[j, ] <- c(lo[i], hi[i])
  }
  group <- integer(K)
  group[fixed] <- K + seq_along(fixed)
  rest <- setdiff(seq_len(K), fixed)
  key <- paste(lims[rest, 1L], lims[rest, 2L])
  group[rest] <- match(key, unique(key))
  bp$lims <- lims
  bp$group <- group
  bp$pin <- list(held = fixed, position = pos, min_obs = min_obs)
  term@build_info <- bp
  term
}

#' Confinement Intervals of the Break-Points of a Term
#'
#' @description
#' Returns the interval in which each break-point of a built [seg()],
#' [jump()] or [jseg()] term can lie: the interval between the 5th and the
#' 95th percentile of the covariate for every break-point, or the narrower
#' interval that [seg_pin()] gave to a break-point that is free beside a held
#' one.
#'
#' @param term A built break-point term.
#'
#' @return A numeric matrix with one row for each break-point and the columns
#'   `lower` and `upper`.
#'
#' @seealso [seg_pin()].
#'
#' @examples
#' set.seed(1)
#' d <- data.frame(x = sort(runif(100, 0, 10)))
#' seg_limits(term_build(seg(x, npsi = 2), d))
#'
#' @export
seg_limits <- function(term) {
  if (!S7::S7_inherits(term, SegTerm)) {
    stop("'term' must be a break-point term.", call. = FALSE)
  }
  .assert_built(term)
  bp <- term@build_info
  K <- bp$npsi
  out <- if (is.null(bp$lims)) matrix(rep(bp$lim, each = K), K, 2L) else bp$lims
  dimnames(out) <- list(paste0("psi", seq_len(K)), c("lower", "upper"))
  out
}
