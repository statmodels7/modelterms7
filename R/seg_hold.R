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
#' family accepts the positions only where its own objective improves.
#'
#' A held term whose break-point is developed, \eqn{\psi_i = w_i'p} with
#' \eqn{w_i} the row of the sub-design, is polished in the coefficients
#' \eqn{p}. Along a line \eqn{p + tv} the position of observation \eqn{i}
#' crosses \eqn{x_i} at \eqn{t_i = (x_i - \psi_i)/(w_i'v)}, so the profile
#' is evaluated once between each pair of consecutive crossings, and for a
#' [jseg()] whose change of slope is not developed like its position it is
#' then minimized inside the best interval, where it varies smoothly. Where
#' the sub-design partitions the observations into groups (`psi ~ g`,
#' `by = ~ 0 + g`) the lines move one group's position at a time, which is a
#' search over every interval of that group, and the sweep also starts from
#' each group's own minimum, the minimum of the profile on that group's
#' observations alone; the better of the two results is kept. With a
#' continuous covariate in the sub-design the lines are the coordinate
#' directions and eight fixed directions in each plane of two coordinates,
#' so the result is the best point found along those lines and not a global
#' minimum. A developed term is polished only when held (see [seg_hold()]).
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
  if (S7::S7_inherits(term, SegTerm) && isTRUE(term@blueprint$developed) &&
      !identical(term@blueprint$kind, "seg") &&
      is.null(term@blueprint$smooth)) {
    return(.seg_polish_developed(term, y, weights, sweeps))
  }
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

# The polish of a held term whose break-points carry a development: exact
# line searches in the coefficients p_k of each position, psi_i = w_i'p_k.
# See seg_polish_exact().
.seg_polish_developed <- function(term, y, weights = NULL, sweeps = 10) {
  bp <- term@blueprint
  if (!isTRUE(bp$held)) {
    stop("a break-point term with a development is polished held; see seg_hold().",
         call. = FALSE)
  }
  xv <- bp$xv
  n <- length(xv)
  y <- as.numeric(y)
  if (length(y) != n) {
    stop("'y' must have one value per observation of the build data.",
         call. = FALSE)
  }
  w <- if (is.null(weights)) rep(1, n) else as.numeric(weights)
  if (length(w) != n || anyNA(w) || any(w < 0)) {
    stop("'weights' must be non-negative, one per observation.", call. = FALSE)
  }
  sw <- sqrt(w)
  ys <- sw * y
  K <- bp$npsi
  jseg <- identical(bp$kind, "jseg")
  lim <- bp$lim

  # the profile: least squares of y on an intercept and the term's columns
  rss_of <- function(tm) {
    X <- as.matrix(tm@X)
    keep <- colSums(abs(X)) > 0
    Z <- cbind(1, X[, keep, drop = FALSE]) * sw
    out <- tryCatch(sum(qr.resid(qr(Z), ys)^2), error = function(e) Inf)
    if (is.finite(out)) out else Inf
  }
  at_pk <- function(tm, k, p) {
    b <- tm@blueprint
    b$pk[[k]] <- p
    asm <- .seg_assemble(b, b$xv, b$coef)
    X <- asm$X
    colnames(X) <- tm@coef_names
    b$value <- asm$value
    b$psi <- asm$psi
    b$pk <- asm$pk
    tm@X <- X
    tm@blueprint <- b
    tm
  }
  design_of <- function(k) {
    z <- bp$Z[[paste0("psi", k)]]
    if (is.null(z)) matrix(1, n, 1L) else as.matrix(z)
  }

  # one exact line search for break-point k along v, the positions moving
  # only inside `range`; NULL where no candidate improves
  line <- function(tm, k, v, range, cur) {
    Zk <- design_of(k)
    p0 <- as.numeric(tm@blueprint$pk[[k]])
    s <- as.numeric(Zk %*% v)
    psi0 <- as.numeric(Zk %*% p0)
    mv <- abs(s) > 1e-12 * max(abs(s))
    if (!any(mv)) return(NULL)
    tb <- sort(unique((xv[mv] - psi0[mv]) / s[mv]))
    if (length(tb) < 2L) return(NULL)
    cand <- (tb[-1L] + tb[-length(tb)]) / 2
    inside <- vapply(cand, function(t) {
      q <- psi0[mv] + t * s[mv]
      all(q > range[1L] & q < range[2L])
    }, logical(1))
    if (!any(inside)) return(NULL)
    ix <- which(inside)
    r <- vapply(cand[ix], function(t) rss_of(at_pk(tm, k, p0 + t * v)),
                numeric(1))
    i <- which.min(r)
    if (!length(i) || !is.finite(r[i])) return(NULL)
    tbest <- cand[ix[i]]
    rbest <- r[i]
    # a jseg whose change of slope is not developed like its position varies
    # inside the interval: (x - psi_i)_+ moves with t on the active rows
    if (jseg) {
      j <- ix[i]
      o <- tryCatch(stats::optimize(
        function(t) rss_of(at_pk(tm, k, p0 + t * v)),
        c(tb[j], tb[j + 1L])), error = function(e) NULL)
      if (!is.null(o) && is.finite(o$objective) && o$objective < rbest) {
        tbest <- o$minimum
        rbest <- o$objective
      }
    }
    if (rbest < cur - 1e-10 * (cur + 1)) {
      list(tm = at_pk(tm, k, p0 + tbest * v), rss = rbest)
    } else NULL
  }

  # the groups a sub-design defines, where its distinct rows are as many as
  # its columns and of full rank: then each group's position is free
  groups_of <- function(Zk) {
    key <- apply(Zk, 1L, function(r) paste(format(r, digits = 17), collapse = "\r"))
    first <- !duplicated(key)
    U <- Zk[first, , drop = FALSE]
    if (nrow(U) != ncol(Zk) || qr(U)$rank < ncol(Zk)) return(NULL)
    list(U = U, grp = match(key, key[first]))
  }
  # the directions for break-point k and the range each may move in
  plan_of <- function(k) {
    Zk <- design_of(k)
    gr <- groups_of(Zk)
    if (!is.null(gr)) {
      dirs <- lapply(seq_len(ncol(Zk)), function(j) {
        e <- numeric(ncol(Zk))
        e[j] <- 1
        v <- solve(gr$U, e)
        # U v = e_j: only the rows of group j move
        rows <- which(gr$grp == j)
        xg <- xv[rows]
        qg <- as.numeric(stats::quantile(xg, c(0.05, 0.95), names = FALSE))
        list(v = v, range = c(max(lim[1L], qg[1L]), min(lim[2L], qg[2L])),
             rows = rows)
      })
      return(list(groups = gr, dirs = dirs))
    }
    # with a continuous covariate the positions are confined by the clamp of
    # .seg_positions() alone, so every crossing is a candidate
    q <- ncol(Zk)
    a <- sqrt(colMeans(Zk^2))
    a[a == 0] <- 1
    dirs <- list()
    free <- c(-Inf, Inf)
    for (j in seq_len(q)) {
      e <- numeric(q)
      e[j] <- 1 / a[j]
      dirs[[length(dirs) + 1L]] <- list(v = e, range = free)
    }
    if (q > 1L) {
      for (i in seq_len(q - 1L)) for (j in (i + 1L):q) {
        for (th in (1:7) * pi / 8) {
          if (abs(th - pi / 2) < 1e-12) next
          e <- numeric(q)
          e[i] <- cos(th) / a[i]
          e[j] <- sin(th) / a[j]
          dirs[[length(dirs) + 1L]] <- list(v = e, range = free)
        }
      }
    }
    list(groups = NULL, dirs = dirs)
  }

  plans <- lapply(seq_len(K), plan_of)
  descend <- function(tm) {
    cur <- rss_of(tm)
    for (s in seq_len(as.integer(sweeps))) {
      moved <- FALSE
      for (k in seq_len(K)) {
        for (d in plans[[k]]$dirs) {
          res <- line(tm, k, d$v, d$range, cur)
          if (!is.null(res)) {
            tm <- res$tm
            cur <- res$rss
            moved <- TRUE
          }
        }
      }
      if (!moved) break
    }
    list(tm = tm, rss = cur)
  }

  # each group's own minimum: the profile on that group's observations alone,
  # the other break-points held where they are
  own_start <- function(tm) {
    for (k in seq_len(K)) {
      pl <- plans[[k]]
      if (is.null(pl$groups)) next
      psi <- tm@blueprint$psi
      qk <- numeric(length(pl$dirs))
      for (j in seq_along(pl$dirs)) {
        d <- pl$dirs[[j]]
        rows <- d$rows
        xg <- xv[rows]
        u <- sort(unique(xg))
        mid <- (u[-1L] + u[-length(u)]) / 2
        mid <- mid[mid > d$range[1L] & mid < d$range[2L]]
        cur_k <- psi[rows[1L], k]
        if (!length(mid)) {
          qk[j] <- cur_k
          next
        }
        fixed <- matrix(1, length(xg), 1L)
        if (bp$linear) fixed <- cbind(fixed, xg)
        for (kk in setdiff(seq_len(K), k)) {
          pk2 <- psi[rows, kk]
          fixed <- cbind(fixed, if (jseg) pmax(xg - pk2, 0), as.numeric(xg > pk2))
        }
        swg <- sw[rows]
        yg <- ys[rows]
        r <- vapply(mid, function(q) {
          Zl <- cbind(fixed, if (jseg) pmax(xg - q, 0), as.numeric(xg > q)) * swg
          out <- tryCatch(sum(qr.resid(qr(Zl), yg)^2), error = function(e) Inf)
          if (is.finite(out)) out else Inf
        }, numeric(1))
        qk[j] <- if (any(is.finite(r))) mid[which.min(r)] else cur_k
      }
      tm <- at_pk(tm, k, as.numeric(solve(pl$groups$U, qk)))
    }
    tm
  }

  best <- descend(term)
  if (any(vapply(plans, function(pl) !is.null(pl$groups), logical(1)))) {
    alt <- descend(own_start(term))
    if (alt$rss < best$rss) best <- alt
  }
  if (best$rss < rss_of(term) - 1e-10 * (rss_of(term) + 1)) best$tm else term
}
