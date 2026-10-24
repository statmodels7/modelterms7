#' @include gas_joint_core.R
NULL

#' @name term_continue.GasJoint
#'
#' @title Continuing a Joint Score-Driven Filter Past the Series
#'
#' @description
#' Computes the levels that the joint recursion reaches at rows that come
#' after the observed ones.
#'
#' @details
#' Past the data each driving score is at its conditional mean, zero, so the
#' continuation is the deterministic recursion of each level,
#' \deqn{f_{k,n+h} = \omega_k + \sum_i \Big(\kappa_{k,i}\, u_{k,n+h-i}
#'   + \sum_{l \in C_k} \kappa_{k,i,l}\, u_{l,n+h-i}\Big)
#'   + \sum_j \xi_{k,j} f_{k,n+h-j},}
#' in which the loadings contribute only while \eqn{n+h-i} is an observed
#' time. A cross loading therefore reads the last observed scores of the
#' other parameters of its group in the first steps, after which the levels
#' evolve each on its own. The rules on the new rows are those of
#' [term_continue()] for one term: each group continues its own series, the
#' rows of a group are taken in time order and must come after every observed
#' time of that group, and the terms must have been built with `time`.
#'
#' @param term A built [GasJoint()].
#' @param psi The parameters, as for [term_filter()].
#' @param f_past The levels at the observed rows, a matrix with one column per
#'   filtered parameter.
#' @param s_past The driving scores at the observed rows, a matrix like
#'   `f_past`.
#' @param newdata The rows onto which to continue.
#' @param ... Ignored.
#' @param deriv `NULL`, or a list with `df_past` and `ds_past`, lists with one
#'   matrix per filtered parameter (one row per observed row, `m` columns),
#'   and `dpsi`, the derivative of the joint parameter vector, one row per
#'   entry of [term_params()] and `m` columns.
#'
#' @return A matrix of continued levels with one row per row of `newdata` and
#'   one column per filtered parameter. With `deriv`, it carries the attribute
#'   `"jacobian"`, a list with one `nrow(newdata)` by `m` matrix per filtered
#'   parameter.
#'
#' @seealso [term_continue()], [term_filter()]
#'
#' @keywords internal
S7::method(term_continue, GasJoint) <- function(term, psi, f_past, s_past,
                                                newdata, ..., deriv = NULL) {
  K <- length(term@params)
  nn <- nrow(newdata)
  out <- matrix(0, nn, K, dimnames = list(NULL, term@params))
  if (!nn) return(out)
  tm1 <- term@terms[[1L]]
  bp1 <- tm1@blueprint
  if (is.null(tm1@time)) {
    stop("continuing a score-driven term past the series needs its `time`: ",
         "without one, a row's place in the series is its position in the ",
         "data frame, which gives no information about new rows.",
         call. = FALSE)
  }
  gf <- if (is.null(tm1@by)) factor(rep(1L, nn), levels = "1") else
    factor(as.character(eval(tm1@by, newdata, baseenv())),
           levels = bp1$levels)
  tv <- eval(tm1@time, newdata, baseenv())
  if (anyNA(gf)) {
    stop(sprintf(paste0("rows %s name a group that was not in the data of ",
                        "the fit, so there is no state to continue."),
                 paste(utils::head(which(is.na(gf)), 5L), collapse = ", ")),
         call. = FALSE)
  }
  if (anyNA(tv)) stop("`time` must not be missing.", call. = FALSE)
  f_past <- matrix(as.numeric(f_past), ncol = K)
  s_past <- matrix(as.numeric(s_past), ncol = K)
  # the coefficients at every new row with their derivative rows in the
  # joint vector; constant for an undeveloped term, read from the new rows
  # through the subformulas for a developed one
  cf <- .gas_joint_new_coefs(term, psi, newdata)
  M <- length(term_params(term))
  m <- if (!is.null(deriv)) ncol(deriv$dpsi) else 0L
  lift <- function(g) as.numeric(g %*% deriv$dpsi)
  if (!is.null(deriv)) jac <- lapply(seq_len(K), function(k) matrix(0, nn, m))

  gi <- as.integer(gf)
  for (lv in sort(unique(gi))) {
    rows_old <- bp1$order[[as.character(lv)]]
    if (is.null(rows_old)) {
      stop(paste0("a group of the new data has no observed rows from which",
                  " to continue."), call. = FALSE)
    }
    new <- which(gi == lv)
    new <- new[order(tv[new])]
    last <- max(bp1$times[rows_old])
    if (any(tv[new] <= last)) {
      bad <- new[tv[new] <= last]
      stop(sprintf(paste0("rows %s fall inside the observed series (their ",
                          "group ends at time %s).\n  A row whose response ",
                          "is known is evaluated by running the filter on ",
                          "it."),
                   paste(utils::head(bad, 5L), collapse = ", "),
                   format(last)), call. = FALSE)
    }
    n0 <- length(rows_old)
    nh <- length(new)
    f <- rbind(f_past[rows_old, , drop = FALSE], matrix(0, nh, K))
    sc <- rbind(s_past[rows_old, , drop = FALSE], matrix(0, nh, K))
    if (!is.null(deriv)) {
      dF <- lapply(seq_len(K), function(k)
        rbind(deriv$df_past[[k]][rows_old, , drop = FALSE], matrix(0, nh, m)))
      dS <- lapply(seq_len(K), function(k)
        rbind(deriv$ds_past[[k]][rows_old, , drop = FALSE], matrix(0, nh, m)))
    }
    for (h in seq_len(nh)) {
      t <- n0 + h
      r <- new[[h]]
      for (k in seq_len(K)) {
        c_k <- cf[[k]]
        dcoef <- !is.null(deriv)
        ft <- c_k$om[r]
        dft <- if (dcoef) lift(c_k$dom[r, ])
        for (i in seq_len(c_k$p)) {
          if (t - i < 1L || t - i > n0) next
          ft <- ft + c_k$A[r, i] * sc[t - i, k]
          if (dcoef) {
            dft <- dft + sc[t - i, k] * lift(c_k$dA[[i]][r, ]) +
              c_k$A[r, i] * dS[[k]][t - i, ]
          }
        }
        for (x in seq_along(c_k$xlag)) {
          i <- c_k$xlag[x]
          if (t - i < 1L || t - i > n0) next
          l <- c_k$xpar[x]
          ft <- ft + c_k$cross[[x]][r] * sc[t - i, l]
          if (dcoef) {
            dft <- dft + sc[t - i, l] * lift(c_k$dcross[[x]][r, ]) +
              c_k$cross[[x]][r] * dS[[l]][t - i, ]
          }
        }
        for (j in seq_len(c_k$q)) {
          # a lag before the series is the starting level, which the first
          # level equals
          lag <- max(t - j, 1L)
          ft <- ft + c_k$B[r, j] * f[lag, k]
          if (dcoef) {
            dft <- dft + f[lag, k] * lift(c_k$dB[[j]][r, ]) +
              c_k$B[r, j] * dF[[k]][lag, ]
          }
        }
        f[t, k] <- ft
        out[r, k] <- ft
        if (dcoef) {
          dF[[k]][t, ] <- dft
          jac[[k]][r, ] <- dft
        }
      }
    }
  }
  if (!is.null(deriv)) {
    names(jac) <- term@params
    attr(out, "jacobian") <- jac
  }
  out
}

# The coefficients of the joint recursion at the rows of `newdata`: per
# term, the level, the loadings and the autoregressive coefficients at each
# new row (constant for a term without a subformula), the cross loadings,
# and the derivative row of each in the joint parameter vector, on the scale
# of `psi` that the filter reads.
.gas_joint_new_coefs <- function(term, psi, newdata) {
  psi <- .gas_joint_psi(term, psi)
  nn <- nrow(newdata)
  M <- length(term_params(term))
  out <- vector("list", length(term@params))
  for (k in seq_along(term@params)) {
    tm <- term@terms[[term@params[k]]]
    nm <- term_params(tm)
    v <- unlist(psi[[term@params[k]]][nm])
    p <- tm@p
    q <- tm@q
    off <- term@offsets[k]
    xn <- .gas_cross_params(tm)
    nb <- length(nm) - length(xn)
    lifted <- function(D) {
      G <- matrix(0, nn, M)
      G[, off + seq_len(ncol(D))] <- D
      G
    }
    if (length(tm@blueprint) && !is.null(tm@blueprint$sub)) {
      nd <- .gas_sub_new(tm, v[seq_len(nb)], newdata)
      om <- nd$om
      A <- nd$A
      B <- nd$B
      dom <- lifted(nd$dom)
      dA <- lapply(nd$dA, lifted)
      dB <- lapply(nd$dB, lifted)
    } else {
      base <- .gas_base_params(p, q)
      cf <- .gas_coefs(v[base], p, q)
      om <- rep(cf$omega, nn)
      A <- matrix(rep(cf$a, each = nn), nn, p)
      B <- matrix(rep(cf$b, each = nn), nn, q)
      unit_rows <- function(col) {
        G <- matrix(0, nn, M)
        G[, off + col] <- 1
        G
      }
      dom <- unit_rows(1L)
      dA <- lapply(seq_len(p), function(i) unit_rows(cf$i_a[i]))
      dB <- lapply(seq_len(q), function(j) {
        G <- matrix(0, nn, M)
        G[, off + seq_len(length(base))] <- matrix(cf$db[j, ], nn,
                                                   length(base), byrow = TRUE)
        G
      })
    }
    # the cross loadings at every new row, partner by partner and lag by
    # lag as .gas_cross_params() orders them, with their derivative rows
    pk <- term@partners[[k]]
    xp <- .gas_cross_pairs(tm)
    xv <- vector("list", length(xn))
    xd <- vector("list", length(xn))
    developed <- length(tm@blueprint) && !is.null(tm@blueprint$sub)
    for (h in seq_along(xn)) {
      if (developed) {
        xv[[h]] <- nd$cross[[xn[h]]]
        xd[[h]] <- lifted(nd$dcross[[xn[h]]])
      } else {
        xv[[h]] <- rep(v[[xn[h]]], nn)
        G <- matrix(0, nn, M)
        G[, off + nb + h] <- 1
        xd[[h]] <- G
      }
    }
    out[[k]] <- list(p = p, q = q, om = om, A = A, B = B, dom = dom,
                     dA = dA, dB = dB, partners = pk,
                     xpar = match(xp$partner, term@params), xlag = xp$lag,
                     cross = xv, dcross = xd)
  }
  out
}

#' @name term_simulate.GasJoint
#'
#' @title Generating From a Joint Score-Driven Filter
#'
#' @description
#' Runs the joint recursion forward, drawing the response at each step.
#'
#' @details
#' As for one term ([term_simulate()]): the filter calls its driving vector
#' once per observation, in time order within each group, at the filtered
#' predictors that it has just produced, so a callback that draws the
#' response there, stores it and returns the driving vector of the drawn
#' value turns the filter into a generator. The curvature is read at the
#' response already drawn.
#'
#' @param term A built [GasJoint()].
#' @param psi The parameters, as for [term_filter()].
#' @param eta The static part of the filtered predictors, a matrix with one
#'   column per filtered parameter.
#' @param draw A function `(e, i)` of the vector of filtered predictors at
#'   row `i`, returning one response value.
#' @param score,curvature Functions `(y, e, i)` returning the driving vector
#'   and its derivative in the filtered predictors, as the callbacks of
#'   [term_filter()] do, at the response `y`.
#' @param ... Ignored.
#'
#' @return A list with `eta` (the filtered predictors), `y` and `latent` (the
#'   levels), matrices with one column per filtered parameter apart from `y`.
#'
#' @seealso [term_simulate()], [term_filter()]
#'
#' @keywords internal
S7::method(term_simulate, GasJoint) <- function(term, psi, eta, draw,
                                                score, curvature, ...) {
  n <- term@blueprint$n
  yv <- rep(NA_real_, n)
  get <- function(e, i) {
    if (is.na(yv[[i]])) yv[[i]] <<- as.numeric(draw(e, i))
    yv[[i]]
  }
  out <- term_filter(term, eta, yv,
                     function(e, i) score(get(e, i), e, i),
                     function(e, i) curvature(get(e, i), e, i), psi)
  eta <- matrix(as.numeric(eta), n, length(term@params))
  list(eta = out$eta, y = yv, latent = out$eta - eta)
}
