#' @include gas_joint.R
NULL

# ---------------------------------------------------------------- jets
#
# A jet is a quantity with its derivatives in the caller's unknowns: `v`, the
# gradient `d1`, the Hessian `d2`, the third derivative contracted against
# each direction `d3` (a list), and the fourth contracted against both
# directions `d4`. The number of directions sets how far a jet goes.

.gas_jzero <- function(v, m, nd) {
  z <- matrix(0, m, m)
  list(v = v, d1 = numeric(m), d2 = z,
       d3 = if (nd >= 1L) rep(list(z), nd) else NULL,
       d4 = if (nd >= 2L) z else NULL)
}

.gas_jsum <- function(A, B) {
  A$v <- A$v + B$v
  A$d1 <- A$d1 + B$d1
  A$d2 <- A$d2 + B$d2
  if (!is.null(A$d3)) A$d3 <- Map(`+`, A$d3, B$d3)
  if (!is.null(A$d4)) A$d4 <- A$d4 + B$d4
  A
}

# the product of two jets: Leibniz's rule at every order, the fourth through
# .gas_prod4(), which already holds the sixteen terms
.gas_jprod <- function(A, B, dirs) {
  nd <- length(dirs)
  out <- list(v = A$v * B$v,
              d1 = A$v * B$d1 + B$v * A$d1,
              d2 = A$v * B$d2 + B$v * A$d2 + outer(A$d1, B$d1) +
                outer(B$d1, A$d1),
              d3 = NULL, d4 = NULL)
  if (nd >= 1L) {
    out$d3 <- lapply(seq_len(nd), function(d) {
      x <- dirs[[d]]
      A2x <- as.numeric(A$d2 %*% x)
      B2x <- as.numeric(B$d2 %*% x)
      A$v * B$d3[[d]] + B$v * A$d3[[d]] + sum(A$d1 * x) * B$d2 +
        sum(B$d1 * x) * A$d2 + outer(A2x, B$d1) + outer(B$d1, A2x) +
        outer(B2x, A$d1) + outer(A$d1, B2x)
    })
  }
  if (nd >= 2L) {
    out$d4 <- .gas_prod4(
      .gas_jet(A$v, A$d1, A$d2, A$d3[[1L]], A$d3[[2L]], A$d4, dirs),
      .gas_jet(B$v, B$d1, B$d2, B$d3[[1L]], B$d3[[2L]], B$d4, dirs))
  }
  out
}

# the jet of a coordinate that rides the identity link, as a cross loading
# does: its value, and a one at its own column
.gas_jcoord <- function(v, col, m, nd) {
  J <- .gas_jzero(v, m, nd)
  J$d1[col] <- 1
  J
}

# the jets of a term's chart, from .gas_chart_prep()
.gas_jchart <- function(gp, nd) {
  mk <- function(v, d1, d2, d3, d4) {
    list(v = v, d1 = d1, d2 = d2, d3 = if (nd >= 1L) d3 else NULL,
         d4 = if (nd >= 2L) d4 else NULL)
  }
  p <- length(gp$a)
  q <- length(gp$b)
  list(
    omega = mk(gp$omega, gp$om_u, gp$om_uu, gp$t_om, gp$q_om),
    a = lapply(seq_len(p), function(i) mk(
      gp$a[[i]], gp$a_u[[i]], gp$a_uu[[i]],
      if (nd >= 1L) lapply(seq_len(nd), function(d) gp$t_a[[d]][[i]]),
      if (nd >= 2L) gp$q_a[[i]])),
    b = lapply(seq_len(q), function(j) mk(
      gp$b[[j]], gp$b_u[[j]], gp$b_uu[[j]],
      if (nd >= 1L) lapply(seq_len(nd), function(d) gp$t_b[[d]][[j]]),
      if (nd >= 2L) gp$q_b[[j]])),
    f0 = mk(gp$f0, gp$f0_u, gp$f0_uu, gp$f0_3, gp$f0_4))
}

# The jet of a driving score from the jets of the predictors of every
# equation and the driving score's derivatives in those predictors at the
# row: U1 a vector, U2 a matrix, U3 and U4 arrays over the predictors. This
# is the chain rule for a function of several predictors, order by order;
# with Vm the gradients side by side,
#   d1 = Vm U1,
#   d2 = sum_a U1_a Phi_a + Vm U2 Vm',
# and the third and fourth orders differentiate these along the directions.
.gas_jdrive <- function(v, E, U1, U2, U3, U4, dirs) {
  nd <- length(dirs)
  P <- length(E)
  m <- length(E[[1L]]$d1)
  Vm <- matrix(vapply(E, function(e) e$d1, numeric(m)), m, P)
  out <- list(v = v, d1 = as.numeric(Vm %*% U1), d3 = NULL, d4 = NULL)
  d2 <- Vm %*% U2 %*% t(Vm)
  for (a in seq_len(P)) d2 <- d2 + U1[a] * E[[a]]$d2
  out$d2 <- d2
  if (nd == 0L) return(out)
  # the contractions along each direction
  dx <- lapply(dirs, function(x) as.numeric(crossprod(Vm, x)))
  Px <- lapply(dirs, function(x)
    matrix(vapply(E, function(e) as.numeric(e$d2 %*% x), numeric(m)), m, P))
  U2x <- lapply(dx, function(z) as.numeric(U2 %*% z))
  U3x <- lapply(dx, function(z) {
    M <- matrix(0, P, P)
    for (c in seq_len(P)) M <- M + U3[, , c] * z[c]
    M
  })
  out$d3 <- lapply(seq_len(nd), function(d) {
    r <- Vm %*% U3x[[d]] %*% t(Vm) + Px[[d]] %*% U2 %*% t(Vm) +
      Vm %*% U2 %*% t(Px[[d]])
    for (a in seq_len(P)) {
      r <- r + U1[a] * E[[a]]$d3[[d]] + U2x[[d]][a] * E[[a]]$d2
    }
    r
  })
  if (nd == 1L) return(out)
  v1 <- dirs[[1L]]
  w1 <- dirs[[2L]]
  dv <- dx[[1L]]
  dw <- dx[[2L]]
  Pv <- Px[[1L]]
  Pw <- Px[[2L]]
  # v' Phi_c w, and each predictor's third derivative along w applied to v
  phivw <- as.numeric(crossprod(Pv, w1))
  Psiwv <- matrix(vapply(E, function(e) as.numeric(e$d3[[2L]] %*% v1),
                         numeric(m)), m, P)
  U4vw <- matrix(0, P, P)
  U3phi <- matrix(0, P, P)
  for (c in seq_len(P)) {
    U3phi <- U3phi + U3[, , c] * phivw[c]
    for (e in seq_len(P)) U4vw <- U4vw + U4[, , c, e] * (dv[c] * dw[e])
  }
  U3wv <- as.numeric(U3x[[2L]] %*% dv)
  U2phi <- as.numeric(U2 %*% phivw)
  r <- Vm %*% (U4vw + U3phi) %*% t(Vm) +
    Pw %*% U3x[[1L]] %*% t(Vm) + Vm %*% U3x[[1L]] %*% t(Pw) +
    Pv %*% U3x[[2L]] %*% t(Vm) + Vm %*% U3x[[2L]] %*% t(Pv) +
    Psiwv %*% U2 %*% t(Vm) + Vm %*% U2 %*% t(Psiwv) +
    Pv %*% U2 %*% t(Pw) + Pw %*% U2 %*% t(Pv)
  for (a in seq_len(P)) {
    r <- r + (U3wv[a] + U2phi[a]) * E[[a]]$d2 +
      U2x[[1L]][a] * E[[a]]$d3[[2L]] + U2x[[2L]][a] * E[[a]]$d3[[1L]] +
      U1[a] * E[[a]]$d4
  }
  out$d4 <- r
  out
}

# The second order of the joint filter on the compiled route: the term charts
# lifted to the unknowns by .gas_chart_prep() (no direction), flattened to
# plain vectors, and gas_joint_curv_cpp() running the recursion. The result has
# the structure that the R recursion returns.
.gas_joint_core_cpp <- function(term, psi, cf, g, seed, fp, u_at, U1, U2, m,
                                zoff, threads) {
  K <- length(term@params)
  P <- length(seed)
  cols <- function(lst) {
    if (length(lst)) as.numeric(do.call(cbind, lst)) else numeric(0)
  }
  cube <- function(lst) as.numeric(unlist(lst))
  charts <- vector("list", K)
  for (k in seq_len(K)) {
    tm <- term@terms[[term@params[k]]]
    base <- .gas_base_params(tm@p, tm@q)
    links <- term_links(tm)
    v <- unlist(psi[[term@params[k]]][base])
    zeta <- vapply(base, function(j)
      linkfunctions7::linkfun(links[[j]], v[[j]]), numeric(1))
    gp <- .gas_chart_prep(zeta, tm@p, tm@q, links,
                          zoff + term@offsets[k] + seq_along(base),
                          seq_len(m), list())
    charts[[k]] <- list(
      omega = gp$omega, om_u = gp$om_u, om_uu = as.numeric(gp$om_uu),
      f0 = gp$f0, f0_u = gp$f0_u, f0_uu = as.numeric(gp$f0_uu),
      a = as.numeric(unlist(gp$a)), a_u = cols(gp$a_u), a_uu = cube(gp$a_uu),
      b = as.numeric(unlist(gp$b)), b_u = cols(gp$b_u), b_uu = cube(gp$b_uu))
  }
  out <- gas_joint_curv_cpp(seed, charts, cf, as.integer(fp), g,
                            matrix(as.numeric(u_at), nrow(g), K),
                            as.numeric(U1), as.numeric(U2),
                            term@blueprint$order, as.integer(m),
                            as.integer(zoff), as.integer(threads))
  names(out$jacobian) <- term@params
  out
}

# The body that term_curvature(), term_third() and term_fourth() share for a
# joint filter: the levels' jets propagated through the recursion, the
# driving scores' jets formed from the predictors' at every row.
.gas_joint_core <- function(term, psi, g, seed, blocks, direction = NULL,
                            threads = getOption("modelterms7.joint_threads",
                                                1L)) {
  bp <- term@blueprint
  K <- length(term@params)
  n <- bp$n
  dirs <- .gas_dirs(direction)
  nd <- length(dirs)
  if (!is.list(seed) || !length(seed)) {
    stop("'seed' must be a list with one matrix of static rows per equation.",
         call. = FALSE)
  }
  seed <- lapply(seed, as.matrix)
  P <- length(seed)
  m <- ncol(seed[[1L]])
  M <- length(term_params(term))
  if (m < M) {
    stop(sprintf("'seed' has %d columns and the filter has %d parameters.",
                 m, M), call. = FALSE)
  }
  fp <- blocks$filtered
  if (length(fp) != K || anyNA(fp) || any(fp < 1L | fp > P)) {
    stop("'blocks$filtered' must give the equation of each filtered parameter.",
         call. = FALSE)
  }
  g <- matrix(as.numeric(g), n, K)
  dirs <- lapply(dirs, function(d) {
    d <- as.numeric(d)
    if (length(d) != m) {
      stop(sprintf("'direction' must have one value per unknown (%d).", m),
           call. = FALSE)
    }
    d
  })
  u_at <- matrix(blocks$u, n, K)
  U1 <- blocks$U1
  U2 <- blocks$U2
  U3 <- if (nd >= 1L) blocks$U3 else NULL
  U4 <- if (nd >= 2L) blocks$U4 else NULL
  if ((nd >= 1L && is.null(U3)) || (nd >= 2L && is.null(U4))) {
    stop("'blocks' must carry U3 for the third order and U4 for the fourth.",
         call. = FALSE)
  }

  # the chart jets of each term, lifted to the unknowns: the filter's own
  # parameters are the LAST M columns, in the order of term_params()
  zoff <- m - M
  psi <- .gas_joint_psi(term, psi)
  cf <- .gas_joint_coefs(term, psi)

  # the second order (no direction) on the compiled route; the R recursion
  # below stays reachable with options(modelterms7.joint_core_r = TRUE)
  if (nd == 0L && !isTRUE(getOption("modelterms7.joint_core_r", FALSE))) {
    return(.gas_joint_core_cpp(term, psi, cf, g, seed, fp, u_at, U1, U2, m,
                               zoff, as.integer(threads)))
  }
  charts <- vector("list", K)
  for (k in seq_len(K)) {
    tm <- term@terms[[term@params[k]]]
    base <- .gas_base_params(tm@p, tm@q)
    links <- term_links(tm)
    v <- unlist(psi[[term@params[k]]][base])
    zeta <- vapply(base, function(j)
      linkfunctions7::linkfun(links[[j]], v[[j]]), numeric(1))
    gp <- .gas_chart_prep(zeta, tm@p, tm@q, links,
                          zoff + term@offsets[k] + seq_along(base),
                          seq_len(m), dirs)
    charts[[k]] <- .gas_jchart(gp, nd)
  }

  jac <- lapply(seq_len(K), function(k) matrix(0, n, m))
  dphi <- lapply(seq_len(K), function(k)
    lapply(seq_len(nd), function(d) matrix(0, n, m)))
  dpsi <- if (nd >= 2L) lapply(seq_len(K), function(k) matrix(0, n, m))
  W <- matrix(0, m, m)
  zero <- .gas_jzero(0, m, nd)
  for (rows in bp$order) {
    T_ <- length(rows)
    fj <- lapply(seq_len(K), function(k) vector("list", T_))
    uj <- lapply(seq_len(K), function(k) vector("list", T_))
    for (t in seq_len(T_)) {
      row <- rows[t]
      for (k in seq_len(K)) {
        ch <- charts[[k]]
        c_k <- cf[[k]]
        J <- ch$omega
        for (i in seq_len(c_k$p)) {
          if (t - i < 1L) next
          J <- .gas_jsum(J, .gas_jprod(ch$a[[i]], uj[[k]][[t - i]], dirs))
          for (j in seq_along(c_k$partners)) {
            l <- c_k$partners[j]
            cj <- .gas_jcoord(c_k$cross[j, i], zoff + c_k$col_cross[j, i],
                              m, nd)
            J <- .gas_jsum(J, .gas_jprod(cj, uj[[l]][[t - i]], dirs))
          }
        }
        for (j in seq_len(c_k$q)) {
          lagj <- if (t - j >= 1L) fj[[k]][[t - j]] else ch$f0
          J <- .gas_jsum(J, .gas_jprod(ch$b[[j]], lagj, dirs))
        }
        fj[[k]][[t]] <- J
      }
      # the predictors of every equation at this row: the static rows, and
      # for a filtered equation the level's jet on top
      E <- lapply(seq_len(P), function(a) {
        e <- zero
        e$d1 <- as.numeric(seed[[a]][row, ])
        e
      })
      for (k in seq_len(K)) {
        a <- fp[k]
        J <- fj[[k]][[t]]
        E[[a]]$d1 <- E[[a]]$d1 + J$d1
        E[[a]]$d2 <- J$d2
        if (nd >= 1L) E[[a]]$d3 <- J$d3
        if (nd >= 2L) E[[a]]$d4 <- J$d4
        jac[[k]][row, ] <- E[[a]]$d1
        for (d in seq_len(nd)) {
          dphi[[k]][[d]][row, ] <- as.numeric(J$d2 %*% dirs[[d]])
        }
        if (nd >= 2L) {
          dpsi[[k]][row, ] <- as.numeric(J$d3[[1L]] %*% dirs[[2L]])
        }
        top <- if (nd >= 2L) J$d4 else if (nd == 1L) J$d3[[1L]] else J$d2
        W <- W + g[row, k] * top
      }
      for (l in seq_len(K)) {
        uj[[l]][[t]] <- .gas_jdrive(
          u_at[row, l], E, U1[row, l, ], matrix(U2[row, l, , ], P, P),
          if (nd >= 1L) array(U3[row, l, , , ], c(P, P, P)),
          if (nd >= 2L) array(U4[row, l, , , , ], c(P, P, P, P)), dirs)
      }
    }
  }
  # the entries and their transposes collect the same terms in a different
  # order, so the sum is symmetric to a rounding and is made so exactly
  W <- (W + t(W)) / 2
  names(jac) <- names(dphi) <- term@params
  if (nd >= 2L) names(dpsi) <- term@params
  out <- list(jacobian = jac, curvature = W)
  if (nd >= 1L) out$dphi <- dphi
  if (nd >= 2L) out$dpsi <- dpsi
  out
}

#' @title The Curvature of a Joint Score-Driven Filter
#' @name term_curvature.GasJoint
#' @description
#' Computes the Jacobian of every filtered predictor in a caller's unknowns
#' and the second derivative of the filtered predictors contracted against
#' the caller's weights, propagated through the joint recursion.
#' @details
#' The level of each parameter and each driving score are carried as jets: a
#' value with its derivatives. The jet of a level is a sum of products of the
#' chart's jets with the jets of earlier driving scores and levels, and the
#' jet of a driving score at a row is formed from the jets of the predictors
#' of every equation by the chain rule, with the driving score's derivatives
#' in those predictors supplied as data in `blocks`. The static rows of every
#' equation enter through `seed`, and only a filtered equation's predictor has
#' a second derivative.
#' @param term A [GasJoint()].
#' @param eta,y,score,curvature Unused: the recursion reads the data in
#'   `blocks`. They are kept for the signature of the generic.
#' @param psi The terms' parameters, as for [term_filter()].
#' @param g The weights, a matrix with one column per filtered parameter.
#' @param seed A list with one matrix per distribution parameter: the
#'   derivative of that equation's static predictor in the unknowns, one row
#'   per observation. The joint filter's parameters are the last columns, in
#'   the order of [term_params()].
#' @param blocks A list with `filtered` (the position in `seed` of the
#'   equation of each filtered parameter), `u` (the driving scores at the
#'   fitted predictors, a matrix with one column per filtered parameter),
#'   `U1` and `U2` (the first and second derivatives of each driving score in
#'   every predictor, arrays of dimension `c(n, K, P)` and `c(n, K, P, P)`),
#'   and, for the third and fourth orders, `U3` and `U4`, with one more
#'   dimension of length `P` each.
#' @param ... Unused.
#' @return A list with `jacobian` (one matrix per filtered parameter) and
#'   `curvature`.
#' @keywords internal
S7::method(term_curvature, GasJoint) <- function(term, eta, y, score,
                                                 curvature, psi, g, seed,
                                                 blocks, ...) {
  .gas_joint_core(term, psi, g, seed, blocks)
}

#' @title Third Derivatives of a Joint Score-Driven Filter
#' @name term_third.GasJoint
#' @description
#' Differentiates the second derivative of the filtered predictors once more
#' along one direction, through the jets of [term_curvature()].
#' @param term A [GasJoint()].
#' @param eta,y,score,curvature Unused, as for [term_curvature()].
#' @param psi,g,seed,blocks As for [term_curvature()]; `blocks` carries `U3`.
#' @param direction The direction against which to contract.
#' @param ... Unused.
#' @return A list with `jacobian`, `dphi` (per filtered parameter, a list
#'   with one matrix per direction holding the second derivative of the
#'   predictor applied to the direction) and `curvature`.
#' @keywords internal
S7::method(term_third, GasJoint) <- function(term, eta, y, score, curvature,
                                             psi, g, seed, blocks, direction,
                                             ...) {
  .gas_joint_core(term, psi, g, seed, blocks, direction = direction)
}

#' @title The Fourth Derivative of a Joint Score-Driven Filter
#' @name term_fourth.GasJoint
#' @description
#' Computes the fourth derivative of the filtered predictors contracted
#' against the caller's weights and two directions, through the jets of
#' [term_curvature()].
#' @param term A [GasJoint()].
#' @param eta,y,score,curvature Unused, as for [term_curvature()].
#' @param psi,g,seed,blocks As for [term_curvature()]; `blocks` carries `U3`
#'   and `U4`.
#' @param directions A list of two directions.
#' @param ... Unused.
#' @return A list with `jacobian`, `dphi`, `dpsi` (per filtered parameter,
#'   the third derivative along the first direction applied to the second)
#'   and `curvature`.
#' @keywords internal
S7::method(term_fourth, GasJoint) <- function(term, eta, y, score, curvature,
                                              psi, g, seed, blocks,
                                              directions, ...) {
  if (!is.list(directions) || length(directions) != 2L) {
    stop("'directions' must be a list of two vectors, one per direction.",
         call. = FALSE)
  }
  .gas_joint_core(term, psi, g, seed, blocks, direction = directions)
}
