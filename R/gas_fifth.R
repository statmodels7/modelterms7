#' @include structural.R gas.R gas_submodels.R gas_joint.R gas_joint_core.R
NULL

# The pieces of the score-driven recursion's fifth order that do not depend on
# how the model's derivatives reach the recursion: the fifth derivative of the
# Levinson-Durbin map, the chart's derivatives of any order up to five in a
# set of directions, and the product rule of any order over a set of
# directions.

#' The Fifth Derivative of the Levinson-Durbin Map, in Three Directions
#'
#' @description
#' Differentiates the fourth derivatives of [gas_levinson4()] once more and
#' contracts them against a third direction, which gives one matrix per
#' coefficient.
#'
#' @details
#' The map is bilinear at every step of the recursion, so differentiating it
#' once more adds no new kind of term. Contracting three of the five slots of
#' the full fifth derivative against \eqn{v}, \eqn{w} and \eqn{x} leaves
#' \deqn{F^{(k)}_i = F^{(k-1)}_i - \rho_k F^{(k-1)}_{k-i}
#'   - v_k Q^{(k-1)}_{k-i}[w,x] - w_k Q^{(k-1)}_{k-i}[v,x]
#'   - x_k Q^{(k-1)}_{k-i}[v,w]
#'   - e_k\big(Q^{(k-1)}_{k-i}[v,w]x\big)^{\top}
#'   - \big(Q^{(k-1)}_{k-i}[v,w]x\big)e_k^{\top},}
#' where \eqn{Q[\cdot,\cdot]} are the fourth-order contractions of
#' [gas_levinson4()], all three of which are carried, with the third-order
#' contractions along each direction that they need. The result is symmetric
#' in the three directions, and identically zero for \eqn{q \le 4}.
#'
#' @param pacf A numeric vector of partial autocorrelations in
#'   \eqn{(-1, 1)}.
#' @param v,w,x The three directions, each as long as `pacf`.
#'
#' @return A list of one `q` by `q` matrix per coefficient.
#'
#' @seealso [gas_levinson4()]
#'
#' @keywords internal
gas_levinson5 <- function(pacf, v, w, x) {
  q <- length(pacf)
  if (q == 0L) return(list())
  dirs <- list(v, w, x)
  pairs <- list(c(1L, 2L), c(1L, 3L), c(2L, 3L))
  phi <- numeric(0)
  jac <- matrix(0, 0, q)
  hes <- list()
  tt <- list(list(), list(), list())
  qq <- list(list(), list(), list())
  ff <- list()
  zero <- function(k) replicate(k, matrix(0, q, q), simplify = FALSE)
  for (k in seq_len(q)) {
    new <- numeric(k)
    njac <- matrix(0, k, q)
    nhes <- zero(k)
    ntt <- list(zero(k), zero(k), zero(k))
    nqq <- list(zero(k), zero(k), zero(k))
    nff <- zero(k)
    new[k] <- pacf[k]
    njac[k, k] <- 1
    if (k > 1L) {
      rev_idx <- rev(seq_len(k - 1L))
      new[seq_len(k - 1L)] <- phi - pacf[k] * phi[rev_idx]
      njac[seq_len(k - 1L), ] <- jac - pacf[k] * jac[rev_idx, , drop = FALSE]
      njac[seq_len(k - 1L), k] <- njac[seq_len(k - 1L), k] - phi[rev_idx]
      for (i in seq_len(k - 1L)) {
        r <- rev_idx[i]
        sym_k <- function(m, u) {
          m[k, ] <- m[k, ] - u
          m[, k] <- m[, k] - u
          m
        }
        nhes[[i]] <- sym_k(hes[[i]] - pacf[k] * hes[[r]], jac[r, ])
        for (d in 1:3) {
          ntt[[d]][[i]] <- sym_k(tt[[d]][[i]] - pacf[k] * tt[[d]][[r]] -
                                   dirs[[d]][[k]] * hes[[r]],
                                 as.numeric(hes[[r]] %*% dirs[[d]]))
        }
        for (pp in 1:3) {
          a <- pairs[[pp]][1L]
          b <- pairs[[pp]][2L]
          nqq[[pp]][[i]] <- sym_k(qq[[pp]][[i]] - pacf[k] * qq[[pp]][[r]] -
                                    dirs[[a]][[k]] * tt[[b]][[r]] -
                                    dirs[[b]][[k]] * tt[[a]][[r]],
                                  as.numeric(tt[[a]][[r]] %*% dirs[[b]]))
        }
        nff[[i]] <- sym_k(ff[[i]] - pacf[k] * ff[[r]] -
                            v[[k]] * qq[[3L]][[r]] - w[[k]] * qq[[2L]][[r]] -
                            x[[k]] * qq[[1L]][[r]],
                          as.numeric(qq[[1L]][[r]] %*% x))
      }
    }
    phi <- new
    jac <- njac
    hes <- nhes
    tt <- ntt
    qq <- nqq
    ff <- nff
  }
  ff
}


#' The Levinson-Durbin Map's Derivative of Any Order, Contracted
#'
#' @description
#' Returns, per coefficient, the derivative of order `length(dirs) + 2` of the
#' map, contracted against the directions in `dirs` and leaving two slots
#' free, as a `q` by `q` matrix. With no direction this is the Hessian.
#'
#' @param rho The partial autocorrelations.
#' @param dirs A list of zero to three vectors as long as `rho`.
#'
#' @return A list of one matrix per coefficient.
#'
#' @keywords internal
.gas_levinson_contract <- function(rho, dirs) {
  switch(length(dirs) + 1L,
         gas_levinson2(rho)$hessian,
         gas_levinson3(rho, dirs[[1L]]),
         gas_levinson4(rho, dirs[[1L]], dirs[[2L]]),
         gas_levinson5(rho, dirs[[1L]], dirs[[2L]], dirs[[3L]]))
}


#' The Chart's Derivatives of Any Order up to Five, in a Set of Directions
#'
#' @description
#' Returns the chart's derivative of order `length(dirs) + 2` in the term's
#' base coordinates, contracted against the directions in `dirs` and leaving
#' two slots free: one matrix for the level, one per score loading and one per
#' autoregressive coefficient. With zero, one and two directions it equals
#' the second, third and fourth orders of [.gas_chart_derivs()],
#' [.gas_chart_derivs3()] and [.gas_chart_derivs4()].
#'
#' @details
#' The level and the loadings ride scalar links, so their derivative is a
#' single diagonal entry, the link's derivative of the order times the
#' directions' components there. The autoregressive coefficients are the
#' Levinson-Durbin map \eqn{\xi} read at \eqn{\rho = h(z)}, with \eqn{h}
#' acting componentwise. Faa di Bruno's formula over the set partitions
#' \eqn{\pi} of the slots (the directions and the two free slots) gives
#' \deqn{D^n(\xi\circ h) = \sum_{\pi}\xi^{(|\pi|)}\big[h^{(|\beta|)}
#'   \circ \textstyle\prod_{s\in\beta}s\big]_{\beta\in\pi},}
#' since the inner map is diagonal. A block holding no free slot pushes a
#' vector, a block holding one free slot scales that slot, and a block holding
#' both puts a diagonal. Each \eqn{\xi^{(j)}} is read from
#' [.gas_levinson_contract()] against the pushed vectors.
#'
#' @param zeta The term's base parameters on the unconstrained scale.
#' @param p,q The score and autoregressive orders.
#' @param links The links, as [term_links()] gives them.
#' @param dirs A list of zero to three directions in the coordinates of
#'   `zeta`.
#'
#' @return A list with `omega`, `a` and `b`: a matrix, a list of `max(p, 1)`
#'   matrices and a list of `q` matrices over the base coordinates.
#'
#' @seealso [gas_levinson5()]
#'
#' @keywords internal
.gas_chart_contract <- function(zeta, p, q, links, dirs) {
  np <- length(zeta)
  nm <- names(zeta)
  nd <- length(dirs)
  ord <- nd + 2L
  i_om <- 1L
  i_a <- if (p > 0L) 1L + seq_len(p) else integer(0)
  i_pa <- if (q > 0L) 1L + p + seq_len(q) else integer(0)
  dlk <- list(linkfunctions7::dlinkinv, linkfunctions7::d2linkinv,
              linkfunctions7::d3linkinv, linkfunctions7::d4linkinv,
              linkfunctions7::d5linkinv)
  scalar <- function(i) {
    h <- matrix(0, np, np)
    val <- dlk[[ord]](links[[nm[i]]], zeta[[i]])
    for (d in dirs) val <- val * d[[i]]
    h[i, i] <- val
    h
  }
  omega <- scalar(i_om)
  a <- lapply(seq_len(max(p, 1L)), function(i) {
    if (p > 0L) scalar(i_a[i]) else matrix(0, np, np)
  })
  b <- list()
  if (q > 0L) {
    lk <- links[nm[i_pa]]
    z <- zeta[i_pa]
    rho <- .gas_links_at(linkfunctions7::linkinv, lk, z)
    k <- lapply(1:5, function(o) .gas_links_at(dlk[[o]], lk, z))
    dv <- lapply(dirs, function(d) d[i_pa])
    jac <- gas_levinson2(rho)$jacobian
    # the slots: the directions 1..nd, then the free slots nd + 1, nd + 2
    fa <- nd + 1L
    fb <- nd + 2L
    acc <- replicate(q, matrix(0, q, q), simplify = FALSE)
    for (part in numericals7::set_partitions(ord)) {
      fixed <- list()
      ca <- cb <- cab <- NULL
      for (blk in part) {
        dirs_in <- blk[blk <= nd]
        pushed <- k[[length(blk)]]
        for (s in dirs_in) pushed <- pushed * dv[[s]]
        has_a <- fa %in% blk
        has_b <- fb %in% blk
        if (has_a && has_b) cab <- pushed
        else if (has_a) ca <- pushed
        else if (has_b) cb <- pushed
        else fixed[[length(fixed) + 1L]] <- pushed
      }
      j <- length(part)
      for (cf in seq_len(q)) {
        if (!is.null(cab)) {
          # both free slots in one block: one free slot of xi^(j) left
          vec <- if (j == 1L) jac[cf, ] else {
            mm <- .gas_levinson_contract(rho, fixed[-length(fixed)])[[cf]]
            as.numeric(mm %*% fixed[[length(fixed)]])
          }
          acc[[cf]] <- acc[[cf]] + diag(vec * cab, nrow = q)
        } else {
          mm <- .gas_levinson_contract(rho, fixed)[[cf]]
          acc[[cf]] <- acc[[cf]] + (ca * mm) %*% diag(cb, nrow = q)
        }
      }
    }
    for (cf in seq_len(q)) {
      h <- matrix(0, np, np)
      h[i_pa, i_pa] <- acc[[cf]]
      b[[cf]] <- h
    }
  }
  list(omega = omega, a = a, b = b)
}


#' A Jet over a Set of Directions, and Its Product Rule
#'
#' @description
#' `.gas_jetn()` holds the derivatives of one factor of the recursion along
#' every subset of up to three directions, each with two free slots, and
#' `.gas_prodn()` returns the derivative of the product of two such factors
#' along all the directions, with two free slots.
#'
#' @details
#' A jet has the value, the first derivative and one matrix per subset
#' \eqn{S} of the directions, the derivative of order \eqn{|S| + 2} along
#' \eqn{S} with two free slots. The derivative along \eqn{S} with one free
#' slot is the matrix of \eqn{S} without one of its directions applied to
#' that direction, and with no free slot it is that vector applied to another
#' direction. The product rule deals the directions and the two free slots to
#' the two factors:
#' \deqn{D_S(AB) = \sum_{T\subseteq S}\big(A_{T,2}B_{S\setminus T,0}
#'   + A_{T,0}B_{S\setminus T,2} + A_{T,1}B_{S\setminus T,1}^{\top}
#'   + B_{S\setminus T,1}A_{T,1}^{\top}\big),}
#' which is 16 terms with two directions, as [.gas_prod4()] writes them, and
#' 32 with three.
#'
#' @param v0 The value.
#' @param d1 The first derivative, a vector.
#' @param mats A list of matrices named by subset: `"s"` for the second
#'   derivative, `"s1"`, `"s2"`, `"s12"` and so on, the digits being the
#'   positions of the directions in increasing order.
#' @param dirs The directions.
#' @param A,B Two jets over the same directions.
#'
#' @return `.gas_jetn()` a jet; `.gas_prodn()` a matrix.
#'
#' @keywords internal
.gas_jetn <- function(v0, d1, mats, dirs) {
  list(v = v0, d1 = d1, mats = mats, dirs = dirs)
}

#' @rdname dot-gas_jetn
.gas_subsets <- function(nd) {
  out <- list(integer(0))
  for (s in seq_len(nd)) out <- c(out, lapply(out, function(t) c(t, s)))
  out[order(lengths(out))]
}

#' @rdname dot-gas_jetn
.gas_jet_get <- function(A, S, free) {
  key <- paste0("s", paste(sort(S), collapse = ""))
  if (free == 2L) return(A$mats[[key]])
  if (free == 1L) {
    if (!length(S)) return(A$d1)
    last <- S[length(S)]
    return(as.numeric(.gas_jet_get(A, S[-length(S)], 2L) %*% A$dirs[[last]]))
  }
  if (!length(S)) return(A$v)
  last <- S[length(S)]
  sum(.gas_jet_get(A, S[-length(S)], 1L) * A$dirs[[last]])
}

#' @rdname dot-gas_jetn
.gas_prodn <- function(A, B) {
  nd <- length(A$dirs)
  all <- seq_len(nd)
  out <- 0
  for (Tset in .gas_subsets(nd)) {
    R <- setdiff(all, Tset)
    a1 <- .gas_jet_get(A, Tset, 1L)
    b1 <- .gas_jet_get(B, R, 1L)
    out <- out + .gas_jet_get(A, Tset, 2L) * .gas_jet_get(B, R, 0L) +
      .gas_jet_get(A, Tset, 0L) * .gas_jet_get(B, R, 2L) +
      outer(a1, b1) + outer(b1, a1)
  }
  out
}


#' @rdname dot-gas_jetn
#' @param S The subset of the directions along which to differentiate, by
#'   position.
.gas_prodS <- function(A, B, S) {
  out <- 0
  for (Tset in .gas_subsets_of(S)) {
    R <- setdiff(S, Tset)
    a1 <- .gas_jet_get(A, Tset, 1L)
    b1 <- .gas_jet_get(B, R, 1L)
    out <- out + .gas_jet_get(A, Tset, 2L) * .gas_jet_get(B, R, 0L) +
      .gas_jet_get(A, Tset, 0L) * .gas_jet_get(B, R, 2L) +
      outer(a1, b1) + outer(b1, a1)
  }
  out
}

#' @rdname dot-gas_jetn
.gas_subsets_of <- function(S) {
  lapply(.gas_subsets(length(S)), function(i) S[i])
}

#' @rdname dot-gas_jetn
.gas_jet_key <- function(S) paste0("s", paste(sort(S), collapse = ""))

#' The Product of Two Jets, as a Jet
#'
#' @description
#' Returns the jet of \eqn{AB}: its value, its first derivative and its
#' matrix along every subset of the directions.
#'
#' @param A,B Two jets over the same directions.
#'
#' @return A jet.
#'
#' @keywords internal
.gas_jet_times <- function(A, B) {
  nd <- length(A$dirs)
  mats <- list()
  for (S in .gas_subsets(nd)) mats[[.gas_jet_key(S)]] <- .gas_prodS(A, B, S)
  .gas_jetn(A$v * B$v, A$v * B$d1 + B$v * A$d1, mats, A$dirs)
}

#' @rdname dot-gas_jet_times
#' @param w A number, the weight of the second jet.
.gas_jet_add <- function(A, B, w = 1) {
  mats <- A$mats
  for (k in names(mats)) mats[[k]] <- mats[[k]] + w * B$mats[[k]]
  .gas_jetn(A$v + w * B$v, A$d1 + w * B$d1, mats, A$dirs)
}

#' @rdname dot-gas_jet_times
#' @param v0,d1 A value and a first derivative.
#' @param m The number of coordinates.
#' @param dirs The directions.
.gas_jet_const <- function(v0, d1, m, dirs) {
  mats <- list()
  for (S in .gas_subsets(length(dirs))) {
    mats[[.gas_jet_key(S)]] <- matrix(0, m, m)
  }
  .gas_jetn(v0, d1, mats, dirs)
}


#' The Jet of a Function of Several Predictors
#'
#' @description
#' Composes the derivatives of a function \eqn{u(e_1, \dots, e_P)} of the
#' predictors of every equation with the jets of those predictors, by Faa di
#' Bruno's formula over the set partitions of the slots: the directions of a
#' subset and the two free slots. `.gas_score_jet()` is the case of one
#' filtered equation, whose predictor alone has derivatives above the first.
#'
#' @details
#' A partition with \eqn{k} blocks contributes
#' \eqn{\sum_{r_1,\dots,r_k} u_{r_1\dots r_k}\prod_\beta
#' \partial^{|\beta|}e_{r_\beta}[\beta]}, where each factor is the
#' derivative of one predictor along the slots of its block: a number for a
#' block without a free slot, a row for a block with one, a matrix for a
#' block with both. The coefficients \eqn{u_{r_1\dots r_k}} are read from
#' `ell`, the arrays of the driving quantity's derivatives in the predictors
#' at the observation (`H`, `D3`, ..., `D6` of [term_fifth()], or `U1` to
#' `U5` of a joint filter).
#'
#' @param ell A list whose element \eqn{k} is the array of the derivatives
#'   of order \eqn{k}, with \eqn{k} dimensions of `P` each.
#' @param s0 The value.
#' @param E A list of `P` jets, one per equation's predictor.
#' @param ej The jet of the filtered predictor.
#' @param V A `P` by `m` matrix of the predictors' first derivatives (the
#'   row `ap` is replaced by the jet's).
#' @param ap The index of the filtered equation.
#'
#' @return A jet.
#'
#' @keywords internal
.gas_compose <- function(ell, s0, E) {
  dirs <- E[[1L]]$dirs
  nd <- length(dirs)
  P <- length(E)
  m <- length(E[[1L]]$d1)
  mats <- list()
  for (S in .gas_subsets(nd)) {
    ns <- length(S)
    slots <- c(S, -1L, -2L)
    acc <- matrix(0, m, m)
    for (part in numericals7::set_partitions(ns + 2L)) {
      k <- length(part)
      cvec <- list()
      free <- list()
      for (bi in seq_along(part)) {
        sl <- slots[part[[bi]]]
        dset <- sl[sl > 0L]
        nfree <- sum(sl < 0L)
        if (nfree == 0L) {
          cvec[[bi]] <- vapply(E, function(e) .gas_jet_get(e, dset, 0L),
                               numeric(1))
        } else if (nfree == 1L) {
          rows <- t(vapply(E, function(e) .gas_jet_get(e, dset, 1L),
                           numeric(m)))
          free[[length(free) + 1L]] <- list(block = bi, rows = rows,
                                            a = -1L %in% sl)
        } else {
          free[[length(free) + 1L]] <- list(
            block = bi, mats = lapply(E, function(e) .gas_jet_get(e, dset, 2L)))
        }
      }
      fb <- vapply(free, function(f) f$block, integer(1))
      sc <- setdiff(seq_len(k), fb)
      arr <- aperm(array(ell[[k]], rep(P, k)), c(sc, fb))
      for (bi in sc) {
        arr <- array(arr, c(P, length(arr) / P))
        arr <- as.numeric(crossprod(cvec[[bi]], arr))
      }
      if (length(free) == 1L) {
        for (r in seq_len(P)) {
          if (arr[r] != 0) acc <- acc + arr[r] * free[[1L]]$mats[[r]]
        }
      } else {
        cmat <- matrix(arr, P, P)
        if (free[[1L]]$a) {
          acc <- acc + crossprod(free[[1L]]$rows, cmat %*% free[[2L]]$rows)
        } else {
          acc <- acc + crossprod(free[[2L]]$rows, t(cmat) %*% free[[1L]]$rows)
        }
      }
    }
    mats[[.gas_jet_key(S)]] <- acc / 2 + t(acc) / 2
  }
  Vm <- t(vapply(E, function(e) e$d1, numeric(m)))
  .gas_jetn(s0, as.numeric(crossprod(Vm, as.numeric(ell[[1L]]))), mats, dirs)
}

#' @rdname dot-gas_compose
.gas_score_jet <- function(ell, s0, ej, V, ap) {
  m <- ncol(V)
  E <- lapply(seq_len(nrow(V)), function(r) {
    if (r == ap) ej else .gas_jet_const(0, as.numeric(V[r, ]), m, ej$dirs)
  })
  .gas_compose(ell, s0, E)
}


#' The Chart's Jets in the Caller's Unknowns
#'
#' @description
#' Returns the jets of the level, of each score loading, of each
#' autoregressive coefficient and of the starting level, over the given
#' directions, in the caller's unknowns.
#'
#' @details
#' The chart's derivatives are those of [.gas_chart_contract()] in the base
#' coordinates, lifted to the unknowns through the columns `zbase` that hold
#' the base coordinates. The starting level \eqn{f_0 = \omega/(1 - S)}, with
#' \eqn{S} the sum of the autoregressive coefficients, is read off the fixed
#' point \eqn{f_0 = \omega + S f_0}: along each subset of the directions,
#' taken by increasing size, its matrix is \eqn{(1 - S)^{-1}} times the
#' level's plus the product rule of \eqn{S f_0} without the term that
#' carries every derivative on \eqn{f_0}.
#'
#' @param zeta The base parameters on the unconstrained scale.
#' @param p,q The orders.
#' @param links The links.
#' @param zbase The columns of the unknowns that hold the base coordinates.
#' @param m The number of unknowns.
#' @param dirs The directions, in the unknowns.
#'
#' @return A list with `omega`, `a` (a list), `b` (a list) and `f0`, jets.
#'
#' @keywords internal
.gas_chart_jets <- function(zeta, p, q, links, zbase, m, dirs) {
  nb <- length(zbase)
  A <- matrix(0, nb, m)
  A[cbind(seq_len(nb), zbase)] <- 1
  lift <- function(v) as.numeric(crossprod(A, v))
  lift2 <- function(h) crossprod(A, h %*% A)
  ch <- .gas_chart_derivs(zeta, p, q, links)
  nd <- length(dirs)
  subs <- .gas_subsets(nd)
  per <- lapply(subs, function(S) {
    .gas_chart_contract(zeta, p, q, links,
                        lapply(dirs[S], function(d) d[zbase]))
  })
  jet_of <- function(v0, d1, pick) {
    mats <- list()
    for (k in seq_along(subs)) {
      mats[[.gas_jet_key(subs[[k]])]] <- lift2(pick(per[[k]]))
    }
    .gas_jetn(v0, lift(d1), mats, dirs)
  }
  omj <- jet_of(ch$omega, ch$d_omega, function(r) r$omega)
  aj <- lapply(seq_len(p), function(i)
    jet_of(ch$a[[i]], ch$d_a[[i]], function(r) r$a[[i]]))
  bj <- lapply(seq_len(q), function(j)
    jet_of(ch$b[[j]], ch$d_b[[j]], function(r) r$b[[j]]))
  Sj <- .gas_jet_const(0, numeric(m), m, dirs)
  for (j in seq_len(q)) Sj <- .gas_jet_add(Sj, bj[[j]])
  if (abs(1 - Sj$v) < 1e-10) {
    stop("the autoregressive polynomial is at the unit root; the filter has no starting level.",
         call. = FALSE)
  }
  cst <- 1 / (1 - Sj$v)
  f0 <- omj$v * cst
  f0j <- .gas_jet_const(f0, cst * (omj$d1 + Sj$d1 * f0), m, dirs)
  for (S in subs) {
    key <- .gas_jet_key(S)
    f0j$mats[[key]] <- cst * (omj$mats[[key]] + .gas_prodS(Sj, f0j, S))
  }
  list(omega = omj, a = aj, b = bj, f0 = f0j)
}


#' The Score-Driven Recursion's Derivatives to the Fifth Order, R Route
#'
#' @description
#' Propagates the jets of the filtered predictor over three directions
#' through the recursion and returns the contractions of [term_fifth()].
#' `.gas_fifth_scalar()` serves a term without subformulas, whose chart is the
#' same at every observation, and `.gas_fifth_sub()` a term with subformulas,
#' whose chart is read per observation; both run `.gas_fifth_run()`.
#'
#' @details
#' At each observation the level is
#' \eqn{f_t = \omega + \sum_i a_i s_{t-i} + \sum_j b_j f_{t-j}}, and its jet
#' is the sum of the products of jets ([.gas_jet_times()]). The score's jet
#' is composed from the family's derivative arrays by [.gas_score_jet()]. The
#' orders below the fifth are carried along every subset of the directions,
#' so this route is a second, independent computation of what
#' [term_curvature()], [term_third()] and [term_fourth()] return. The
#' starting level of a group is read off the fixed point at the group's first
#' observation, as at the lower orders.
#'
#' @inheritParams term_fifth
#' @param dirs The three directions.
#' @param chart_at A function of a row returning the jets `omega`, `a` and
#'   `b` there.
#' @param ch The chart's jets at the group's first observation.
#' @param m The number of unknowns.
#'
#' @return `.gas_start_jet()` the jet of the starting level; the others as
#'   [term_fifth()].
#'
#' @keywords internal
.gas_fifth_run <- function(term, eta, score, g, seed, blocks_data, dirs,
                           chart_at) {
  bp <- term@build_info
  m <- ncol(seed)
  n <- bp$n
  p <- term@p
  q <- term@q
  static <- .gas_static(term)
  d0 <- if (static) .gas_static_start(term, eta) else NULL
  bd <- blocks_data
  npar <- length(bd$Vs)
  ap <- as.integer(bd$ap)
  arrs <- list(bd$H, bd$D3, bd$D4, bd$D5, bd$D6)
  ell_at <- function(row) {
    lapply(seq_len(5L), function(k) array(arrs[[k]][row, ], rep(npar, k)))
  }
  V_at <- function(row) {
    t(vapply(seq_len(npar), function(r) {
      if (r == ap) numeric(m) else as.numeric(bd$Vs[[r]][row, ])
    }, numeric(m)))
  }
  pairs <- list(c(1L, 2L), c(1L, 3L), c(2L, 3L))
  D <- matrix(0, n, m)
  dphi <- lapply(1:3, function(d) matrix(0, n, m))
  dpsi <- lapply(1:3, function(d) matrix(0, n, m))
  dom <- matrix(0, n, m)
  W <- matrix(0, m, m)
  zero <- .gas_jet_const(0, numeric(m), m, dirs)
  for (l in seq_along(bp$order)) {
    rows <- bp$order[[l]]
    k <- length(rows)
    f0j <- if (static) {
      .gas_jet_const(d0[[l]], -seed[rows[1L], ], m, dirs)
    } else {
      .gas_start_jet(chart_at(rows[1L]), m, dirs)
    }
    fj <- vector("list", k)
    sj <- vector("list", k)
    for (t in seq_len(k)) {
      row <- rows[t]
      if (static && t == 1L) {
        ft <- f0j
      } else {
        ch <- chart_at(row)
        ft <- ch$omega
        for (i in seq_len(p)) {
          lag <- t - i
          ft <- .gas_jet_add(ft, .gas_jet_times(ch$a[[i]],
                                                if (lag >= 1L) sj[[lag]] else zero))
        }
        for (j in seq_len(q)) {
          lag <- t - j
          ft <- .gas_jet_add(ft, .gas_jet_times(ch$b[[j]],
                                                if (lag >= 1L) fj[[lag]] else f0j))
        }
      }
      fj[[t]] <- ft
      ej <- ft
      ej$v <- ft$v + eta[row]
      ej$d1 <- ft$d1 + seed[row, ]
      sj[[t]] <- .gas_score_jet(ell_at(row), score(ej$v, row), ej,
                                V_at(row), ap)
      D[row, ] <- ej$d1
      for (d in 1:3) dphi[[d]][row, ] <- .gas_jet_get(ej, d, 1L)
      for (pp in 1:3) dpsi[[pp]][row, ] <- .gas_jet_get(ej, pairs[[pp]], 1L)
      dom[row, ] <- .gas_jet_get(ej, 1:3, 1L)
      W <- W + g[row] * ft$mats$s123
    }
  }
  W <- (W + t(W)) / 2
  list(jacobian = D, dphi = dphi, dpsi = dpsi, domega = dom, curvature = W)
}

#' @rdname dot-gas_fifth_run
.gas_start_jet <- function(ch, m, dirs) {
  Sj <- .gas_jet_const(0, numeric(m), m, dirs)
  for (bj in ch$b) Sj <- .gas_jet_add(Sj, bj)
  if (abs(1 - Sj$v) < 1e-10) {
    stop("the autoregressive polynomial is at the unit root; the filter has no starting level.",
         call. = FALSE)
  }
  cst <- 1 / (1 - Sj$v)
  omj <- ch$omega
  f0 <- omj$v * cst
  f0j <- .gas_jet_const(f0, cst * (omj$d1 + Sj$d1 * f0), m, dirs)
  for (S in .gas_subsets(length(dirs))) {
    key <- .gas_jet_key(S)
    f0j$mats[[key]] <- cst * (omj$mats[[key]] + .gas_prodS(Sj, f0j, S))
  }
  f0j
}

#' @rdname dot-gas_fifth_run
.gas_fifth_scalar <- function(term, eta, y, score, psi, g, seed, blocks_data,
                              dirs) {
  nm <- term_params(term)
  links <- term_links(term)
  psiv <- unlist(psi[nm])
  base <- .gas_base_params(term@p, term@q)
  nb <- length(base)
  zeta <- vapply(base, function(j)
    linkfunctions7::linkfun(links[[j]], psiv[[j]]), numeric(1))
  m <- ncol(seed)
  np <- length(nm)
  zcol <- m - np + seq_len(np)
  ch <- .gas_chart_jets(zeta, term@p, term@q, links, zcol[seq_len(nb)], m,
                        dirs)
  .gas_fifth_run(term, eta, score, g, seed, blocks_data, dirs,
                 function(row) ch)
}

#' @rdname dot-gas_fifth_run
.gas_fifth_sub <- function(term, eta, y, score, psi, g, seed, blocks_data,
                           dirs) {
  nm <- term_params(term)
  m <- ncol(seed)
  zcol <- m - length(nm) + seq_along(nm)
  chart_at <- .gas_sub_chart_fun(term, unlist(psi[nm]), zcol, m, dirs)
  .gas_fifth_run(term, eta, score, g, seed, blocks_data, dirs, chart_at)
}

#' @rdname dot-gas_fifth_run
#' @param tm A built score-driven term, with or without subformulas.
#' @param psiv Its parameters in the order of [term_params()]: the parameter
#'   scale for a parameter without a subformula, the coefficients for a
#'   developed one.
#' @param cols The columns of those parameters among the unknowns.
#' @details `.gas_sub_chart_fun()` returns the chart at a row for a term
#'   whose values are read per observation: each base parameter and each cross
#'   loading is \eqn{h(z^\top u)}, with \eqn{z} its design row among the
#'   unknowns (a unit vector for a parameter without a subformula), and the
#'   autoregressive coefficients come from [.gas_map_jets()]. The result at a
#'   row carries `x`, the cross loadings' jets in the order of the term's
#'   pairs (partner, lag).
.gas_sub_chart_fun <- function(tm, psiv, cols, m, dirs) {
  bp <- tm@build_info
  p <- tm@p
  q <- tm@q
  lay <- .gas_sub_layout(tm)
  slots <- names(lay$idx)
  sub <- bp$sub
  u <- as.numeric(psiv)
  for (j in slots) {
    if (is.null(sub[[j]])) {
      u[lay$idx[[j]]] <- linkfunctions7::linkfun(.gas_param_link(tm, j),
                                                 u[lay$idx[[j]]])
    }
  }
  dl <- list(linkfunctions7::linkinv, linkfunctions7::dlinkinv,
             linkfunctions7::d2linkinv, linkfunctions7::d3linkinv,
             linkfunctions7::d4linkinv, linkfunctions7::d5linkinv)
  rowinfo <- lapply(stats::setNames(slots, slots), function(j) {
    idx <- lay$idx[[j]]
    lk <- .gas_param_link(tm, j)
    if (!is.null(sub[[j]])) {
      Z <- as.matrix(sub[[j]]$Z)
      e <- as.numeric(Z %*% u[idx])
    } else {
      Z <- matrix(1, bp$n, 1L)
      e <- rep(u[[idx]], bp$n)
    }
    list(cols = cols[idx], Z = Z,
         K = matrix(vapply(dl, function(f) f(lk, e), numeric(bp$n)), bp$n))
  })
  scalar_jet <- function(info, r) {
    z <- numeric(m)
    z[info$cols] <- info$Z[r, ]
    K <- info$K[r, ]
    mats <- list()
    for (S in .gas_subsets(length(dirs))) {
      f <- K[3L + length(S)]
      for (s in S) f <- f * sum(z * dirs[[s]])
      mats[[.gas_jet_key(S)]] <- f * outer(z, z)
    }
    .gas_jetn(K[1L], K[2L] * z, mats, dirs)
  }
  cn <- .gas_cross_params(tm)
  function(r) {
    a <- lapply(seq_len(p), function(i)
      scalar_jet(rowinfo[[paste0("kappa", i)]], r))
    b <- list()
    if (q > 0L) {
      pj <- paste0("pacf", seq_len(q))
      zs <- t(vapply(pj, function(j) {
        z <- numeric(m)
        z[rowinfo[[j]]$cols] <- rowinfo[[j]]$Z[r, ]
        z
      }, numeric(m)))
      K <- t(vapply(pj, function(j) rowinfo[[j]]$K[r, ], numeric(6L)))
      b <- .gas_map_jets(K, zs, dirs)
    }
    list(omega = scalar_jet(rowinfo$omega, r), a = a, b = b,
         x = lapply(cn, function(j) scalar_jet(rowinfo[[j]], r)))
  }
}


#' The Autoregressive Coefficients' Jets from Rank-One Partial Autocorrelations
#'
#' @description
#' Returns the jets of the coefficients \eqn{\xi_j(\rho)} of the
#' Levinson-Durbin map, where each partial autocorrelation is
#' \eqn{\rho_k = h_k(z_k^\top u)}.
#'
#' @details
#' The derivative of order \eqn{n} of \eqn{\rho_k} is
#' \eqn{h_k^{(n)} z_k^{\otimes n}}, so Faa di Bruno's formula over the set
#' partitions of the slots reads, for each block, a vector over the partial
#' autocorrelations (a block without a free slot), a \eqn{q \times m} matrix
#' (a block with one free slot) or a weighted sum of \eqn{z_k z_k^\top} (a
#' block with both). The map's derivatives come from
#' [.gas_levinson_contract()]. With \eqn{z_k} a unit vector this is
#' [.gas_chart_contract()] lifted to the unknowns.
#'
#' @param K A `q` by 6 matrix: the inverse links and their derivatives of
#'   orders one to five, per partial autocorrelation.
#' @param zs A `q` by `m` matrix of the rows \eqn{z_k}.
#' @param dirs The directions.
#'
#' @return A list of `q` jets.
#'
#' @keywords internal
.gas_map_jets <- function(K, zs, dirs) {
  q <- nrow(K)
  m <- ncol(zs)
  rho <- K[, 1L]
  ld <- gas_levinson2(rho)
  jac <- ld$jacobian
  zd <- lapply(dirs, function(d) as.numeric(zs %*% d))
  mats <- replicate(q, list(), simplify = FALSE)
  for (S in .gas_subsets(length(dirs))) {
    ns <- length(S)
    slots <- c(S, -1L, -2L)
    acc <- replicate(q, matrix(0, m, m), simplify = FALSE)
    for (part in numericals7::set_partitions(ns + 2L)) {
      fixed <- list()
      Ca <- Cb <- cab <- NULL
      for (blk in part) {
        sl <- slots[blk]
        co <- K[, 1L + length(blk)]
        for (s in sl[sl > 0L]) co <- co * zd[[s]]
        if (all(c(-1L, -2L) %in% sl)) cab <- co
        else if (-1L %in% sl) Ca <- co * zs
        else if (-2L %in% sl) Cb <- co * zs
        else fixed[[length(fixed) + 1L]] <- co
      }
      J <- length(part)
      for (cf in seq_len(q)) {
        if (!is.null(cab)) {
          vec <- if (J == 1L) jac[cf, ] else {
            mm <- .gas_levinson_contract(rho, fixed[-length(fixed)])[[cf]]
            as.numeric(mm %*% fixed[[length(fixed)]])
          }
          acc[[cf]] <- acc[[cf]] + crossprod(zs, (vec * cab) * zs)
        } else {
          mm <- .gas_levinson_contract(rho, fixed)[[cf]]
          acc[[cf]] <- acc[[cf]] + crossprod(Ca, mm %*% Cb)
        }
      }
    }
    for (cf in seq_len(q)) mats[[cf]][[.gas_jet_key(S)]] <- acc[[cf]]
  }
  lapply(seq_len(q), function(cf) {
    .gas_jetn(ld$phi[cf], as.numeric(crossprod(zs, jac[cf, ] * K[, 2L])),
              mats[[cf]], dirs)
  })
}


#' @title The Score-Driven Recursion's Fifth Derivative
#' @name term_fifth.GasTerm
#'
#' @description
#' Computes the fifth derivative of the filtered predictor in the caller's
#' unknowns, contracted against the caller's weights and three directions,
#' together with the lower orders' contractions along the same directions.
#'
#' @details
#' The jets of the level, of the score and of the chart are carried along
#' every subset of the three directions, and each product of the recursion is
#' differentiated by the product rule over those subsets. The chart's
#' derivatives come from [.gas_chart_contract()] and [gas_levinson5()], the
#' score's from the family's arrays by [.gas_score_jet()].
#'
#' @inheritParams term_fifth
#' @return As [term_fifth()].
#' @keywords internal
S7::method(term_fifth, GasTerm) <- function(term, eta, y, score, psi, g, seed,
                                            blocks_data, directions, ...) {
  bp <- term@build_info
  if (!length(bp)) {
    stop("the term has not been built; call term_build(term, data) first.",
         call. = FALSE)
  }
  if (!is.list(directions) || length(directions) != 3L) {
    stop("'directions' must be a list of three vectors, one per direction.",
         call. = FALSE)
  }
  need <- c("H", "D3", "D4", "D5", "D6", "Vs", "ap")
  if (!is.list(blocks_data) || !all(need %in% names(blocks_data))) {
    stop("'blocks_data' must be a list with H, D3, D4, D5, D6, Vs and ap.",
         call. = FALSE)
  }
  seed <- as.matrix(seed)
  m <- ncol(seed)
  if (nrow(seed) != bp$n) {
    stop(sprintf("'seed' must have one row per observation (%d).", bp$n),
         call. = FALSE)
  }
  g <- as.numeric(g)
  if (length(g) != bp$n) {
    stop(sprintf("'g' must give one value per observation (%d).", bp$n),
         call. = FALSE)
  }
  dirs <- lapply(directions, function(d) {
    d <- as.numeric(d)
    if (length(d) != m) {
      stop(sprintf("each direction must have one value per unknown (%d).", m),
           call. = FALSE)
    }
    d
  })
  if (!is.null(bp$sub)) {
    return(.gas_fifth_sub(term, eta, y, score, psi, g, seed, blocks_data,
                          dirs))
  }
  .gas_fifth_scalar(term, eta, y, score, psi, g, seed, blocks_data, dirs)
}


#' @title The Fifth Derivative of a Joint Score-Driven Filter
#' @name term_fifth.GasJoint
#'
#' @description
#' Computes the fifth derivative of the filtered predictors contracted
#' against the caller's weights and three directions, with the lower orders'
#' contractions along the same directions, through jets over every subset of
#' the directions.
#'
#' @details
#' The levels' jets are sums of products of the chart's jets with the jets
#' of earlier driving scores and levels, as at the lower orders of
#' [term_fourth()]. The jet of each driving score is composed from the jets
#' of the predictors of every equation by [.gas_compose()], with the
#' driving score's derivatives in those predictors given as data in
#' `blocks_data`: `filtered`, `u`, `U1`, `U2`, `U3`, `U4` as for the lower
#' orders, and `U5`, an array with the observation, the driving score and five
#' predictor indices.
#'
#' @inheritParams term_fifth
#' @param g The weights, one column per filtered parameter.
#' @param seed A list of one matrix of static rows per equation.
#'
#' @return A list with `jacobian`, `dphi`, `dpsi` and `domega`, each with one
#'   element per filtered parameter (`dphi` and `dpsi` lists of three), and
#'   `curvature`.
#'
#' @keywords internal
S7::method(term_fifth, GasJoint) <- function(term, eta, y, score, psi, g, seed,
                                             blocks_data, directions, ...) {
  if (!is.list(directions) || length(directions) != 3L) {
    stop("'directions' must be a list of three vectors, one per direction.",
         call. = FALSE)
  }
  if (is.null(blocks_data$U5)) {
    stop("'blocks_data' must carry U5 for the fifth order.", call. = FALSE)
  }
  if (.gas_joint_has_sub(term)) {
    return(.gas_joint_fifth(term, psi, g, seed, blocks_data, directions,
                            eta = eta, sub = TRUE))
  }
  .gas_joint_fifth(term, psi, g, seed, blocks_data, directions, eta = eta)
}

#' @rdname term_fifth.GasJoint
#' @param sub Whether the members carry subformulas, whose charts are then
#'   read per observation.
.gas_joint_fifth <- function(term, psi, g, seed, blocks_data, directions,
                             eta = NULL, sub = FALSE) {
  bp <- term@build_info
  K <- length(term@params)
  n <- bp$n
  seed <- lapply(seed, as.matrix)
  P <- length(seed)
  m <- ncol(seed[[1L]])
  M <- length(term_params(term))
  blocks <- blocks_data
  fp <- blocks$filtered
  g <- matrix(as.numeric(g), n, K)
  dirs <- lapply(directions, function(d) {
    d <- as.numeric(d)
    if (length(d) != m) {
      stop(sprintf("each direction must have one value per unknown (%d).", m),
           call. = FALSE)
    }
    d
  })
  u_at <- matrix(blocks$u, n, K)
  zoff <- m - M
  psi <- .gas_joint_psi(term, psi)
  cf <- .gas_joint_start(term, .gas_joint_coefs(term, psi), eta)
  zero <- .gas_jet_const(0, numeric(m), m, dirs)
  # each member's chart at a row: the constant jets of a member without
  # subformulas, the per-row jets of one with them, nothing for a driver
  charts <- lapply(seq_len(K), function(k) {
    tm <- term@terms[[term@params[k]]]
    if (.gas_is_driver(tm)) return(function(row) NULL)
    cols <- zoff + term@offsets[k] + seq_along(term_params(tm))
    if (sub) {
      return(.gas_sub_chart_fun(tm, unlist(psi[[term@params[k]]][term_params(tm)]),
                                cols, m, dirs))
    }
    base <- .gas_base_params(tm@p, tm@q)
    links <- term_links(tm)
    v <- unlist(psi[[term@params[k]]][base])
    zeta <- vapply(base, function(j)
      linkfunctions7::linkfun(links[[j]], v[[j]]), numeric(1))
    ch <- .gas_chart_jets(zeta, tm@p, tm@q, links, cols[seq_along(base)], m,
                          dirs)
    function(row) ch
  })
  # the cross loadings of each member as pairs (partner, lag): the layout's
  # order with subformulas, the coefficients' order without
  xpairs <- lapply(seq_len(K), function(k) {
    tm <- term@terms[[term@params[k]]]
    if (sub && !.gas_is_driver(tm)) {
      xp <- .gas_cross_pairs(tm)
      list(par = as.integer(match(xp$partner, term@params)),
           lag = as.integer(xp$lag))
    } else {
      list(par = cf[[k]]$xpar, lag = cf[[k]]$xlag)
    }
  })
  pairs <- list(c(1L, 2L), c(1L, 3L), c(2L, 3L))
  nm0 <- function() lapply(seq_len(K), function(k) matrix(0, n, m))
  jac <- nm0()
  dom <- nm0()
  dphi <- lapply(seq_len(K), function(k) lapply(1:3, function(d) matrix(0, n, m)))
  dpsi <- lapply(seq_len(K), function(k) lapply(1:3, function(d) matrix(0, n, m)))
  W <- matrix(0, m, m)
  arrs <- list(blocks$U1, blocks$U2, blocks$U3, blocks$U4, blocks$U5)
  for (gi in seq_along(bp$order)) {
    rows <- bp$order[[gi]]
    T_ <- length(rows)
    fj <- lapply(seq_len(K), function(k) vector("list", T_))
    uj <- lapply(seq_len(K), function(k) vector("list", T_))
    f0 <- lapply(seq_len(K), function(k) {
      if (isTRUE(cf[[k]]$direct)) {
        return(.gas_jet_const(cf[[k]]$d0[[gi]],
                              -as.numeric(seed[[fp[k]]][rows[1L], ]), m, dirs))
      }
      ch <- charts[[k]](rows[1L])
      if (is.null(ch)) zero else .gas_start_jet(ch, m, dirs)
    })
    for (t in seq_len(T_)) {
      row <- rows[t]
      for (k in seq_len(K)) {
        c_k <- cf[[k]]
        if (isTRUE(c_k$direct) && t == 1L) {
          fj[[k]][[t]] <- f0[[k]]
          next
        }
        ch <- charts[[k]](row)
        J <- if (is.null(ch)) zero else ch$omega
        for (i in seq_len(c_k$p)) {
          if (t - i < 1L) next
          J <- .gas_jet_add(J, .gas_jet_times(ch$a[[i]], uj[[k]][[t - i]]))
        }
        xs <- xpairs[[k]]
        for (h in seq_along(xs$par)) {
          if (t - xs$lag[h] < 1L) next
          if (sub) {
            cj <- ch$x[[h]]
          } else {
            d1 <- numeric(m)
            d1[zoff + c_k$col_cross[h]] <- 1
            cj <- .gas_jet_const(c_k$cross[h], d1, m, dirs)
          }
          J <- .gas_jet_add(J, .gas_jet_times(cj, uj[[xs$par[h]]][[t - xs$lag[h]]]))
        }
        for (j in seq_len(c_k$q)) {
          lagj <- if (t - j >= 1L) fj[[k]][[t - j]] else f0[[k]]
          J <- .gas_jet_add(J, .gas_jet_times(ch$b[[j]], lagj))
        }
        fj[[k]][[t]] <- J
      }
      E <- lapply(seq_len(P), function(a) {
        .gas_jet_const(0, as.numeric(seed[[a]][row, ]), m, dirs)
      })
      for (k in seq_len(K)) {
        a <- fp[k]
        J <- fj[[k]][[t]]
        E[[a]]$d1 <- E[[a]]$d1 + J$d1
        E[[a]]$mats <- J$mats
        jac[[k]][row, ] <- E[[a]]$d1
        for (d in 1:3) dphi[[k]][[d]][row, ] <- .gas_jet_get(E[[a]], d, 1L)
        for (pp in 1:3) {
          dpsi[[k]][[pp]][row, ] <- .gas_jet_get(E[[a]], pairs[[pp]], 1L)
        }
        dom[[k]][row, ] <- .gas_jet_get(E[[a]], 1:3, 1L)
        W <- W + g[row, k] * J$mats$s123
      }
      for (l in seq_len(K)) {
        ell <- lapply(seq_len(5L), function(o) {
          x <- arrs[[o]]
          idx <- c(list(row, l), rep(list(seq_len(P)), o))
          array(do.call(`[`, c(list(x), idx, list(drop = FALSE))), rep(P, o))
        })
        uj[[l]][[t]] <- .gas_compose(ell, u_at[row, l], E)
      }
    }
  }
  W <- (W + t(W)) / 2
  names(jac) <- names(dphi) <- names(dpsi) <- names(dom) <- term@params
  list(jacobian = jac, dphi = dphi, dpsi = dpsi, domega = dom, curvature = W)
}
