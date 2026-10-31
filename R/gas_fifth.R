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


#' The Score's Jet at One Observation, from the Family's Derivative Arrays
#'
#' @description
#' Composes the family's derivatives of the driving score in the predictors
#' with the jets of the predictors, by Faa di Bruno's formula over the set
#' partitions of the slots: the directions of a subset and the two free
#' slots.
#'
#' @details
#' The score \eqn{s = \ell_p(e_1, \dots, e_Q)} is read at the predictors of
#' every equation. The predictor \eqn{e_p} of the filter's own equation has
#' derivatives of every order; the others are linear in the unknowns, with
#' first derivative \eqn{V_r}. A block of a partition with two or more slots
#' therefore lands on \eqn{e_p} only, and a block of one slot on any
#' equation. The coefficient of a partition with \eqn{k} blocks is
#' \eqn{\ell_{p r_1 \dots r_k}}, read from `H`, `D3`, ..., `D6`.
#'
#' @param ell A list whose element \eqn{k} is the array of
#'   \eqn{\ell_{p r_1\dots r_k}} at the observation, with \eqn{k} dimensions
#'   of `npar` each.
#' @param s0 The score's value.
#' @param ej The jet of the filter's own predictor.
#' @param V A `npar` by `m` matrix of the predictors' first derivatives (the
#'   row `ap` is replaced by the jet's).
#' @param ap The index of the filter's own parameter.
#'
#' @return The score's jet.
#'
#' @keywords internal
.gas_score_jet <- function(ell, s0, ej, V, ap) {
  dirs <- ej$dirs
  nd <- length(dirs)
  npar <- nrow(V)
  V[ap, ] <- ej$d1
  m <- ncol(V)
  mats <- list()
  for (S in .gas_subsets(nd)) {
    ns <- length(S)
    slots <- c(S, -1L, -2L)
    acc <- matrix(0, m, m)
    for (part in numericals7::set_partitions(ns + 2L)) {
      k <- length(part)
      # per block: a coefficient vector over the equations for the scalar
      # blocks, and the free part for the blocks holding a free slot
      cvec <- list()
      free <- list()
      for (bi in seq_along(part)) {
        sl <- slots[part[[bi]]]
        dset <- sl[sl > 0L]
        nfree <- sum(sl < 0L)
        if (length(sl) == 1L) {
          if (nfree == 0L) {
            cvec[[bi]] <- as.numeric(V %*% dirs[[dset]])
          } else {
            free[[length(free) + 1L]] <- list(block = bi, rows = V,
                                              a = -1L %in% sl)
          }
        } else {
          val <- .gas_jet_get(ej, dset, nfree)
          if (nfree == 0L) {
            u <- numeric(npar)
            u[ap] <- val
            cvec[[bi]] <- u
          } else if (nfree == 1L) {
            rows <- matrix(0, npar, m)
            rows[ap, ] <- val
            free[[length(free) + 1L]] <- list(block = bi, rows = rows,
                                              a = -1L %in% sl)
          } else {
            free[[length(free) + 1L]] <- list(block = bi, mat = val)
          }
        }
      }
      # contract the coefficient array over the scalar blocks
      arr <- ell[[k]]
      if (k == 1L) arr <- array(arr, npar)
      fb <- vapply(free, function(f) f$block, integer(1))
      sc <- setdiff(seq_len(k), fb)
      arr <- aperm(array(arr, rep(npar, k)), c(sc, fb))
      for (bi in sc) {
        arr <- array(arr, c(npar, length(arr) / npar))
        arr <- as.numeric(crossprod(cvec[[bi]], arr))
      }
      if (length(free) == 1L) {
        # the two free slots in one block of two or more slots: e_p only
        acc <- acc + arr[ap] * free[[1L]]$mat
      } else {
        # the first remaining index is the first free block's; the matrix
        # is laid out with the slot a in its rows
        cmat <- matrix(arr, npar, npar)
        if (free[[1L]]$a) {
          acc <- acc + crossprod(free[[1L]]$rows, cmat %*% free[[2L]]$rows)
        } else {
          acc <- acc + crossprod(free[[2L]]$rows, t(cmat) %*% free[[1L]]$rows)
        }
      }
    }
    mats[[.gas_jet_key(S)]] <- acc / 2 + t(acc) / 2
  }
  d1 <- as.numeric(crossprod(V, as.numeric(ell[[1L]])))
  .gas_jetn(s0, d1, mats, dirs)
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
#' through the recursion, for a term without subformulas, and returns the
#' contractions of [term_fifth()].
#'
#' @details
#' At each observation the level is
#' \eqn{f_t = \omega + \sum_i a_i s_{t-i} + \sum_j b_j f_{t-j}}, and its jet
#' is the sum of the products of jets ([.gas_jet_times()]). The score's jet
#' is composed from the family's derivative arrays by [.gas_score_jet()]. The
#' orders below the fifth are carried along every subset of the directions,
#' so this route is a second, independent computation of what
#' [term_curvature()], [term_third()] and [term_fourth()] return.
#'
#' @inheritParams term_fifth
#' @param dirs The three directions.
#'
#' @return As [term_fifth()].
#'
#' @keywords internal
.gas_fifth_scalar <- function(term, eta, y, score, psi, g, seed, blocks_data,
                              dirs) {
  bp <- term@build_info
  nm <- term_params(term)
  links <- term_links(term)
  psiv <- unlist(psi[nm])
  base <- .gas_base_params(term@p, term@q)
  nb <- length(base)
  zeta <- vapply(base, function(j)
    linkfunctions7::linkfun(links[[j]], psiv[[j]]), numeric(1))
  seed <- as.matrix(seed)
  m <- ncol(seed)
  np <- length(nm)
  n <- bp$n
  zcol <- m - np + seq_len(np)
  p <- term@p
  q <- term@q
  ch <- .gas_chart_jets(zeta, p, q, links, zcol[seq_len(nb)], m, dirs)
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
    f0j <- ch$f0
    if (static) f0j <- .gas_jet_const(d0[[l]], -seed[rows[1L], ], m, dirs)
    fj <- vector("list", k)
    sj <- vector("list", k)
    for (t in seq_len(k)) {
      row <- rows[t]
      if (static && t == 1L) {
        ft <- f0j
      } else {
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
