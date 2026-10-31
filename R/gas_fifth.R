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
