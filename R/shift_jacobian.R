#' @include term_classes.R generics.R structural.R regime.R marginal.R
NULL

#' @title The Posterior Shift of a Term Mixed Over States, and Its Jacobian
#'
#' @description
#' Computes, for a structural term whose contribution is a likelihood mixed
#' over latent states, the shift that the term adds to its equation's
#' predictor at each observation, averaged over the smoothed posterior of the
#' states, together with its exact derivative in the whole of a caller's
#' unknown vector.
#'
#' @details
#' A term of the likelihood shape reports no predictor of its own. The
#' predictor of its equation at observation \eqn{i} is the static part plus
#' the posterior mean of the shift,
#' \deqn{\bar s_i = \sum_c \gamma_{ic}\, s_{ic},}
#' with \eqn{\gamma_{ic}} the posterior probability of component \eqn{c}
#' given the whole series ([term_posterior()]) and \eqn{s_{ic}} its shift
#' ([term_levels()]). The derivative in the unknown vector \eqn{u} is
#' \deqn{\frac{\partial \bar s_i}{\partial u} = \sum_c
#'   \frac{\partial \gamma_{ic}}{\partial u}\, s_{ic}
#'   + \gamma_{ic}\, \frac{\partial s_{ic}}{\partial u},}
#' and the posterior moves with every unknown: with the coefficients of
#' every equation, through the emissions, and with the term's own
#' parameters, through the emissions, the shifts and the prior of the
#' states. The delta method for a prediction reads this derivative.
#'
#' Where the posterior is the smoothed marginal of a hidden Markov chain (a
#' regime chain, and the side chain of a marginal step), the derivatives
#' are carried by forward mode through the normalized forward and backward
#' recursions and combined by the quotient rule; the scale of each step
#' multiplies a whole vector and cancels in the normalization of the
#' posterior, so its own derivative is not needed. Where the posterior is
#' a group's distribution over quadrature nodes (the continuous marginal
#' kinds), \eqn{\partial \gamma_c = \gamma_c (\partial a_c - \sum_{c'}
#' \gamma_{c'} \partial a_{c'})} with \eqn{a_c} the log of the node's
#' weight times the group's emissions.
#'
#' The derivatives in the term's own parameters are taken on the
#' unconstrained scale of [term_links()], as [term_hessian()] takes them.
#'
#' @param term A built term of the likelihood shape.
#' @param eta The static part of the predictor that the term shifts.
#' @param y The response.
#' @param logdens A function of a predictor value and a row index returning
#'   the log-density there, as [term_loglik()] takes it.
#' @param grad A function of the same two arguments returning a matrix with
#'   one row per observation and one column per distribution parameter, the
#'   derivative of the log-density in each predictor.
#' @param psi The term's parameters on the parameter scale, named as
#'   [term_params()].
#' @param seed A list with one matrix per distribution parameter, each with
#'   one row per observation and one column per unknown, giving the
#'   derivative of that predictor in the caller's unknowns, with zeros in the
#'   columns that the term's own parameters occupy.
#' @param cols The columns of the unknown vector that the term's own
#'   parameters occupy, in the order of [term_params()].
#' @param level The index, among the distribution parameters, of the one that
#'   the term shifts.
#' @param ... Passed to methods.
#'
#' @return A list with `shift`, the `n` posterior shifts, and `jacobian`, an
#'   `n` by `m` matrix with `m` the length of the unknown vector.
#'
#' @examples
#' set.seed(1)
#' dd <- data.frame(t = 1:30, y = c(rnorm(15), rnorm(15, 3)))
#' term <- term_build(regime(2, time = t), dd)
#' psi <- list(level1 = 0, gap2 = 3, alr1.1 = 2, alr2.1 = -2)
#' out <- term_shift_jacobian(
#'   term, eta = rep(0, 30), y = dd$y,
#'   logdens = function(e, i) dnorm(dd$y[i], e, log = TRUE),
#'   grad = function(e, i) matrix(dd$y[i] - e, ncol = 1L),
#'   psi = psi, seed = list(matrix(c(rep(1, 30), rep(0, 30 * 4)), 30, 5)),
#'   cols = 2:5, level = 1L)
#' dim(out$jacobian)
#'
#' @seealso [term_posterior()], [term_levels()], [term_hessian()]
#' @export
#' @aliases term_shift_jacobian.structural_term
term_shift_jacobian <- S7::new_generic("term_shift_jacobian", "term",
  function(term, eta, y, logdens, grad, psi, seed, cols, level, ...)
    S7::S7_dispatch())

S7::method(term_shift_jacobian, structural_term) <- function(term, eta, y,
                                                             logdens, grad,
                                                             psi, seed, cols,
                                                             level, ...) {
  stop(sprintf("the term class '%s' does not implement term_shift_jacobian().",
               attr(S7::S7_class(term), "name")), call. = FALSE)
}

# the checks every method makes of the caller's seed, columns and level
.shift_inputs <- function(seed, n, cols, level, nm) {
  if (!is.list(seed) || !length(seed)) {
    stop("'seed' must be a list with one matrix per distribution parameter.",
         call. = FALSE)
  }
  seed <- lapply(seed, as.matrix)
  m <- ncol(seed[[1L]])
  if (any(vapply(seed, nrow, 1L) != n) ||
      any(vapply(seed, ncol, 1L) != m)) {
    stop(sprintf("every element of 'seed' must be %d by %d.", n, m),
         call. = FALSE)
  }
  cols <- as.integer(cols)
  level <- as.integer(level)
  if (length(cols) != length(nm) || any(cols < 1L) || any(cols > m)) {
    stop(sprintf("'cols' must give the %d columns that the term's parameters occupy.",
                 length(nm)), call. = FALSE)
  }
  if (length(level) != 1L || level < 1L || level > length(seed)) {
    stop("'level' must index one of the distribution parameters.",
         call. = FALSE)
  }
  list(seed = seed, m = m, cols = cols, level = level)
}

# the emission score of one component in the unknowns, one row per
# observation: the static rows weighted by the family's score in each
# predictor, without the component's own shift
.shift_emission <- function(g, seed) {
  gm <- 0 * seed[[1L]]
  for (q in seq_along(seed)) gm <- gm + g[, q] * seed[[q]]
  gm
}

#' @title The Posterior Shift of a Regime Term, and Its Jacobian
#' @name term_shift_jacobian.RegimeTerm
#' @description
#' Computes the shift and its Jacobian by forward mode through the normalized
#' forward and backward recursions of [term_posterior()], over each group in
#' time order.
#' @param term A built [RegimeTerm()].
#' @param eta,y,logdens,grad,psi,seed,cols,level As for
#'   [term_shift_jacobian()].
#' @param ... Unused.
#' @return A list with `shift` and `jacobian`.
#' @keywords internal
S7::method(term_shift_jacobian, RegimeTerm) <- function(term, eta, y, logdens,
                                                        grad, psi, seed, cols,
                                                        level, ...) {
  bp <- term@blueprint
  if (!length(bp)) {
    stop("the term has not been built; call term_build(term, data) first.",
         call. = FALSE)
  }
  nm <- term_params(term)
  v <- unlist(psi[nm])
  if (length(v) != length(nm) || anyNA(v)) {
    stop(sprintf("'psi' must supply %s.", paste(nm, collapse = ", ")),
         call. = FALSE)
  }
  k <- term@k
  n <- bp$n
  ins <- .shift_inputs(seed, n, cols, level, nm)
  seed <- ins$seed
  m <- ins$m
  cols <- ins$cols
  level <- ins$level
  npar <- length(seed)

  ch <- .regime_chart(v, k, nm, term@chain)
  np <- ch$np
  P <- ch$P
  # the level of each regime and its derivative in the unknowns
  dmu <- matrix(0, k, m)
  dmu[, cols] <- ch$dmu
  idx <- seq_len(n)
  LF <- matrix(0, n, k)
  GM <- vector("list", k)
  for (j in seq_len(k)) {
    e <- eta + ch$mu[j]
    a <- as.numeric(logdens(e, idx))
    if (length(a) != n) {
      stop(paste0("'logdens' must return one value per observation when ",
                  "given the whole index vector."), call. = FALSE)
    }
    LF[, j] <- a
    g <- as.matrix(grad(e, idx))
    if (nrow(g) != n || ncol(g) != npar) {
      stop(sprintf("'grad' must return an %d by %d matrix.", n, npar),
           call. = FALSE)
    }
    gm <- .shift_emission(g, seed)
    gm[, cols] <- gm[, cols] + outer(g[, level], ch$dmu[j, ])
    GM[[j]] <- gm
  }
  # the chain's own derivative, column l of P in the term's coordinates
  Dl <- lapply(seq_len(k), function(l) matrix(ch$dPz[, l, ], k, np))

  shift <- numeric(n)
  jac <- matrix(0, n, m)
  for (rows in bp$order) {
    mr <- length(rows)
    W <- exp(LF[rows, , drop = FALSE] - apply(LF[rows, , drop = FALSE], 1L, max))
    # forward, normalized at every step
    A <- matrix(0, mr, k)
    dA <- array(0, c(mr, k, m))
    a <- ch$delta
    da <- matrix(0, k, m)
    da[, cols] <- t(do.call(rbind, ch$ddelta))
    for (t in seq_len(mr)) {
      if (t == 1L) {
        pred <- a
        dpred <- da
      } else {
        pred <- as.numeric(a %*% P)
        dpred <- crossprod(P, da)
        for (l in seq_len(k)) {
          dpred[l, cols] <- dpred[l, cols] + as.numeric(a %*% Dl[[l]])
        }
      }
      atil <- W[t, ] * pred
      datil <- W[t, ] * (dpred + pred * do.call(rbind, lapply(seq_len(k),
        function(l) GM[[l]][rows[t], ])))
      cs <- sum(atil)
      dcs <- colSums(datil)
      a <- atil / cs
      da <- (datil - outer(a, dcs)) / cs
      A[t, ] <- a
      dA[t, , ] <- da
    }
    # backward, normalized the same way
    B <- matrix(1, mr, k)
    dB <- array(0, c(mr, k, m))
    b <- rep(1, k)
    db <- matrix(0, k, m)
    if (mr > 1L) {
      for (t in seq.int(mr - 1L, 1L)) {
        r1 <- rows[t + 1L]
        u1 <- W[t + 1L, ] * b
        du1 <- W[t + 1L, ] * (db + b * do.call(rbind, lapply(seq_len(k),
          function(l) GM[[l]][r1, ])))
        num <- as.numeric(P %*% u1)
        dnum <- P %*% du1
        for (l in seq_len(k)) {
          dnum[, cols] <- dnum[, cols] + u1[l] * Dl[[l]]
        }
        s <- sum(num)
        if (s > 0) {
          b <- num / s
          db <- (dnum - outer(b, colSums(dnum))) / s
        } else {
          b <- num
          db <- dnum
        }
        B[t, ] <- b
        dB[t, , ] <- db
      }
    }
    # the smoothed posterior and the shift it weights
    for (t in seq_len(mr)) {
      g <- A[t, ] * B[t, ]
      dg <- dA[t, , ] * B[t, ] + A[t, ] * dB[t, , ]
      dg <- matrix(dg, k, m)
      sg <- sum(g)
      gam <- g / sg
      dgam <- (dg - outer(gam, colSums(dg))) / sg
      shift[rows[t]] <- sum(gam * ch$mu)
      jac[rows[t], ] <- as.numeric(crossprod(dgam, ch$mu)) +
        as.numeric(crossprod(gam, dmu))
    }
  }
  list(shift = shift, jacobian = jac)
}

#' @title The Posterior Shift of a Marginal Break-Point Term, and Its Jacobian
#' @name term_shift_jacobian.MarginalBreakTerm
#' @description
#' Computes the shift and its Jacobian. For the step kind, the function applies
#' forward mode through the forward and backward recursions of the side chain
#' over the \eqn{2^K} side patterns, and the interval masses and the survivals
#' carry their derivatives in the prior's parameters. For the continuous kinds,
#' it computes the derivative of each group's posterior over its quadrature
#' nodes, and the node weights, the node positions and the shifts carry their
#' derivatives in the term's own parameters.
#' @param term A built [MarginalBreakTerm()].
#' @param eta,y,logdens,grad,psi,seed,cols,level As for
#'   [term_shift_jacobian()].
#' @param ... Unused.
#' @return A list with `shift` and `jacobian`.
#' @keywords internal
S7::method(term_shift_jacobian, MarginalBreakTerm) <- function(term, eta, y,
                                                               logdens, grad,
                                                               psi, seed, cols,
                                                               level, ...) {
  bp <- .marg_built(term)
  v <- .marg_check_psi(term, psi)
  nm <- term_params(term)
  n <- bp$n
  ins <- .shift_inputs(seed, n, cols, level, nm)
  # the chain onto the unconstrained scale, one factor per own parameter
  links <- term_links(term)
  J1 <- vapply(nm, function(j)
    linkfunctions7::dlinkinv(links[[j]],
                             linkfunctions7::linkfun(links[[j]], v[[j]])),
    numeric(1))
  if (term@kind == "jump") {
    .marg_jump_shift_jacobian(term, eta, logdens, grad, v, ins, nm, J1)
  } else {
    .marg_seg_shift_jacobian(term, eta, logdens, grad, v, ins, nm, J1)
  }
}

.marg_seg_shift_jacobian <- function(term, eta, logdens, grad, v, ins, nm,
                                     J1) {
  bp <- term@blueprint
  n <- bp$n
  seed <- ins$seed
  m <- ins$m
  cols <- ins$cols
  level <- ins$level
  npar <- length(seed)
  col_of <- function(j) cols[match(j, nm)]
  own <- c(if (term@linear) "beta", "gamma1",
           if (term@kind == "jseg") "delta1")
  cm <- col_of(.marg_mname(1L))
  ct <- col_of(.marg_tname(1L))
  shift <- numeric(n)
  jac <- matrix(0, n, m)
  for (rs in bp$groups) {
    ng <- length(rs)
    nd <- .marg_nodes(bp$x[rs], v[[.marg_mname(1L)]], v[[.marg_tname(1L)]])
    C <- length(nd$p)
    sh <- .marg_seg_shift(term, bp$x[rs], nd, v)
    ei <- rep(eta[rs], C) + as.numeric(sh$shift)
    ii <- rep(rs, C)
    LD <- matrix(as.numeric(logdens(ei, ii)), ng, C)
    G <- as.matrix(grad(ei, ii))
    if (nrow(G) != ng * C || ncol(G) != npar) {
      stop(sprintf("'grad' must return a matrix of %d columns, one row per value.",
                   npar), call. = FALSE)
    }
    # the derivative of each node's shift at each observation, in the term's
    # own columns, on the unconstrained scale: the coefficients through the
    # covariate, the hinge and the step, the prior's location and scale
    # through the motion of the node
    dsh <- array(0, c(ng, C, m))
    if (term@linear) {
      dsh[, , col_of("beta")] <- matrix(bp$x[rs], ng, C) * J1[["beta"]]
    }
    dsh[, , col_of("gamma1")] <- sh$hinge * J1[["gamma1"]]
    if (term@kind == "jseg") dsh[, , col_of("delta1")] <- sh$step * J1[["delta1"]]
    dsh[, , cm] <- sweep(sh$dshift_dpsi, 2L, nd$dpsi_m, `*`) *
      J1[[.marg_mname(1L)]]
    dsh[, , ct] <- sweep(sh$dshift_dpsi, 2L, nd$dpsi_t, `*`) *
      J1[[.marg_tname(1L)]]
    # the posterior over the nodes and the derivative of its log weights
    la <- nd$lw + colSums(LD)
    gam <- exp(la - .marg_lse(la))
    dla <- matrix(0, C, m)
    dla[, cm] <- nd$glw_m * J1[[.marg_mname(1L)]]
    dla[, ct] <- nd$glw_t * J1[[.marg_tname(1L)]]
    for (c in seq_len(C)) {
      rr <- (c - 1L) * ng + seq_len(ng)
      gm <- matrix(0, ng, m)
      for (q in seq_len(npar)) gm <- gm + G[rr, q] * seed[[q]][rs, , drop = FALSE]
      gm <- gm + G[rr, level] * matrix(dsh[, c, ], ng, m)
      dla[c, ] <- dla[c, ] + colSums(gm)
    }
    dgam <- gam * sweep(dla, 2L, colSums(gam * dla), `-`)
    shift[rs] <- as.numeric(sh$shift %*% gam)
    for (t in seq_len(ng)) {
      jac[rs[t], ] <- as.numeric(crossprod(dgam, sh$shift[t, ])) +
        as.numeric(crossprod(gam, matrix(dsh[t, , ], C, m)))
    }
  }
  list(shift = shift, jacobian = jac)
}

.marg_jump_shift_jacobian <- function(term, eta, logdens, grad, v, ins, nm,
                                      J1) {
  bp <- term@blueprint
  n <- bp$n
  seed <- ins$seed
  m <- ins$m
  cols <- ins$cols
  level <- ins$level
  npar <- length(seed)
  K <- term@npsi
  P <- 2^K
  idx <- seq_len(n)
  pb <- .marg_bits(K, v[paste0("delta", seq_len(K))])
  dcol <- vapply(seq_len(K), function(k)
    cols[match(paste0("delta", k), nm)], integer(1))
  LD <- matrix(0, n, P)
  GM <- vector("list", P)
  for (p in seq_len(P)) {
    e <- eta + pb$shifts[p]
    LD[, p] <- as.numeric(logdens(e, idx))
    g <- as.matrix(grad(e, idx))
    if (nrow(g) != n || ncol(g) != npar) {
      stop(sprintf("'grad' must return an %d by %d matrix.", n, npar),
           call. = FALSE)
    }
    gm <- .shift_emission(g, seed)
    for (k in seq_len(K)) {
      if (pb$bits[p, k]) {
        gm[, dcol[k]] <- gm[, dcol[k]] + g[, level] * J1[[paste0("delta", k)]]
      }
    }
    GM[[p]] <- gm
  }
  # the shift of each pattern and its derivative in the deltas' columns
  dshift <- matrix(0, P, m)
  for (k in seq_len(K)) {
    dshift[, dcol[k]] <- pb$bits[, k] * J1[[paste0("delta", k)]]
  }
  tix <- .marg_hmm_idx(K)
  shift <- numeric(n)
  jac <- matrix(0, n, m)
  for (rs in bp$groups) {
    ng <- length(rs)
    hm <- .marg_hmm_masses(term, bp$x[rs], v)
    pcol <- lapply(hm$pcols, function(pc) cols[match(pc, nm)])
    # the masses' and the survivals' derivatives onto the unconstrained scale
    dM <- lapply(seq_len(K), function(k) {
      sweep(hm$dM[[k]], 2L, J1[match(hm$pcols[[k]], nm)], `*`)
    })
    dSV <- lapply(seq_len(K), function(k) {
      sweep(hm$dSV[[k]], 2L, J1[match(hm$pcols[[k]], nm)], `*`)
    })
    em <- function(t) {
      lf <- LD[rs[t], ]
      w <- exp(lf - max(lf))
      list(w = w, g = do.call(rbind, lapply(seq_len(P),
                                            function(p) GM[[p]][rs[t], ])))
    }
    # forward: the flips of step t, coordinate by coordinate as in
    # .marg_hmm_forward(), then the emission
    A <- matrix(0, ng, P)
    dA <- array(0, c(ng, P, m))
    alpha <- numeric(P)
    alpha[1L] <- 1
    dal <- matrix(0, P, m)
    for (t in seq_len(ng)) {
      for (k in seq_len(K)) {
        i0 <- tix[[k]]$i0
        i1 <- tix[[k]]$i1
        q <- hm$M[t, k]
        dal[i1, ] <- dal[i1, ] + q * dal[i0, , drop = FALSE]
        dal[i1, pcol[[k]]] <- dal[i1, pcol[[k]]] + alpha[i0] %o% dM[[k]][t, ]
        alpha[i1] <- alpha[i1] + q * alpha[i0]
      }
      e <- em(t)
      dal <- e$w * (dal + alpha * e$g)
      alpha <- alpha * e$w
      cs <- max(alpha)
      if (cs > 0) {
        alpha <- alpha / cs
        dal <- dal / cs
      }
      A[t, ] <- alpha
      dA[t, , ] <- dal
    }
    # backward: the survival product of the final states, then each step's
    # emission and flips in the reverse order of .marg_hmm_backward()
    B <- matrix(0, ng, P)
    dB <- array(0, c(ng, P, m))
    SVf <- hm$SV[ng + 1L, ]
    beta <- rep(1, P)
    dbe <- matrix(0, P, m)
    for (k in seq_len(K)) {
      i0 <- tix[[k]]$i0
      dbe[i0, ] <- SVf[k] * dbe[i0, , drop = FALSE]
      dbe[i0, pcol[[k]]] <- dbe[i0, pcol[[k]]] +
        beta[i0] %o% dSV[[k]][ng + 1L, ]
      beta[i0] <- beta[i0] * SVf[k]
    }
    B[ng, ] <- beta
    dB[ng, , ] <- dbe
    if (ng > 1L) {
      for (t in seq.int(ng, 2L)) {
        e <- em(t)
        dtmp <- e$w * (dbe + beta * e$g)
        tmp <- beta * e$w
        for (k in seq_len(K)) {
          i0 <- tix[[k]]$i0
          i1 <- tix[[k]]$i1
          q <- hm$M[t, k]
          dtmp[i0, ] <- dtmp[i0, ] + q * dtmp[i1, , drop = FALSE]
          dtmp[i0, pcol[[k]]] <- dtmp[i0, pcol[[k]]] + tmp[i1] %o% dM[[k]][t, ]
          tmp[i0] <- tmp[i0] + q * tmp[i1]
        }
        cs <- max(tmp)
        if (cs > 0) {
          tmp <- tmp / cs
          dtmp <- dtmp / cs
        }
        beta <- tmp
        dbe <- dtmp
        B[t - 1L, ] <- beta
        dB[t - 1L, , ] <- dbe
      }
    }
    for (t in seq_len(ng)) {
      g <- A[t, ] * B[t, ]
      dg <- matrix(dA[t, , ] * B[t, ] + A[t, ] * dB[t, , ], P, m)
      sg <- sum(g)
      gam <- g / sg
      dgam <- (dg - outer(gam, colSums(dg))) / sg
      shift[rs[t]] <- sum(gam * pb$shifts)
      jac[rs[t], ] <- as.numeric(crossprod(dgam, pb$shifts)) +
        as.numeric(crossprod(gam, dshift))
    }
  }
  list(shift = shift, jacobian = jac)
}
