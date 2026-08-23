#' @include gas.R
NULL

# The submodel route of a score-driven term: every quantity of the filter
# read per observation. A parameter carrying a subformula is
# psi_{j,t} = g_j^-1(z_t' gamma_j), so the recursion's coefficients vary
# with t and the scalar kernel does not apply; the functions here run the
# general recursion -- the filter, its exact derivative, the reverse pass
# and the second-order propagation -- in R, sharing their shape with
# gas_filter_r() and the scalar term_curvature() method.
#
# Everything is carried on each group's ACTIVE SET: a scalar coordinate
# reaches every group, and a development's coordinate reaches the groups
# where its column is not identically zero, so with grouping indicators
# (a random intercept per panel unit) the per-observation work is constant
# in the number of groups. Forming the full square instead was measured at
# about twelve minutes at five hundred groups, which is the trap the
# deviations machinery this route replaced had already closed.

# the names and positions of the term's parameters with the developed ones
# expanded in place
.gas_sub_layout <- function(term) {
  base <- .gas_base_params(term@p, term@q)
  sub <- term@blueprint$sub
  nm <- character(0)
  idx <- list()
  pos <- 0L
  for (j in base) {
    if (!is.null(sub[[j]])) {
      k <- length(sub[[j]]$coef_names)
      nm <- c(nm, paste(j, sub[[j]]$coef_names, sep = "."))
      idx[[j]] <- pos + seq_len(k)
      pos <- pos + k
    } else {
      nm <- c(nm, j)
      idx[[j]] <- pos + 1L
      pos <- pos + 1L
    }
  }
  list(names = nm, idx = idx)
}

S7::method(term_components, GasTerm) <- function(term, ...) {
  bp <- term@blueprint
  if (!length(bp)) return(list())
  # a structural term contributes no design columns, so what divides here is
  # the PARAMETER vector: `index` gives positions in term_params(), which is
  # the vector its state, its readable quantities and its variance matrix are
  # all indexed by. For an additive term the same field gives columns of the
  # block, and in both cases it is the term's own coefficients.
  lay <- .gas_sub_layout(term)
  sub <- bp$sub
  out <- lapply(names(lay$idx), function(j) {
    s <- if (is.null(sub[[j]])) list() else sub[[j]]$terms
    list(name = j, index = lay$idx[[j]], subs = s,
         sub_index = component_sub_index(lay$idx[[j]], s))
  })
  stats::setNames(out, names(lay$idx))
}

# the link a parameter's own chart rides, override or default; the
# COORDINATES of a developed parameter are unconstrained, and this is the
# link applied inside the development
.gas_param_link <- function(term, j) {
  if (!is.null(term@links[[j]])) return(term@links[[j]])
  if (startsWith(j, "pacf")) linkfunctions7::rhobit_link()
  else if (startsWith(j, "alpha")) linkfunctions7::log_link()
  else linkfunctions7::identity_link()
}

# The per-observation value of every base parameter and its derivatives in
# the term's parameters. A developed parameter reads its chart at
# z_t' gamma; a scalar one is constant, and the scale of its coordinate is
# the caller's -- the parameter scale for the filter, whose jacobian
# statmodels7 chains through the links outside, and the unconstrained
# scale for the curvature, which differentiates the chart itself. W holds
# d value / d coordinate row by row over the parameter's own coordinates,
# and k2 the diagonal chart curvature the second order needs.
.gas_sub_values <- function(term, u, scale = c("parameter", "zeta"),
                            second = FALSE, third = FALSE) {
  scale <- match.arg(scale)
  base <- .gas_base_params(term@p, term@q)
  sub <- term@blueprint$sub
  lay <- .gas_sub_layout(term)
  n <- term@blueprint$n
  out <- list()
  for (j in base) {
    idx <- lay$idx[[j]]
    if (!is.null(sub[[j]])) {
      lk <- .gas_param_link(term, j)
      Z <- as.matrix(sub[[j]]$Z)
      eta <- as.numeric(Z %*% u[idx])
      out[[j]] <- list(
        idx = idx, developed = TRUE, Z = Z,
        v = linkfunctions7::linkinv(lk, eta),
        W = linkfunctions7::dlinkinv(lk, eta) * Z,
        k2 = if (second) linkfunctions7::d2linkinv(lk, eta) else NULL,
        k3 = if (third) linkfunctions7::d3linkinv(lk, eta) else NULL)
    } else if (identical(scale, "parameter")) {
      out[[j]] <- list(idx = idx, developed = FALSE, Z = NULL,
                       v = rep(u[[idx]], n), W = matrix(1, n, 1L),
                       k2 = if (second) rep(0, n) else NULL,
                       k3 = if (third) rep(0, n) else NULL)
    } else {
      lk <- .gas_param_link(term, j)
      out[[j]] <- list(
        idx = idx, developed = FALSE, Z = NULL,
        v = rep(linkfunctions7::linkinv(lk, u[[idx]]), n),
        W = matrix(linkfunctions7::dlinkinv(lk, u[[idx]]), n, 1L),
        k2 = if (second) rep(linkfunctions7::d2linkinv(lk, u[[idx]]), n)
          else NULL,
        k3 = if (third) rep(linkfunctions7::d3linkinv(lk, u[[idx]]), n)
          else NULL)
    }
  }
  out
}

# each group's active coordinates among the term's parameters: a scalar
# coordinate reaches every group, and a development's coordinate reaches
# the groups where its column is not identically zero
.gas_sub_active <- function(term, vals) {
  base <- .gas_base_params(term@p, term@q)
  lapply(term@blueprint$order, function(rows) {
    act <- integer(0)
    for (j in base) {
      vj <- vals[[j]]
      if (!vj$developed) {
        act <- c(act, vj$idx)
      } else {
        nz <- colSums(abs(vj$Z[rows, , drop = FALSE])) > 0
        act <- c(act, vj$idx[nz])
      }
    }
    sort(act)
  })
}

# The autoregressive coefficients per observation and their derivative
# rows over the partial autocorrelations' own coordinates: the
# Levinson-Durbin map at that observation's values, chained onto whatever
# coordinates each one has. Constant, and computed once, when no partial
# autocorrelation is developed. DB spans only the pacf coordinates
# (pidx), every other entry being exactly zero.
.gas_sub_b <- function(term, vals) {
  q <- term@q
  n <- term@blueprint$n
  if (q == 0L) return(list(B = NULL, DB = NULL, pidx = integer(0)))
  pj <- paste0("pacf", seq_len(q))
  pidx <- sort(unlist(lapply(pj, function(j) vals[[j]]$idx)))
  npb <- length(pidx)
  B <- matrix(0, n, q)
  DB <- replicate(q, matrix(0, n, npb), simplify = FALSE)
  varying <- any(vapply(pj, function(j) vals[[j]]$developed, logical(1)))
  if (!varying) {
    ld <- gas_levinson(vapply(pj, function(j) vals[[j]]$v[1L], numeric(1)))
    for (j in seq_len(q)) {
      B[, j] <- ld$phi[j]
      row <- numeric(npb)
      for (k in seq_len(q)) {
        vk <- vals[[pj[k]]]
        at <- match(vk$idx, pidx)
        row[at] <- row[at] + ld$jacobian[j, k] * vk$W[1L, ]
      }
      DB[[j]] <- matrix(row, n, npb, byrow = TRUE)
    }
  } else {
    for (r in seq_len(n)) {
      ld <- gas_levinson(vapply(pj, function(j) vals[[j]]$v[r], numeric(1)))
      for (j in seq_len(q)) {
        B[r, j] <- ld$phi[j]
        for (k in seq_len(q)) {
          vk <- vals[[pj[k]]]
          at <- match(vk$idx, pidx)
          DB[[j]][r, at] <- DB[[j]][r, at] + ld$jacobian[j, k] * vk$W[r, ]
        }
      }
    }
  }
  list(B = B, DB = DB, pidx = pidx)
}

# The per-group preparation the compiled route reads: the derivative row
# of every coefficient scattered onto the group's active columns once,
# vectorized (one matrix assignment per parameter per group), so the
# kernel's loop allocates nothing per step. The starting level and its
# derivative come from the group's first observation, as in the R route.
.gas_sub_prep <- function(term, vals, bd, acts) {
  bp <- term@blueprint
  p <- term@p
  q <- term@q
  aj <- if (p > 0L) paste0("alpha", seq_len(p)) else character(0)
  scatter <- function(j, rows, act) {
    vj <- vals[[j]]
    at <- match(vj$idx, act)
    ok <- !is.na(at)
    M <- matrix(0, length(rows), length(act))
    if (any(ok)) M[, at[ok]] <- vj$W[rows, ok, drop = FALSE]
    M
  }
  groups <- vector("list", length(bp$order))
  for (l in seq_along(bp$order)) {
    rows <- bp$order[[l]]
    act <- acts[[l]]
    dom <- scatter("omega", rows, act)
    da <- lapply(aj, scatter, rows = rows, act = act)
    dbl <- list()
    if (q > 0L) {
      pat <- match(bd$pidx, act)
      ok <- !is.na(pat)
      dbl <- lapply(seq_len(q), function(j) {
        M <- matrix(0, length(rows), length(act))
        if (any(ok)) M[, pat[ok]] <- bd$DB[[j]][rows, ok, drop = FALSE]
        M
      })
    }
    r1 <- rows[1L]
    sb <- if (q > 0L) sum(bd$B[r1, ]) else 0
    if (abs(1 - sb) < 1e-10) {
      stop("the autoregressive polynomial is at the unit root; the filter has no starting level.",
           call. = FALSE)
    }
    om1 <- vals$omega$v[r1]
    f0 <- om1 / (1 - sb)
    df0 <- dom[1L, ] / (1 - sb)
    if (q > 0L) {
      dbsum <- Reduce(`+`, lapply(dbl, function(M) M[1L, ]))
      df0 <- df0 + om1 * dbsum / (1 - sb)^2
    }
    groups[[l]] <- list(rows = rows, act = act, dom = dom, da = da,
                        dbl = dbl, f0 = f0, df0 = df0)
  }
  A <- matrix(0, bp$n, p)
  for (i in seq_len(p)) A[, i] <- vals[[aj[i]]]$v
  B <- matrix(0, bp$n, q)
  for (j in seq_len(q)) B[, j] <- bd$B[, j]
  list(groups = groups, om = vals$omega$v, A = A, B = B)
}

# The per-group preparation the compiled CURVATURE route reads
# (gas_curvature_sub_cpp): every per-observation quantity of the
# second-order recursion scattered once onto the group's active columns --
# the chained rows W, the raw chart rows Z with their curvature k2 (a
# scalar coordinate scatters a one, so k2 * (1 * 1) is its diagonal
# entry exactly), the Levinson-chained db rows, the CONSTANT hb blocks of
# the autoregressive chart, and the starting level with its two
# derivatives. Every expression here is the R recursion's own, so the
# prep is shared arithmetic and not a second derivation.
.gas_curv_prep <- function(term, vals, bd, acts, lay, m, np, zcol,
                           ld2_const, base, aj, pj, p, q) {
  bp <- term@blueprint
  groups <- vector("list", length(bp$order))
  for (l in seq_along(bp$order)) {
    rows <- bp$order[[l]]
    act_t <- acts[[l]]
    act <- c(seq_len(m - np), zcol[act_t])
    mk <- length(act)
    pos <- lapply(base, function(j) (m - np) + match(vals[[j]]$idx, act_t))
    names(pos) <- base
    ppos <- (m - np) + match(bd$pidx, act_t)

    wrow_mat <- function(j) {
      vj <- vals[[j]]
      at <- pos[[j]]
      ok <- !is.na(at)
      M <- matrix(0, length(rows), mk)
      if (any(ok)) M[, at[ok]] <- vj$W[rows, ok, drop = FALSE]
      M
    }
    zrow_mat <- function(j) {
      vj <- vals[[j]]
      at <- pos[[j]]
      ok <- !is.na(at)
      M <- matrix(0, length(rows), mk)
      if (any(ok)) {
        M[, at[ok]] <- if (vj$developed) vj$Z[rows, ok, drop = FALSE] else 1
      }
      M
    }
    db_mat <- function(j) {
      M <- matrix(0, length(rows), mk)
      ok <- !is.na(ppos)
      if (any(ok)) M[, ppos[ok]] <- bd$DB[[j]][rows, ok, drop = FALSE]
      M
    }
    row_of <- function(j, r) {
      out <- numeric(mk)
      at <- pos[[j]]
      ok <- !is.na(at)
      out[at[ok]] <- vals[[j]]$W[r, ok]
      out
    }
    h_of <- function(j, r) {
      h <- matrix(0, mk, mk)
      vj <- vals[[j]]
      at <- pos[[j]]
      ok <- !is.na(at)
      if (!any(ok)) return(h)
      if (vj$developed) {
        z <- vj$Z[r, ok]
        h[at[ok], at[ok]] <- vj$k2[r] * outer(z, z)
      } else {
        h[at[ok], at[ok]] <- vj$k2[r]
      }
      h
    }
    db_of <- function(j, r) {
      out <- numeric(mk)
      ok <- !is.na(ppos)
      out[ppos[ok]] <- bd$DB[[j]][r, ok]
      out
    }

    r1 <- rows[1L]
    hb <- if (q > 0L) {
      wrows <- lapply(pj, function(j) row_of(j, r1))
      lapply(seq_len(q), function(j) {
        h <- matrix(0, mk, mk)
        for (k in seq_len(q)) {
          for (ll in seq_len(q)) {
            if (ld2_const$hessian[[j]][k, ll] != 0) {
              h <- h + ld2_const$hessian[[j]][k, ll] *
                outer(wrows[[k]], wrows[[ll]])
            }
          }
          if (ld2_const$jacobian[j, k] != 0) {
            h <- h + ld2_const$jacobian[j, k] * h_of(pj[k], r1)
          }
        }
        h
      })
    } else list()

    sbv <- if (q > 0L) sum(bd$B[r1, ]) else 0
    if (abs(1 - sbv) < 1e-10) {
      stop("the autoregressive polynomial is at the unit root; the filter has no starting level.",
           call. = FALSE)
    }
    om1 <- vals$omega$v[r1]
    om1_u <- row_of("omega", r1)
    om1_uu <- h_of("omega", r1)
    dbsum <- if (q > 0L) {
      Reduce(`+`, lapply(seq_len(q), function(j) db_of(j, r1)))
    } else numeric(mk)
    f0 <- om1 / (1 - sbv)
    f0_u <- om1_u / (1 - sbv) + om1 * dbsum / (1 - sbv)^2
    f0_uu <- om1_uu / (1 - sbv) +
      (outer(om1_u, dbsum) + outer(dbsum, om1_u)) / (1 - sbv)^2 +
      2 * om1 * outer(dbsum, dbsum) / (1 - sbv)^3
    if (q > 0L) {
      f0_uu <- f0_uu + om1 * Reduce(`+`, hb) / (1 - sbv)^2
    }

    groups[[l]] <- list(
      rows = rows, act = act,
      wom = wrow_mat("omega"), zom = zrow_mat("omega"),
      k2om = vals$omega$k2[rows],
      wa = lapply(aj, wrow_mat), za = lapply(aj, zrow_mat),
      k2a = lapply(aj, function(j) vals[[j]]$k2[rows]),
      db = lapply(seq_len(q), db_mat), hb = hb,
      f0 = f0, f0u = f0_u, f0uu = f0_uu)
  }
  A <- matrix(0, bp$n, p)
  for (i in seq_len(p)) A[, i] <- vals[[aj[i]]]$v
  B <- matrix(0, bp$n, q)
  for (j in seq_len(q)) B[, j] <- bd$B[, j]
  list(groups = groups, om = vals$omega$v, A = A, B = B)
}

# the filter of the submodel route: the general recursion with the exact
# derivative in the term's parameters propagated beside the state, carried
# per group on that group's active coordinates. The recursion runs
# compiled; .gas_filter_sub_r() below is the same loop in R, kept so the
# kernel has a twin that shares none of its code.
.gas_filter_sub <- function(term, eta, y, score, curvature, u,
                            fast = NULL, threads = 1L) {
  vals <- .gas_sub_values(term, u, "parameter")
  bd <- .gas_sub_b(term, vals)
  acts <- .gas_sub_active(term, vals)
  prep <- .gas_sub_prep(term, vals, bd, acts)
  res <- gas_filter_sub_cpp(eta, prep$groups, term@p, term@q,
                            prep$om, prep$A, prep$B, length(u),
                            score, curvature,
                            fast = fast, threads = as.integer(threads))
  colnames(res$jacobian) <- term_params(term)
  res
}

# the R twin of the compiled route above, compared at a tolerance in the
# tests (a compiler is free to contract a multiply-add, so identity is a
# statement about one compiler rather than about the arithmetic)
.gas_filter_sub_r <- function(term, eta, y, score, curvature, u) {
  bp <- term@blueprint
  p <- term@p
  q <- term@q
  np <- length(u)
  base <- .gas_base_params(p, q)
  vals <- .gas_sub_values(term, u, "parameter")
  bd <- .gas_sub_b(term, vals)
  acts <- .gas_sub_active(term, vals)
  aj <- if (p > 0L) paste0("alpha", seq_len(p)) else character(0)

  eta_out <- numeric(bp$n)
  jac <- matrix(0, bp$n, np)
  cv <- numeric(bp$n)
  for (l in seq_along(bp$order)) {
    rows <- bp$order[[l]]
    act <- acts[[l]]
    na <- length(act)
    # each base parameter's coordinates among the active ones, NA where a
    # coordinate is inactive here (its column is zero on these rows, so
    # its derivative contributes exactly nothing)
    pos <- lapply(base, function(j) match(vals[[j]]$idx, act))
    names(pos) <- base
    row_of <- function(j, r) {
      out <- numeric(na)
      at <- pos[[j]]
      ok <- !is.na(at)
      out[at[ok]] <- vals[[j]]$W[r, ok]
      out
    }
    db_of <- function(j, r) {
      out <- numeric(na)
      at <- match(bd$pidx, act)
      ok <- !is.na(at)
      out[at[ok]] <- bd$DB[[j]][r, ok]
      out
    }

    m <- length(rows)
    f <- numeric(m)
    s <- numeric(m)
    df <- matrix(0, m, na)
    ds <- matrix(0, m, na)

    # the starting level from the group's first observation: the
    # stationary level its own parameter values imply
    r1 <- rows[1L]
    sb <- if (q > 0L) sum(bd$B[r1, ]) else 0
    if (abs(1 - sb) < 1e-10) {
      stop("the autoregressive polynomial is at the unit root; the filter has no starting level.",
           call. = FALSE)
    }
    om1 <- vals$omega$v[r1]
    f0 <- om1 / (1 - sb)
    df0 <- row_of("omega", r1) / (1 - sb)
    if (q > 0L) {
      dbsum <- Reduce(`+`, lapply(seq_len(q), function(j) db_of(j, r1)))
      df0 <- df0 + om1 * dbsum / (1 - sb)^2
    }

    for (t in seq_len(m)) {
      r <- rows[t]
      ft <- vals$omega$v[r]
      dft <- row_of("omega", r)
      if (p > 0L) {
        for (i in seq_len(p)) {
          s_lag <- if (t - i >= 1L) s[t - i] else 0
          ds_lag <- if (t - i >= 1L) ds[t - i, ] else numeric(na)
          va <- vals[[aj[i]]]
          ft <- ft + va$v[r] * s_lag
          dft <- dft + va$v[r] * ds_lag + s_lag * row_of(aj[i], r)
        }
      }
      if (q > 0L) {
        for (j in seq_len(q)) {
          f_lag <- if (t - j >= 1L) f[t - j] else f0
          df_lag <- if (t - j >= 1L) df[t - j, ] else df0
          ft <- ft + bd$B[r, j] * f_lag
          dft <- dft + bd$B[r, j] * df_lag + f_lag * db_of(j, r)
        }
      }
      f[t] <- ft
      df[t, ] <- dft
      e_t <- eta[r] + ft
      s[t] <- score(e_t, r)
      cv[r] <- curvature(e_t, r)
      ds[t, ] <- cv[r] * dft
    }
    eta_out[rows] <- eta[rows] + f
    jac[rows, act] <- df
  }
  colnames(jac) <- term_params(term)
  list(eta = eta_out, jacobian = jac, curv = cv)
}

# the reverse pass of the submodel route, the coefficients read at the
# time whose equation they enter. The forward pass returns the curvature
# it read at every predictor, so the reverse pass looks it up instead of
# evaluating the callback a second time at the same points.
.gas_adjoint_sub <- function(term, eta, y, score, curvature, u, g,
                             fast = NULL, threads = 1L) {
  bp <- term@blueprint
  p <- term@p
  q <- term@q
  vals <- .gas_sub_values(term, u, "parameter")
  bd <- .gas_sub_b(term, vals)
  aj <- if (p > 0L) paste0("alpha", seq_len(p)) else character(0)
  cvv <- .gas_filter_sub(term, eta, y, score, curvature, u,
                         fast = fast, threads = threads)$curv

  deta <- numeric(bp$n)
  dscore <- numeric(bp$n)
  for (rows in bp$order) {
    k <- length(rows)
    fb <- numeric(k)
    sb <- numeric(k)
    for (t in rev(seq_len(k))) {
      r <- rows[t]
      eb <- g[r] + sb[t] * cvv[r]
      fb[t] <- fb[t] + eb
      deta[r] <- eb
      dscore[r] <- sb[t]
      if (q > 0L) {
        for (j in seq_len(q)) {
          if (t - j >= 1L) fb[t - j] <- fb[t - j] + bd$B[r, j] * fb[t]
        }
      }
      if (p > 0L) {
        for (i in seq_len(p)) {
          if (t - i >= 1L) sb[t - i] <- sb[t - i] + vals[[aj[i]]]$v[r] * fb[t]
        }
      }
    }
  }
  list(deta = deta, dscore = dscore)
}

# the second-order propagation of the submodel route: the forward jacobian
# of the predictor in the caller's unknowns and the contracted second
# derivative, the chart differentiated per observation, everything carried
# on each group's active set
.gas_curvature_sub <- function(term, eta, y, score, curvature, psiv, g,
                               seed, blocks, direction = NULL,
                               score_values = NULL, curvature_values = NULL,
                               blocks_data = NULL, threads = 1L) {
  bp <- term@blueprint
  p <- term@p
  q <- term@q
  third <- !is.null(direction)
  nm <- term_params(term)
  np <- length(nm)
  lay <- .gas_sub_layout(term)
  base <- .gas_base_params(p, q)
  # the caller passes the parameter scale; the chart is differentiated on
  # the unconstrained scale, which for a developed parameter is its own
  # coordinates already and for a scalar one is the link of its chart
  u <- psiv
  for (j in base) {
    if (is.null(bp$sub[[j]])) {
      u[lay$idx[[j]]] <- linkfunctions7::linkfun(.gas_param_link(term, j),
                                                 psiv[[lay$idx[[j]]]])
    }
  }
  vals <- .gas_sub_values(term, u, "zeta", second = TRUE, third = third)
  bd <- .gas_sub_b(term, vals)
  acts <- .gas_sub_active(term, vals)
  aj <- if (p > 0L) paste0("alpha", seq_len(p)) else character(0)
  pj <- if (q > 0L) paste0("pacf", seq_len(q)) else character(0)

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
  if (m < np) {
    stop(sprintf("'seed' has %d columns and the term has %d parameters.",
                 m, np), call. = FALSE)
  }
  zcol <- m - np + seq_len(np)
  if (third) {
    direction <- as.numeric(direction)
    if (length(direction) != m) {
      stop(sprintf("'direction' must have one value per unknown (%d).", m),
           call. = FALSE)
    }
  }
  varying_b <- q > 0L &&
    any(vapply(pj, function(j) vals[[j]]$developed, logical(1)))
  ld2_const <- if (q > 0L && !varying_b) {
    gas_levinson2(vapply(pj, function(j) vals[[j]]$v[1L], numeric(1)))
  } else NULL

  # THE COMPILED ROUTE (piano_parallel.txt, voce 5's blocked stage): second
  # order only, constant autoregressive charts, and a caller that declares
  # its callbacks as LOOKUPS (score_values/curvature_values) and its blocks
  # as data (H, D3, Vs, ap). Anything else keeps the R recursion below,
  # which stays the reference and the third-order implementation. The
  # kernel mirrors this function's arithmetic operation by operation; the
  # twin test compares the two routes, and threads run over GROUPS with
  # each group's W merged in group order on the main thread.
  if (!third && !varying_b && !is.null(score_values) &&
      !is.null(curvature_values) && !is.null(blocks_data)) {
    prep <- .gas_curv_prep(term, vals, bd, acts, lay, m, np, zcol,
                           ld2_const, base, aj, pj, p, q)
    res <- gas_curvature_sub_cpp(eta, prep$groups, p, q, prep$om, prep$A,
                                 prep$B, as.integer(blocks_data$ap),
                                 as.numeric(score_values),
                                 as.numeric(curvature_values),
                                 g, blocks_data$H, blocks_data$D3,
                                 blocks_data$Vs, seed,
                                 as.integer(threads))
    W <- res$curvature
    W <- (W + t(W)) / 2
    return(list(jacobian = res$jacobian, curvature = W))
  }

  # does the caller take the active set? A three-argument callback is the
  # earlier contract and is given the full row instead
  wants_act <- length(formals(blocks)) >= 4L

  D <- matrix(0, bp$n, m)
  W <- matrix(0, m, m)
  dP <- if (third) matrix(0, bp$n, m) else NULL
  for (l in seq_along(bp$order)) {
    rows <- bp$order[[l]]
    # this group's columns among the caller's unknowns: everything that is
    # not the term's, then the term's active coordinates
    act_t <- acts[[l]]
    act <- c(seq_len(m - np), zcol[act_t])
    mk <- length(act)
    pos <- lapply(base, function(j) {
      (m - np) + match(vals[[j]]$idx, act_t)
    })
    names(pos) <- base
    ppos <- (m - np) + match(bd$pidx, act_t)

    row_of <- function(j, r) {
      out <- numeric(mk)
      at <- pos[[j]]
      ok <- !is.na(at)
      out[at[ok]] <- vals[[j]]$W[r, ok]
      out
    }
    db_of <- function(j, r) {
      out <- numeric(mk)
      ok <- !is.na(ppos)
      out[ppos[ok]] <- bd$DB[[j]][r, ok]
      out
    }
    # the chart's second derivative of one value at one observation:
    # diagonal for a scalar coordinate, the outer product of the design
    # row for a developed one, on the active columns
    h_of <- function(j, r) {
      h <- matrix(0, mk, mk)
      vj <- vals[[j]]
      at <- pos[[j]]
      ok <- !is.na(at)
      if (!any(ok)) return(h)
      if (vj$developed) {
        z <- vj$Z[r, ok]
        h[at[ok], at[ok]] <- vj$k2[r] * outer(z, z)
      } else {
        h[at[ok], at[ok]] <- vj$k2[r]
      }
      h
    }
    # the second derivative of one autoregressive coefficient at one
    # observation: the map's curvature over the chained first derivatives,
    # plus each chart's own curvature under the map's jacobian
    hb_of <- function(r) {
      if (q == 0L) return(NULL)
      ld2 <- if (!varying_b) ld2_const else {
        gas_levinson2(vapply(pj, function(j) vals[[j]]$v[r], numeric(1)))
      }
      wrows <- lapply(pj, function(j) row_of(j, r))
      lapply(seq_len(q), function(j) {
        h <- matrix(0, mk, mk)
        for (k in seq_len(q)) {
          for (ll in seq_len(q)) {
            if (ld2$hessian[[j]][k, ll] != 0) {
              h <- h + ld2$hessian[[j]][k, ll] *
                outer(wrows[[k]], wrows[[ll]])
            }
          }
          if (ld2$jacobian[j, k] != 0) {
            h <- h + ld2$jacobian[j, k] * h_of(pj[k], r)
          }
        }
        h
      })
    }
    hb_const_l <- if (q > 0L && !varying_b) hb_of(rows[1L]) else NULL

    # the third derivative of one value at one observation, contracted
    # against the direction: a developed coordinate reads its chart's own
    # h''' at z'gamma along the design row, a scalar one is diagonal
    vk <- if (third) direction[act] else NULL
    t_of <- function(j, r) {
      h <- matrix(0, mk, mk)
      vj <- vals[[j]]
      at <- pos[[j]]
      ok <- !is.na(at)
      if (!any(ok)) return(h)
      if (vj$developed) {
        z <- vj$Z[r, ok]
        h[at[ok], at[ok]] <- vj$k3[r] * sum(z * vk[at[ok]]) * outer(z, z)
      } else {
        h[at[ok], at[ok]] <- vj$k3[r] * vk[at[ok]]
      }
      h
    }
    # and of one autoregressive coefficient: the Levinson-Durbin map's own
    # third derivative over the chained first derivatives, plus the two
    # places the product rule puts a chart's curvature and the one it puts
    # its third derivative
    tb_of <- function(r) {
      if (q == 0L) return(NULL)
      ld2 <- if (!varying_b) ld2_const else {
        gas_levinson2(vapply(pj, function(j) vals[[j]]$v[r], numeric(1)))
      }
      wrows <- lapply(pj, function(j) row_of(j, r))
      hrows <- lapply(pj, function(j) h_of(j, r))
      trows <- lapply(pj, function(j) t_of(j, r))
      hv <- lapply(hrows, function(H) as.numeric(H %*% vk))
      wdot <- vapply(wrows, function(W) sum(W * vk), numeric(1))
      ld3 <- gas_levinson3(vapply(pj, function(j) vals[[j]]$v[r], numeric(1)),
                           wdot)
      lapply(seq_len(q), function(j) {
        h <- matrix(0, mk, mk)
        hw <- as.numeric(ld2$hessian[[j]] %*% wdot)
        for (k in seq_len(q)) {
          for (ll in seq_len(q)) {
            if (ld3[[j]][k, ll] != 0) {
              h <- h + ld3[[j]][k, ll] * outer(wrows[[k]], wrows[[ll]])
            }
            if (ld2$hessian[[j]][k, ll] != 0) {
              h <- h + ld2$hessian[[j]][k, ll] *
                (outer(hv[[k]], wrows[[ll]]) + outer(wrows[[k]], hv[[ll]]))
            }
          }
          if (hw[[k]] != 0) h <- h + hw[[k]] * hrows[[k]]
          if (ld2$jacobian[j, k] != 0) {
            h <- h + ld2$jacobian[j, k] * trows[[k]]
          }
        }
        h
      })
    }
    tb_const_l <- if (third && q > 0L && !varying_b) tb_of(rows[1L]) else NULL

    r1 <- rows[1L]
    sbv <- if (q > 0L) sum(bd$B[r1, ]) else 0
    if (abs(1 - sbv) < 1e-10) {
      stop("the autoregressive polynomial is at the unit root; the filter has no starting level.",
           call. = FALSE)
    }
    om1 <- vals$omega$v[r1]
    om1_u <- row_of("omega", r1)
    om1_uu <- h_of("omega", r1)
    dbsum <- if (q > 0L) {
      Reduce(`+`, lapply(seq_len(q), function(j) db_of(j, r1)))
    } else numeric(mk)
    hb1 <- if (q > 0L) (if (varying_b) hb_of(r1) else hb_const_l) else NULL
    f0 <- om1 / (1 - sbv)
    f0_u <- om1_u / (1 - sbv) + om1 * dbsum / (1 - sbv)^2
    f0_uu <- om1_uu / (1 - sbv) +
      (outer(om1_u, dbsum) + outer(dbsum, om1_u)) / (1 - sbv)^2 +
      2 * om1 * outer(dbsum, dbsum) / (1 - sbv)^3
    if (q > 0L) {
      f0_uu <- f0_uu + om1 * Reduce(`+`, hb1) / (1 - sbv)^2
    }
    if (third) {
      # f0 = omega/(1 - S): the same expansion the scalar route carries,
      # with every quantity read at this group's first observation
      tb1 <- if (q > 0L) (if (varying_b) tb_of(r1) else tb_const_l) else NULL
      cst <- 1 / (1 - sbv)
      S_uu <- if (q > 0L) Reduce(`+`, hb1) else matrix(0, mk, mk)
      S_3 <- if (q > 0L) Reduce(`+`, tb1) else matrix(0, mk, mk)
      dS <- sum(dbsum * vk)
      dS2 <- as.numeric(S_uu %*% vk)
      dom <- sum(om1_u * vk)
      dom2 <- as.numeric(om1_uu %*% vk)
      om1_3 <- t_of("omega", r1)
      f0_3 <- om1_3 * cst + om1_uu * (cst^2 * dS) +
        2 * cst^3 * dS * (outer(om1_u, dbsum) + outer(dbsum, om1_u)) +
        cst^2 * (outer(dom2, dbsum) + outer(om1_u, dS2) +
                 outer(dS2, om1_u) + outer(dbsum, dom2)) +
        (2 * dom * cst^3 + 6 * om1 * cst^4 * dS) * outer(dbsum, dbsum) +
        2 * om1 * cst^3 * (outer(dS2, dbsum) + outer(dbsum, dS2)) +
        dom * cst^2 * S_uu + 2 * om1 * cst^3 * dS * S_uu + om1 * cst^2 * S_3
      df0 <- sum(f0_u * vk)
      dphi0 <- as.numeric(f0_uu %*% vk)
    }

    Wl <- matrix(0, mk, mk)
    k <- length(rows)
    f <- numeric(k)
    s <- numeric(k)
    F_ <- vector("list", k)
    Phi <- vector("list", k)
    Sd <- vector("list", k)
    Sdd <- vector("list", k)
    Psi <- if (third) vector("list", k) else NULL
    Sddd <- if (third) vector("list", k) else NULL
    for (t in seq_len(k)) {
      r <- rows[t]
      ft <- vals$omega$v[r]
      Ft <- row_of("omega", r)
      Pt <- h_of("omega", r)
      Tt <- if (third) t_of("omega", r) else NULL
      if (p > 0L) {
        for (i in seq_len(p)) {
          lag <- t - i
          s_l <- if (lag >= 1L) s[lag] else 0
          Sd_l <- if (lag >= 1L) Sd[[lag]] else numeric(mk)
          Sdd_l <- if (lag >= 1L) Sdd[[lag]] else matrix(0, mk, mk)
          va <- vals[[aj[i]]]
          a_u <- row_of(aj[i], r)
          a_uu <- h_of(aj[i], r)
          ft <- ft + va$v[r] * s_l
          Ft <- Ft + va$v[r] * Sd_l + s_l * a_u
          Pt <- Pt + va$v[r] * Sdd_l +
            outer(a_u, Sd_l) + outer(Sd_l, a_u) +
            s_l * a_uu
          if (third) {
            Sddd_l <- if (lag >= 1L) Sddd[[lag]] else matrix(0, mk, mk)
            dSd_l <- if (lag >= 1L) sum(Sd_l * vk) else 0
            dSdd_l <- if (lag >= 1L) as.numeric(Sdd_l %*% vk) else numeric(mk)
            da <- sum(a_u * vk)
            da2 <- as.numeric(a_uu %*% vk)
            Tt <- Tt + da * Sdd_l + va$v[r] * Sddd_l +
              outer(da2, Sd_l) + outer(Sd_l, da2) +
              outer(a_u, dSdd_l) + outer(dSdd_l, a_u) +
              dSd_l * a_uu + s_l * t_of(aj[i], r)
          }
        }
      }
      if (q > 0L) {
        hb <- if (varying_b) hb_of(r) else hb_const_l
        tb <- if (third) (if (varying_b) tb_of(r) else tb_const_l) else NULL
        for (j in seq_len(q)) {
          lag <- t - j
          f_l <- if (lag >= 1L) f[lag] else f0
          F_l <- if (lag >= 1L) F_[[lag]] else f0_u
          Phi_l <- if (lag >= 1L) Phi[[lag]] else f0_uu
          b_u <- db_of(j, r)
          ft <- ft + bd$B[r, j] * f_l
          Ft <- Ft + bd$B[r, j] * F_l + f_l * b_u
          Pt <- Pt + bd$B[r, j] * Phi_l +
            outer(b_u, F_l) + outer(F_l, b_u) + f_l * hb[[j]]
          if (third) {
            Psi_l <- if (lag >= 1L) Psi[[lag]] else f0_3
            dF_l <- if (lag >= 1L) sum(F_l * vk) else df0
            dPhi_l <- if (lag >= 1L) as.numeric(Phi_l %*% vk) else dphi0
            db <- sum(b_u * vk)
            db2 <- as.numeric(hb[[j]] %*% vk)
            Tt <- Tt + db * Phi_l + bd$B[r, j] * Psi_l +
              outer(db2, F_l) + outer(F_l, db2) +
              outer(b_u, dPhi_l) + outer(dPhi_l, b_u) +
              dF_l * hb[[j]] + f_l * tb[[j]]
          }
        }
      }
      f[t] <- ft
      F_[[t]] <- Ft
      Phi[[t]] <- Pt
      if (third) Psi[[t]] <- Tt

      e_t <- eta[r] + ft
      Dt <- seed[r, act] + Ft
      D[r, act] <- Dt
      Wl <- Wl + g[r] * (if (third) Tt else Pt)

      s[t] <- score(e_t, r)
      cv <- curvature(e_t, r)
      bl <- if (wants_act) blocks(e_t, r, Dt, act) else {
        full <- numeric(m)
        full[act] <- Dt
        b3 <- blocks(e_t, r, full)
        list(cross = b3$cross[act], M = b3$M[act, act, drop = FALSE],
             dcurv = b3$dcurv[act],
             N = if (third) b3$N[act, act, drop = FALSE] else NULL)
      }
      Sd[[t]] <- cv * Dt + bl$cross
      Sdd[[t]] <- cv * Pt + bl$M
      if (third) {
        dPhi_t <- as.numeric(Pt %*% vk)
        dP[r, act] <- dPhi_t
        Sddd[[t]] <- sum(bl$dcurv * vk) * Pt + cv * Tt + bl$N +
          outer(dPhi_t, bl$dcurv) + outer(bl$dcurv, dPhi_t)
      }
    }
    W[act, act] <- W[act, act] + Wl
  }
  W <- (W + t(W)) / 2
  if (third) return(list(jacobian = D, dphi = dP, curvature = W))
  list(jacobian = D, curvature = W)
}


#' The Filter's Parameters at Rows Outside the Fitting Data
#'
#' @description
#' A developed parameter's value at each new row, read through each
#' sub-term's own blueprint rather than rebuilt.
#'
#' @details
#' It is \code{\link{term_predict}} on every sub-term, which is what keeps
#' a basis or a set of contrasts from being relearned at other rows, and
#' then the same chart the fit used. The persistence goes through
#' Levinson-Durbin exactly as it does inside the filter.
#'
#' @param term A built score-driven term.
#' @param u The term's free vector.
#' @param newdata The rows to read at.
#'
#' @return A list with \code{om}, \code{A} and \code{B}, one row each per
#'   observation of \code{newdata}.
#'
#' @seealso \code{\link{term_continue}}
#'
#' @keywords internal
.gas_sub_new <- function(term, u, newdata) {
  bp <- term@blueprint
  p <- term@p
  q <- term@q
  nn <- nrow(newdata)
  base <- .gas_base_params(p, q)
  lay <- .gas_sub_layout(term)
  v <- list()
  for (j in base) {
    idx <- lay$idx[[j]]
    lk <- .gas_param_link(term, j)
    sj <- bp$sub[[j]]
    if (is.null(sj)) {
      v[[j]] <- rep(u[[idx]], nn)
    } else {
      Z <- .nl_bind(lapply(sj$terms, term_predict, newdata = newdata))
      v[[j]] <- linkfunctions7::linkinv(lk, as.numeric(as.matrix(Z) %*%
                                                         u[idx]))
    }
  }
  A <- matrix(0, nn, p)
  for (i in seq_len(p)) A[, i] <- v[[paste0("alpha", i)]]
  B <- matrix(0, nn, q)
  if (q > 0L) {
    pac <- vapply(seq_len(q), function(j) v[[paste0("pacf", j)]],
                  numeric(nn))
    pac <- matrix(pac, nn, q)
    for (r in seq_len(nn)) B[r, ] <- gas_levinson(pac[r, ])$phi
  }
  list(om = v$omega, A = A, B = B)
}
