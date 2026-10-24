#' @include gas_joint.R gas_submodels.R
NULL

# The joint score-driven filter with coefficients that vary by observation:
# the filter, with its exact derivative in the joint parameter vector, and
# the reverse pass. It is the joint analogue of .gas_filter_sub() and
# .gas_adjoint_sub() (gas_submodels.R). A member's parameter that carries a
# subformula is read at each row from .gas_sub_values(), the autoregressive
# coefficients from .gas_sub_b(), and the cross loadings, scalar or developed,
# are read per row as well. The
# curvature of the second and higher orders is in gas_joint_sub_core.R.
#
# The joint vector is the concatenation of the members' term_params(), member
# k starting after term@offsets[k]; inside a member the developed layout of
# .gas_sub_layout() comes first and the cross loadings follow, partner by
# partner and lag by lag.

# The per-member inputs of the filter at the parameters psi: the
# coefficients per row, their chart derivatives, the autoregressive
# coefficients per row and the cross loadings per row (n x npk*p, the column
# of partner j at lag i being j + npk (i - 1)) with their names in that order.
.gas_joint_sub_members <- function(term, psi) {
  psi <- .gas_joint_psi(term, psi)
  lapply(seq_along(term@params), function(k) {
    tm <- term@terms[[term@params[k]]]
    nm <- term_params(tm)
    u <- unlist(psi[[term@params[k]]][nm], use.names = FALSE)
    if (length(u) != length(nm) || anyNA(u)) {
      stop(sprintf("'psi$%s' must supply %s.", term@params[k],
                   paste(nm, collapse = ", ")), call. = FALSE)
    }
    p <- tm@p
    q <- tm@q
    n <- tm@blueprint$n
    off <- term@offsets[k]
    nl <- length(.gas_sub_layout(tm)$names)
    vals <- .gas_sub_values(tm, u, "parameter")
    bd <- .gas_sub_b(tm, vals)
    pk <- term@partners[[k]]
    npk <- length(pk)
    cnm <- character(0)
    for (i in seq_len(p)) {
      for (j in seq_len(npk)) {
        cnm <- c(cnm, paste0("kappa", i, "_", tm@cross_with[j]))
      }
    }
    cross <- matrix(0, n, npk * p)
    for (cc in seq_along(cnm)) cross[, cc] <- vals[[cnm[cc]]]$v
    aj <- if (p > 0L) paste0("kappa", seq_len(p)) else character(0)
    A <- matrix(0, n, p)
    for (i in seq_len(p)) A[, i] <- vals[[aj[i]]]$v
    B <- matrix(0, n, q)
    for (j in seq_len(q)) B[, j] <- bd$B[, j]
    list(tm = tm, p = p, q = q, off = off, vals = vals, bd = bd,
         om = vals$omega$v, A = A, B = B, partners = as.integer(pk),
         cross = cross, cross_names = cnm)
  })
}

# The groups of the compiled route. Each group is carried on the union of
# the members' active coordinates (moved by the member's offset, the cross
# loadings' included), and every member's derivative rows are scattered onto
# it once, from the single term's .gas_sub_prep(); the cross loadings' rows
# are scattered here.
.gas_joint_sub_groups <- function(term, mem) {
  bp <- term@blueprint
  K <- length(mem)
  acts <- lapply(mem, function(m) .gas_sub_active(m$tm, m$vals))
  preps <- lapply(seq_len(K), function(k)
    .gas_sub_prep(mem[[k]]$tm, mem[[k]]$vals, mem[[k]]$bd, acts[[k]]))
  lapply(seq_along(bp$order), function(l) {
    rows <- bp$order[[l]]
    jact <- sort(unique(as.integer(
      unlist(lapply(seq_len(K), function(k) mem[[k]]$off + acts[[k]][[l]])))))
    na <- length(jact)
    members <- lapply(seq_len(K), function(k) {
      g <- preps[[k]]$groups[[l]]
      at <- match(mem[[k]]$off + g$act, jact)
      sc <- function(M) {
        out <- matrix(0, nrow(M), na)
        if (length(at)) out[, at] <- M
        out
      }
      df0 <- numeric(na)
      if (length(at)) df0[at] <- g$df0
      dcr <- lapply(mem[[k]]$cross_names, function(nmj) {
        vj <- mem[[k]]$vals[[nmj]]
        pos <- match(mem[[k]]$off + vj$idx, jact)
        ok <- !is.na(pos)
        M <- matrix(0, length(rows), na)
        if (any(ok)) M[, pos[ok]] <- vj$W[rows, ok, drop = FALSE]
        M
      })
      list(dom = sc(g$dom), da = lapply(g$da, sc), dbl = lapply(g$dbl, sc),
           f0 = g$f0, df0 = df0, dcr = dcr)
    })
    list(rows = as.integer(rows), act = as.integer(jact), mem = members)
  })
}

# the arguments and the value of term_filter() of GasJoint
.gas_joint_sub_filter <- function(term, eta, y, score, curvature, psi,
                                  fast = NULL, threads = 1L) {
  bp <- term@blueprint
  K <- length(term@params)
  n <- bp$n
  eta <- matrix(as.numeric(eta), n, K)
  nm <- term_params(term)
  M <- length(nm)
  mem <- .gas_joint_sub_members(term, psi)

  if (!is.null(fast)) {
    cm <- lapply(mem, function(m)
      list(p = m$p, q = m$q, om = m$om, A = m$A, B = m$B,
           partners = m$partners, cross = m$cross))
    out <- gas_joint_sub_filter_cpp(eta, cm, .gas_joint_sub_groups(term, mem),
                                    M, fast, as.integer(threads))
    if (!is.null(out)) {
      jac <- lapply(out$jacobian, function(J) { colnames(J) <- nm; J })
      names(jac) <- term@params
      eta_out <- out$eta
      colnames(eta_out) <- term@params
      return(list(eta = eta_out, jacobian = jac, curv = out$curv))
    }
  }

  # the R route: the loop of term_filter() of GasJoint with the coefficients
  # read per row and their derivative rows dense over the M joint columns
  dense <- function(k, j, r) {
    out <- numeric(M)
    vj <- mem[[k]]$vals[[j]]
    out[mem[[k]]$off + vj$idx] <- vj$W[r, ]
    out
  }
  dense_b <- function(k, j, r) {
    out <- numeric(M)
    out[mem[[k]]$off + mem[[k]]$bd$pidx] <- mem[[k]]$bd$DB[[j]][r, ]
    out
  }
  eta_out <- eta
  jac <- lapply(seq_len(K), function(k) matrix(0, n, M))
  cv <- array(0, c(n, K, K))
  for (rows in bp$order) {
    T_ <- length(rows)
    f <- matrix(0, K, T_)
    u <- matrix(0, K, T_)
    Fd <- array(0, c(K, T_, M))
    Sd <- array(0, c(K, T_, M))
    r1 <- rows[1L]
    f0 <- numeric(K)
    df0 <- matrix(0, K, M)
    for (k in seq_len(K)) {
      m <- mem[[k]]
      sb <- if (m$q > 0L) sum(m$bd$B[r1, ]) else 0
      if (abs(1 - sb) < 1e-10) {
        stop("the autoregressive polynomial is at the unit root; the filter has no starting level.",
             call. = FALSE)
      }
      f0[k] <- m$om[r1] / (1 - sb)
      df0[k, ] <- dense(k, "omega", r1) / (1 - sb)
      if (m$q > 0L) {
        dbsum <- Reduce(`+`, lapply(seq_len(m$q), function(j)
          dense_b(k, j, r1)))
        df0[k, ] <- df0[k, ] + m$om[r1] * dbsum / (1 - sb)^2
      }
    }
    for (t in seq_len(T_)) {
      r <- rows[t]
      for (k in seq_len(K)) {
        m <- mem[[k]]
        ft <- m$om[r]
        Ft <- dense(k, "omega", r)
        for (i in seq_len(m$p)) {
          if (t - i < 1L) next
          a_r <- m$A[r, i]
          ft <- ft + a_r * u[k, t - i]
          Ft <- Ft + a_r * Sd[k, t - i, ] +
            u[k, t - i] * dense(k, paste0("kappa", i), r)
          for (j in seq_along(m$partners)) {
            l <- m$partners[j]
            cc <- j + length(m$partners) * (i - 1L)
            cji <- m$cross[r, cc]
            ft <- ft + cji * u[l, t - i]
            Ft <- Ft + cji * Sd[l, t - i, ] +
              u[l, t - i] * dense(k, m$cross_names[cc], r)
          }
        }
        for (j in seq_len(m$q)) {
          if (t - j >= 1L) {
            f_l <- f[k, t - j]
            F_l <- Fd[k, t - j, ]
          } else {
            f_l <- f0[k]
            F_l <- df0[k, ]
          }
          ft <- ft + m$B[r, j] * f_l
          Ft <- Ft + m$B[r, j] * F_l + f_l * dense_b(k, j, r)
        }
        f[k, t] <- ft
        Fd[k, t, ] <- Ft
      }
      e_t <- eta[r, ] + f[, t]
      u[, t] <- score(e_t, r)
      C <- matrix(curvature(e_t, r), K, K)
      cv[r, , ] <- C
      Sd[, t, ] <- C %*% matrix(Fd[, t, ], K, M)
    }
    eta_out[rows, ] <- eta[rows, , drop = FALSE] + t(f)
    for (k in seq_len(K)) jac[[k]][rows, ] <- Fd[k, , ]
  }
  jac <- lapply(jac, function(J) { colnames(J) <- nm; J })
  names(jac) <- term@params
  colnames(eta_out) <- term@params
  list(eta = eta_out, jacobian = jac, curv = cv)
}

# the arguments and the value the adjoint of GasJoint expects: cv is the
# stored curvature, an n x K x K array, and g the direct derivative in the
# filtered predictors; the reverse pass reads each coefficient at the row
# whose equation it enters
.gas_joint_sub_adjoint <- function(term, cv, g, psi) {
  bp <- term@blueprint
  K <- length(term@params)
  n <- bp$n
  g <- matrix(as.numeric(g), n, K)
  mem <- .gas_joint_sub_members(term, psi)
  cm <- lapply(mem, function(m)
    list(p = m$p, q = m$q, A = m$A, B = m$B, partners = m$partners,
         cross = m$cross))
  out <- gas_joint_sub_adjoint_cpp(as.numeric(cv), g, cm, bp$order)
  colnames(out$deta) <- colnames(out$dscore) <- term@params
  out
}
