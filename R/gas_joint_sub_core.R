#' @include gas_joint_core.R gas_submodels.R
NULL

# The second, third and fourth orders of a joint score-driven filter whose
# coefficients vary by observation: the members' subformulas make an omega, a
# score loading or a partial autocorrelation psi_{j,t} = h_j(z_t' gamma_j),
# and the recursion's coefficients are then read at the row they multiply.
# .gas_joint_core() hands the call here when .gas_joint_has_sub(term) holds;
# the result has the structure that .gas_joint_core() returns: `jacobian`
# and `curvature`, with one direction `dphi`, with two `dphi` and `dpsi`.
#
# The preparation is that of the single term's submodel route
# (.gas_curvature_sub()), done for every member: the values of the base
# parameters at every row with their chart derivatives
# (.gas_sub_values()), and each group's active coordinates
# (.gas_sub_active()). A group's active unknowns are the unknowns that are not
# parameters of the filter, every cross loading, and the coordinates of the
# members that are not identically zero on the group, so that the kernels run
# on a set whose size does not grow with the number of groups. The recursion,
# the jets of the coefficients at every row, the autoregressive coefficients'
# Levinson-Durbin map and the starting levels are in
# src/gas_joint_sub_curv.cpp.
.gas_joint_sub_core <- function(term, psi, g, seed, blocks, direction = NULL,
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
  seed <- lapply(seed, function(s) {
    s <- as.matrix(s)
    storage.mode(s) <- "double"
    s
  })
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
  U1 <- blocks$U1
  U2 <- blocks$U2
  U3 <- if (nd >= 1L) blocks$U3 else NULL
  U4 <- if (nd >= 2L) blocks$U4 else NULL
  if ((nd >= 1L && is.null(U3)) || (nd >= 2L && is.null(U4))) {
    stop("'blocks' must carry U3 for the third order and U4 for the fourth.",
         call. = FALSE)
  }
  u_mat <- matrix(as.numeric(blocks$u), n, K)

  # the filter's parameters are the LAST M unknowns, in the order of
  # term_params()
  zoff <- m - M
  psi <- .gas_joint_psi(term, psi)

  mem <- vector("list", K)
  acts <- vector("list", K)
  lay <- vector("list", K)
  bases <- vector("list", K)
  nlay <- integer(K)
  ncross <- integer(K)
  for (k in seq_len(K)) {
    pk <- term@params[k]
    tm <- term@terms[[pk]]
    p <- tm@p
    q <- tm@q
    base <- .gas_base_params(p, q)
    bases[[k]] <- base
    lay[[k]] <- .gas_sub_layout(tm)
    nl <- length(lay[[k]]$names)
    nlay[k] <- nl
    nm <- term_params(tm)
    v <- unlist(psi[[pk]][nm])
    if (length(v) != length(nm) || anyNA(v)) {
      stop(sprintf("'psi$%s' must supply %s.", pk, paste(nm, collapse = ", ")),
           call. = FALSE)
    }
    v <- as.numeric(v)
    # the term's coordinates on the mixed scale of .gas_sub_values(): the
    # coefficients of a developed parameter as they are, a scalar parameter
    # through the link of its chart
    u <- v[seq_len(nl)]
    sub <- tm@blueprint$sub
    for (j in base) {
      if (is.null(sub[[j]])) {
        ix <- lay[[k]]$idx[[j]]
        u[ix] <- linkfunctions7::linkfun(.gas_param_link(tm, j), u[ix])
      }
    }
    vals <- .gas_sub_values(tm, u, "zeta", second = TRUE, third = nd >= 1L,
                            fourth = nd >= 2L)
    acts[[k]] <- .gas_sub_active(tm, vals)
    # the cross loadings stay scalar and ride the identity link: partner by
    # partner, lag by lag, after the developed layout
    partners <- term@partners[[k]]
    np <- length(partners)
    ncross[k] <- np * p
    cl <- matrix(0, np, p)
    ccol <- matrix(0L, np, p)
    for (j in seq_len(np)) {
      for (i in seq_len(p)) {
        cl[j, i] <- v[nl + (j - 1L) * p + i]
        ccol[j, i] <- term@offsets[k] + nl + (j - 1L) * p + i
      }
    }
    vl <- lapply(base, function(j) {
      vj <- vals[[j]]
      Z <- if (isTRUE(vj$developed)) vj$Z else matrix(1, n, 1L)
      list(v = as.numeric(vj$v),
           W = matrix(as.numeric(vj$W), n),
           Z = matrix(as.numeric(Z), n),
           k2 = as.numeric(vj$k2),
           k3 = if (nd >= 1L) as.numeric(vj$k3) else numeric(0),
           k4 = if (nd >= 2L) as.numeric(vj$k4) else numeric(0))
    })
    mem[[k]] <- list(p = p, q = q, partners = as.integer(partners),
                     cross = as.numeric(cl), col_cross = as.integer(ccol),
                     developed = vapply(base, function(j)
                       isTRUE(vals[[j]]$developed), logical(1)),
                     vals = vl)
  }

  # the active unknowns of each group, and where each value's active
  # coordinates sit among them
  ord <- bp$order
  groups <- vector("list", length(ord))
  for (l in seq_along(ord)) {
    cols <- lapply(seq_len(K), function(k) {
      loc <- c(acts[[k]][[l]], nlay[k] + seq_len(ncross[k]))
      zoff + term@offsets[k] + loc
    })
    act <- sort(unique(c(seq_len(zoff), unlist(cols))))
    gv <- lapply(seq_len(K), function(k) {
      lapply(bases[[k]], function(j) {
        ids <- lay[[k]]$idx[[j]]
        sel <- which(ids %in% acts[[k]][[l]])
        list(at = as.integer(match(zoff + term@offsets[k] + ids[sel], act)),
             sel = as.integer(sel))
      })
    })
    groups[[l]] <- list(rows = as.integer(ord[[l]]), act = as.integer(act),
                        gv = gv)
  }

  out <- if (nd == 0L) {
    gas_joint_sub_curv_cpp(seed, mem, groups, as.integer(fp), g, u_mat,
                           as.numeric(U1), as.numeric(U2), as.integer(m),
                           as.integer(zoff), as.integer(threads))
  } else if (nd == 1L) {
    gas_joint_sub_third_cpp(seed, mem, groups, as.integer(fp), g, u_mat,
                            as.numeric(U1), as.numeric(U2), as.numeric(U3),
                            as.numeric(dirs[[1L]]), as.integer(m),
                            as.integer(zoff), as.integer(threads))
  } else {
    gas_joint_sub_fourth_cpp(seed, mem, groups, as.integer(fp), g, u_mat,
                             as.numeric(U1), as.numeric(U2), as.numeric(U3),
                             as.numeric(U4), as.numeric(dirs[[1L]]),
                             as.numeric(dirs[[2L]]), as.integer(m),
                             as.integer(zoff), as.integer(threads))
  }
  names(out$jacobian) <- term@params
  if (nd >= 1L) names(out$dphi) <- term@params
  if (nd >= 2L) names(out$dpsi) <- term@params
  out
}
