#' @include gas.R gas_submodels.R
NULL

# The compiled route of the score-driven recursion's second, third and fourth
# orders (src/gas_curvature_gen.cpp). Both routes of .gas_curvature_core()
# hand it the same per-group description: every value of the recursion
# (the level, the loadings, and, where they vary by observation, the partial
# autocorrelations) as its per-row value, chained derivative row W, raw
# chart row Z and link derivatives k2..k4; the autoregressive coefficients'
# chained rows and, where they are constant, their second to fourth
# derivative blocks; and the starting level with its derivatives, which the
# R code computes once per group. The per-row recursion runs in C++.

# Whether the caller supplied what the compiled route reads: the score and
# the curvature as LOOKUPS and the model's pieces as data, with the fourth
# derivative for the third order and the fifth for the fourth.
.gas_gen_ok <- function(score_values, curvature_values, blocks_data, nd) {
  !is.null(score_values) && !is.null(curvature_values) &&
    !is.null(blocks_data) &&
    (nd < 1L || !is.null(blocks_data$D4)) &&
    (nd < 2L || !is.null(blocks_data$D5))
}

# One value of the recursion on a group, as the kernel reads it.
.gas_gen_val <- function(v, W, Z, k2, k3 = NULL, k4 = NULL) {
  list(v = as.numeric(v), W = W, Z = Z, k2 = as.numeric(k2),
       k3 = if (is.null(k3)) numeric(0) else as.numeric(k3),
       k4 = if (is.null(k4)) numeric(0) else as.numeric(k4))
}

# Calls the kernel and returns what the R route returns at the same order.
.gas_gen_run <- function(eta, groups, p, q, nd, om, A, B, score_values,
                         curvature_values, g, seed, blocks_data, threads) {
  empty <- matrix(0, 0, 0)
  res <- gas_curvature_gen_cpp(
    as.numeric(eta), groups, p, q, nd, as.numeric(om), A, B,
    as.integer(blocks_data$ap), as.numeric(score_values),
    as.numeric(curvature_values), as.numeric(g), blocks_data$H,
    blocks_data$D3,
    if (is.null(blocks_data$D4)) empty else blocks_data$D4,
    if (is.null(blocks_data$D5)) empty else blocks_data$D5,
    blocks_data$Vs, seed, as.integer(threads))
  # the accumulation is not symmetric to the bit; the R route's convention
  W <- (res$curvature + t(res$curvature)) / 2
  if (nd >= 2L) {
    return(list(jacobian = res$jacobian, dphi = res$dphi, dpsi = res$dpsi,
                curvature = W))
  }
  if (nd == 1L) {
    return(list(jacobian = res$jacobian, dphi = res$dphi[[1L]],
                curvature = W))
  }
  list(jacobian = res$jacobian, curvature = W)
}

# The scalar route's groups: every unknown is active, the chart of each
# scalar parameter is its link's derivatives at its coordinate, and the
# autoregressive coefficients are constant, their derivative blocks read off
# the shared preparation.
.gas_gen_groups_scalar <- function(term, zeta, links, base, zcol, gp, bp, m,
                                   dirs, nd) {
  p <- term@p
  q <- term@q
  third <- nd >= 1L
  fourth <- nd >= 2L
  unit_val <- function(j, rows) {
    i <- match(j, base)
    lk <- links[[j]]
    z <- zeta[[j]]
    mg <- length(rows)
    e <- numeric(m)
    e[zcol[i]] <- 1
    Z <- matrix(e, mg, m, byrow = TRUE)
    W <- linkfunctions7::dlinkinv(lk, z) * Z
    .gas_gen_val(rep(linkfunctions7::linkinv(lk, z), mg), W, Z,
                 rep(linkfunctions7::d2linkinv(lk, z), mg),
                 if (third) rep(linkfunctions7::d3linkinv(lk, z), mg),
                 if (fourth) rep(linkfunctions7::d4linkinv(lk, z), mg))
  }
  aj <- if (p > 0L) paste0("kappa", seq_len(p)) else character(0)
  lapply(bp$order, function(rows) {
    mg <- length(rows)
    out <- list(rows = as.integer(rows), act = seq_len(m),
                om = unit_val("omega", rows),
                al = lapply(aj, unit_val, rows = rows),
                varying_b = FALSE,
                db = lapply(seq_len(q), function(j)
                  matrix(gp$b_u[[j]], mg, m, byrow = TRUE)),
                hb = if (q > 0L) gp$b_uu else list(),
                f0 = gp$f0, f0u = gp$f0_u, f0uu = gp$f0_uu)
    if (third) {
      out$vks <- matrix(unlist(lapply(dirs, function(d) d)), m, nd)
      out$tb <- if (q > 0L) gp$t_b else lapply(seq_len(nd), function(d) list())
      out$f0_3 <- gp$f0_3
      out$df0 <- unlist(gp$df0)
      out$dphi0 <- matrix(unlist(gp$dphi0), m, nd)
    }
    if (fourth) {
      out$qb <- if (q > 0L) gp$q_b else list()
      out$f0_4 <- gp$f0_4
    }
    out
  })
}
