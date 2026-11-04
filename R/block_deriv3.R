#' @include nonlinear.R
NULL

#' How a Term's Block Moves Along Three Directions
#'
#' @description
#' Returns \eqn{\sum_{c,d,e} v_c u_d s_e\,\partial^3 X_{ij}/\partial\beta_c
#' \partial\beta_d\partial\beta_e}, one entry per observation and column of
#' the block: the third derivative of a term's block along three directions.
#' The Jeffreys penalty of a block inside the term's subformula reads it in
#' the derivative of the outer criterion.
#'
#' @details
#' The base method returns zeros, which is exact for a block that does not
#' move. For [nl()], with \eqn{\theta_p = h_p(Z_p\beta_p)} and the block's
#' column \eqn{c} in parameter \eqn{p_1}, the entry is the fourth derivative of
#' the predictor along \eqn{(e_c, v, u, s)}, the sum over the fifteen
#' partitions of the four directions, each block of a partition carrying the
#' derivative of one link and each partition the derivative of \eqn{f} of its
#' number of blocks. Writing \eqn{a_q, b_q, c_q} for the three directions
#' carried onto parameter \eqn{q} (\eqn{Z_q v_q} and so on), \eqn{h^{(k)}_q}
#' for the derivatives of its inverse link and \eqn{f_{q\dots}} for those of
#' \eqn{f},
#' \deqn{\begin{aligned}
#' X^{(3)}[t,c] = Z_{p_1}[t,c]\Big[&f_{p_1}h^{(4)}_{p_1}a_{p_1}b_{p_1}c_{p_1}
#'  + h^{(3)}_{p_1}\big(a_{p_1}b_{p_1}S(c) + a_{p_1}c_{p_1}S(b)
#'    + b_{p_1}c_{p_1}S(a)\big)\\
#'  &+ h''_{p_1}\big(a_{p_1}(Q(b,c) + R(b,c)) + b_{p_1}(Q(a,c) + R(a,c))
#'    + c_{p_1}(Q(a,b) + R(a,b))\big)\\
#'  &+ h'_{p_1}\big(\textstyle\sum_q f_{p_1q}h^{(3)}_qa_qb_qc_q
#'    + P(ab,c) + P(ac,b) + P(bc,a)
#'    + \sum_{q,r,s}f_{p_1qrs}h'_qa_q\,h'_rb_r\,h'_sc_s\big)\Big],
#' \end{aligned}}
#' with \eqn{S(x) = \sum_q f_{p_1q}h'_qx_q},
#' \eqn{Q(x,y) = \sum_q f_{p_1q}h''_qx_qy_q},
#' \eqn{R(x,y) = \sum_{q,r}f_{p_1qr}h'_qx_q\,h'_ry_r} and
#' \eqn{P(xy,z) = \sum_{q,r}f_{p_1qr}h''_qx_qy_q\,h'_rz_r}. A parameter whose
#' development moves with the coefficients has a design of its own that moves,
#' and its third derivative is not written: the method signals an error there.
#'
#' A class whose block moves with its coefficients must implement its own
#' method of this generic, because the zeros of the base method are wrong
#' for such a block. statmodels7 reads this generic in the exact derivatives
#' of an outer criterion and uses those derivatives only for a class with a
#' method of its own; for any other class it computes them by differences.
#' It recognizes a method of its own by comparing the method registered for
#' the class with the base method, so a class that registers the base method
#' itself under its own name is treated as having a method of its own, and
#' the zeros enter the derivatives.
#'
#' For [scs()] and [scte()] the column \eqn{j} of level \eqn{l} is
#' \eqn{G_{il}W_{ij}b'(\gamma_j)}, which depends on the coefficient
#' \eqn{\gamma_j} alone, so the entry is
#' \eqn{G_{il}W_{ij}b^{(4)}(\gamma_j)\,v_ju_js_j}, with \eqn{b} the map from
#' the coefficients to the basis coefficients.
#'
#' @param term A built term.
#' @param coef The coefficients, one per coefficient of the term, or `NULL`
#'   for the coefficients that the term carries.
#' @param v,u,s Numeric vectors as long as the term's coefficients.
#' @param ... Passed to methods.
#'
#' @return A numeric matrix, one row per observation and one column per
#'   coefficient of the term.
#'
#' @examples
#' dd <- data.frame(x = seq(0.2, 3, length.out = 20))
#' dd$y <- 2 * exp(-1.3 * dd$x)
#' b <- term_build(nl(~ a * exp(-r * x), start = list(a = 2, r = 1.3)), dd)
#' dim(term_block_deriv3(b, v = c(1, 0), u = c(0, 1), s = c(1, 1)))
#'
#' @seealso [term_block_deriv2()], [term_block_deriv()]
#'
#' @export
#' @aliases term_block_deriv3.model_term term_block_deriv3.NlTerm
#'   term_block_deriv3.ScsTerm term_block_deriv3.ScteTerm
term_block_deriv3 <- S7::new_generic(
  "term_block_deriv3", "term",
  function(term, coef = NULL, v, u, s, ...) S7::S7_dispatch())

S7::method(term_block_deriv3, model_term) <- function(term, coef = NULL, v, u,
                                                      s, ...) {
  d <- dim(term_matrix(term))
  matrix(0, d[1L], d[2L])
}

S7::method(term_block_deriv3, NlTerm) <- function(term, coef = NULL, v, u, s,
                                                  ...) {
  bp <- term@build_info
  if (!length(bp)) stop("the term is not built.", call. = FALSE)
  cf <- .nl_coef(bp, coef)
  for (arg in list(list(v, "v"), list(u, "u"), list(s, "s"))) {
    if (length(arg[[1L]]) != bp$ncoef) {
      stop(sprintf("'%s' must have length %d, the term's coefficients.",
                   arg[[2L]], bp$ncoef), call. = FALSE)
    }
  }
  params <- bp$params
  np <- length(params)
  if (any(vapply(params, function(p) isTRUE(bp$mobile[[p]]), logical(1)))) {
    stop(paste0("the third derivative of the block is not written for a ",
                "parameter whose\n  development moves with the ",
                "coefficients."), call. = FALSE)
  }
  th <- .nl_theta(bp, cf, order = 3L)
  eta <- lapply(params, function(p) {
    idx <- bp$index[[p]]
    Z <- th$X[[p]]
    if (is.null(Z)) rep(cf[idx], bp$n) else as.numeric(Z %*% cf[idx])
  })
  names(eta) <- params
  h4 <- lapply(params, function(p)
    linkfunctions7::d4linkinv(bp$links[[p]], eta[[p]]))
  names(h4) <- params
  f1 <- .nl_forder(bp, th$value, 1L)
  f2 <- .nl_forder(bp, th$value, 2L)
  f3 <- .nl_forder(bp, th$value, 3L)
  f4 <- .nl_forder(bp, th$value, 4L)
  h1 <- th$deriv
  h2 <- th$deriv2
  h3 <- th$deriv3
  A <- .nl_direction(bp, as.numeric(v), th$X)
  B <- .nl_direction(bp, as.numeric(u), th$X)
  C <- .nl_direction(bp, as.numeric(s), th$X)
  F <- function(fk, ix) as.numeric(fk[[.nl_key(params, ix)]])
  zero <- 0 * bp$one
  out <- matrix(0, bp$n, bp$ncoef)
  for (i1 in seq_len(np)) {
    p1 <- params[i1]
    S <- function(x) {
      acc <- zero
      for (q in seq_len(np)) acc <- acc + F(f2, c(i1, q)) * h1[[q]] * x[[q]]
      acc
    }
    Q <- function(x, y) {
      acc <- zero
      for (q in seq_len(np)) {
        acc <- acc + F(f2, c(i1, q)) * h2[[q]] * x[[q]] * y[[q]]
      }
      acc
    }
    Rf <- function(x, y) {
      acc <- zero
      for (q in seq_len(np)) {
        for (r in seq_len(np)) {
          acc <- acc + F(f3, c(i1, q, r)) * h1[[q]] * x[[q]] * h1[[r]] * y[[r]]
        }
      }
      acc
    }
    P <- function(x, y, z) {
      acc <- zero
      for (q in seq_len(np)) {
        for (r in seq_len(np)) {
          acc <- acc + F(f3, c(i1, q, r)) * h2[[q]] * x[[q]] * y[[q]] *
            h1[[r]] * z[[r]]
        }
      }
      acc
    }
    a1 <- A[[i1]]
    b1 <- B[[i1]]
    c1 <- C[[i1]]
    t3 <- zero
    for (q in seq_len(np)) {
      t3 <- t3 + F(f2, c(i1, q)) * h3[[q]] * A[[q]] * B[[q]] * C[[q]]
    }
    t4 <- zero
    for (q in seq_len(np)) {
      for (r in seq_len(np)) {
        for (k in seq_len(np)) {
          t4 <- t4 + F(f4, c(i1, q, r, k)) * h1[[q]] * A[[q]] * h1[[r]] *
            B[[r]] * h1[[k]] * C[[k]]
        }
      }
    }
    acc <- F(f1, i1) * h4[[i1]] * a1 * b1 * c1 +
      h3[[i1]] * (a1 * b1 * S(C) + a1 * c1 * S(B) + b1 * c1 * S(A)) +
      h2[[i1]] * (a1 * (Q(B, C) + Rf(B, C)) + b1 * (Q(A, C) + Rf(A, C)) +
                    c1 * (Q(A, B) + Rf(A, B))) +
      h1[[i1]] * (t3 + P(A, B, C) + P(A, C, B) + P(B, C, A) + t4)
    idx <- bp$index[[p1]]
    Z <- th$X[[p1]]
    out[, idx] <- if (is.null(Z)) acc else as.matrix(Z * acc)
  }
  out
}
