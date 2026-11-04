#' @include block_deriv3.R
NULL

#' How a Term's Block Moves Along Four Directions
#'
#' @description
#' Returns the fourth derivative of a term's block along four directions,
#' one entry per observation and column of the block, the companion of
#' [term_block_deriv3()] at the next order.
#'
#' @details
#' The base method returns zeros. For [nl()] the entry is the fifth
#' derivative of the predictor along the column and the four directions,
#' written out over the 52 partitions of the five: grouped by the block that
#' holds the column, of size five (one partition), four (four), three (six
#' with the other two directions together and six with them apart), two (four
#' choices, each with the other three directions in one block, in a pair and a
#' single in three ways, or apart) and one (the fifteen partitions of the four
#' directions). Each block of a partition carries the derivative of one
#' inverse link to its size and each partition the derivative of \eqn{f} to
#' its number of blocks. A parameter whose development moves is rejected, as
#' in [term_block_deriv3()]. For [scs()] and [scte()] the entry is
#' \eqn{G_{il}W_{ij}b^{(5)}(\gamma_j)\,v_ju_js_jr_j}, as in
#' [term_block_deriv3()] one order more.
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
#' @inheritParams term_block_deriv3
#' @param v,u,s,r Numeric vectors as long as the term's coefficients.
#'
#' @return A numeric matrix, one row per observation and one column per
#'   coefficient of the term.
#'
#' @examples
#' dd <- data.frame(x = seq(0.2, 3, length.out = 20))
#' dd$y <- 2 * exp(-1.3 * dd$x)
#' b <- term_build(nl(~ a * exp(-r * x), start = list(a = 2, r = 1.3)), dd)
#' dim(term_block_deriv4(b, v = c(1, 0), u = c(0, 1), s = c(1, 1),
#'                       r = c(1, -1)))
#'
#' @seealso [term_block_deriv3()]
#'
#' @export
#' @aliases term_block_deriv4.model_term term_block_deriv4.NlTerm
#'   term_block_deriv4.ScsTerm term_block_deriv4.ScteTerm
term_block_deriv4 <- S7::new_generic(
  "term_block_deriv4", "term",
  function(term, coef = NULL, v, u, s, r, ...) S7::S7_dispatch())

S7::method(term_block_deriv4, model_term) <- function(term, coef = NULL, v, u,
                                                      s, r, ...) {
  d <- dim(term_matrix(term))
  matrix(0, d[1L], d[2L])
}

S7::method(term_block_deriv4, NlTerm) <- function(term, coef = NULL, v, u, s,
                                                  r, ...) {
  bp <- term@build_info
  if (!length(bp)) stop("the term is not built.", call. = FALSE)
  cf <- .nl_coef(bp, coef)
  for (arg in list(list(v, "v"), list(u, "u"), list(s, "s"), list(r, "r"))) {
    if (length(arg[[1L]]) != bp$ncoef) {
      stop(sprintf("'%s' must have length %d, the term's coefficients.",
                   arg[[2L]], bp$ncoef), call. = FALSE)
    }
  }
  params <- bp$params
  np <- length(params)
  if (any(vapply(params, function(p) isTRUE(bp$mobile[[p]]), logical(1)))) {
    stop(paste0("the fourth derivative of the block is not written for a ",
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
  h1 <- th$deriv
  h2 <- th$deriv2
  h3 <- th$deriv3
  h4 <- lapply(params, function(p)
    linkfunctions7::d4linkinv(bp$links[[p]], eta[[p]]))
  h5 <- lapply(params, function(p)
    linkfunctions7::d5linkinv(bp$links[[p]], eta[[p]]))
  names(h4) <- params
  names(h5) <- params
  fk <- lapply(1:5, function(k) .nl_forder(bp, th$value, k))
  F <- function(k, ix) as.numeric(fk[[k]][[.nl_key(params, ix)]])
  A <- .nl_direction(bp, as.numeric(v), th$X)
  B <- .nl_direction(bp, as.numeric(u), th$X)
  C <- .nl_direction(bp, as.numeric(s), th$X)
  D <- .nl_direction(bp, as.numeric(r), th$X)
  zero <- 0 * bp$one
  ip <- seq_len(np)
  out <- matrix(0, bp$n, bp$ncoef)
  for (i1 in ip) {
    # the sums over the parameters of the blocks that do not hold the
    # column, one per shape: a block of size m on parameter q carries
    # h^(m)_q times its directions there
    S1 <- function(x) {
      acc <- zero
      for (q in ip) acc <- acc + F(2, c(i1, q)) * h1[[q]] * x[[q]]
      acc
    }
    S2 <- function(x, y) {
      acc <- zero
      for (q in ip) acc <- acc + F(2, c(i1, q)) * h2[[q]] * x[[q]] * y[[q]]
      acc
    }
    S3 <- function(x, y, z) {
      acc <- zero
      for (q in ip) {
        acc <- acc + F(2, c(i1, q)) * h3[[q]] * x[[q]] * y[[q]] * z[[q]]
      }
      acc
    }
    S4 <- function(x, y, z, t) {
      acc <- zero
      for (q in ip) {
        acc <- acc + F(2, c(i1, q)) * h4[[q]] * x[[q]] * y[[q]] * z[[q]] *
          t[[q]]
      }
      acc
    }
    R11 <- function(x, y) {
      acc <- zero
      for (q in ip) for (k in ip) {
        acc <- acc + F(3, c(i1, q, k)) * h1[[q]] * x[[q]] * h1[[k]] * y[[k]]
      }
      acc
    }
    R21 <- function(x, y, z) {
      acc <- zero
      for (q in ip) for (k in ip) {
        acc <- acc + F(3, c(i1, q, k)) * h2[[q]] * x[[q]] * y[[q]] *
          h1[[k]] * z[[k]]
      }
      acc
    }
    R31 <- function(x, y, z, t) {
      acc <- zero
      for (q in ip) for (k in ip) {
        acc <- acc + F(3, c(i1, q, k)) * h3[[q]] * x[[q]] * y[[q]] * z[[q]] *
          h1[[k]] * t[[k]]
      }
      acc
    }
    R22 <- function(x, y, z, t) {
      acc <- zero
      for (q in ip) for (k in ip) {
        acc <- acc + F(3, c(i1, q, k)) * h2[[q]] * x[[q]] * y[[q]] *
          h2[[k]] * z[[k]] * t[[k]]
      }
      acc
    }
    T111 <- function(x, y, z) {
      acc <- zero
      for (q in ip) for (k in ip) for (m in ip) {
        acc <- acc + F(4, c(i1, q, k, m)) * h1[[q]] * x[[q]] * h1[[k]] *
          y[[k]] * h1[[m]] * z[[m]]
      }
      acc
    }
    T211 <- function(x, y, z, t) {
      acc <- zero
      for (q in ip) for (k in ip) for (m in ip) {
        acc <- acc + F(4, c(i1, q, k, m)) * h2[[q]] * x[[q]] * y[[q]] *
          h1[[k]] * z[[k]] * h1[[m]] * t[[m]]
      }
      acc
    }
    U1111 <- function(x, y, z, t) {
      acc <- zero
      for (q in ip) for (k in ip) for (m in ip) for (o in ip) {
        acc <- acc + F(5, c(i1, q, k, m, o)) * h1[[q]] * x[[q]] * h1[[k]] *
          y[[k]] * h1[[m]] * z[[m]] * h1[[o]] * t[[o]]
      }
      acc
    }
    a <- A[[i1]]
    b <- B[[i1]]
    c1 <- C[[i1]]
    d <- D[[i1]]
    # the block holding the column, of size five
    acc <- F(1, i1) * h5[[i1]] * a * b * c1 * d
    # of size four: one direction alone
    acc <- acc + h4[[i1]] * (a * b * c1 * S1(D) + a * b * d * S1(C) +
                               a * c1 * d * S1(B) + b * c1 * d * S1(A))
    # of size three: the other two together, then apart
    acc <- acc + h3[[i1]] * (
      a * b * (S2(C, D) + R11(C, D)) + a * c1 * (S2(B, D) + R11(B, D)) +
        a * d * (S2(B, C) + R11(B, C)) + b * c1 * (S2(A, D) + R11(A, D)) +
        b * d * (S2(A, C) + R11(A, C)) + c1 * d * (S2(A, B) + R11(A, B)))
    # of size two: the other three in one block, a pair and a single, apart
    two <- function(x, y, z) {
      S3(x, y, z) + R21(x, y, z) + R21(x, z, y) + R21(y, z, x) + T111(x, y, z)
    }
    acc <- acc + h2[[i1]] * (a * two(B, C, D) + b * two(A, C, D) +
                               c1 * two(A, B, D) + d * two(A, B, C))
    # of size one: the fifteen partitions of the four directions
    one <- S4(A, B, C, D) +
      R31(A, B, C, D) + R31(A, B, D, C) + R31(A, C, D, B) + R31(B, C, D, A) +
      R22(A, B, C, D) + R22(A, C, B, D) + R22(A, D, B, C) +
      T211(A, B, C, D) + T211(A, C, B, D) + T211(A, D, B, C) +
      T211(B, C, A, D) + T211(B, D, A, C) + T211(C, D, A, B) +
      U1111(A, B, C, D)
    acc <- acc + h1[[i1]] * one
    idx <- bp$index[[params[i1]]]
    Z <- th$X[[params[i1]]]
    out[, idx] <- if (is.null(Z)) acc else as.matrix(Z * acc)
  }
  out
}
