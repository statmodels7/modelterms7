#' @include term_classes.R generics.R
NULL

#' @title S7 Class for Smooth Terms
#' @name SmoothTerm
#'
#' @description
#' A subclass of \code{\link{additive_term}} for a penalized smooth of one
#' or more covariates: a \pkg{basis7} expansion with a roughness penalty
#' on its coefficients. Constructed by \code{\link{s}} for one covariate
#' and by \code{\link{te}} for several.
#'
#' @inheritParams additive_term
#' @param vars The expressions of the covariates being smoothed.
#' @param by An optional expression the smooth varies with.
#' @param spec The construction settings: the basis, its dimension and
#'   degree, and whether the linear part is carried separately.
#' @param sparse Whether the block is a \code{dgCMatrix}, which only a factor
#'   \code{by} admits. See \code{\link{s}}.
#'
#' @return An object of class \code{SmoothTerm}.
#'
#' @seealso \code{\link{s}}, \code{\link{te}}
#' @examples
#' S7::S7_inherits(s(x), SmoothTerm)
#' @export
SmoothTerm <- S7::new_class(
  name = "SmoothTerm",
  parent = additive_term,
  properties = list(
    vars = S7::class_list,
    by = S7::class_any,
    spec = S7::class_list,
    # WHERE the block is sparse: a factor `by` puts each row in the block of
    # its own level, a density of 1/m. The basis itself is dense by
    # construction -- the Demmler-Reinsch rotation makes it so -- and so is
    # a numeric `by`, which merely multiplies it. NULL until the build
    # settles it; the settled value is what the blueprint carries.
    sparse = S7::class_any
  )
)

#' Penalized Smooth of One Covariate
#'
#' @description
#' A smooth function of a covariate, expanded in a \pkg{basis7} basis and
#' penalized for roughness. The default construction is a cubic B-spline
#' basis under the Demmler-Reinsch reparametrization, which separates the
#' linear effect from the nonlinear deviation and turns the roughness
#' penalty into the identity on the deviation.
#'
#' @details
#' The block has one column for the linear effect, centered and scaled,
#' followed by the reparametrized basis. That ordering is what the penalty
#' reads: it is the quadratic penalty of
#' \eqn{\mathrm{diag}(0, 1, \dots, 1)}, rank deficient by exactly one, so
#' the linear effect is unpenalized and the deviation is shrunk towards
#' zero. As the smoothing parameter grows the fit approaches a straight
#' line rather than a constant, and \code{\link{edf}} falls towards one.
#'
#' The Demmler-Reinsch construction (\cite{demmler1975}, as used for
#' effect selection by \cite{bach2024}) is empirical: it takes the inner
#' product at the observed covariate values, so it is built when the term
#' is, and the transform is stored in the blueprint. Prediction at new
#' points is the parent basis evaluated there and multiplied by the same
#' transform, so the separation of the linear from the nonlinear part
#' holds on new data as it does on old.
#'
#' \subsection{Varying the smooth by another variable}{
#' With \code{by} a factor the term carries one smooth per level, the
#' block being the smooth multiplied by each level's indicator and the
#' penalty the same matrix repeated blockwise, so one smoothing parameter
#' governs every level. With \code{by} numeric the term is a
#' varying-coefficient one: the smooth multiplies that variable, and the
#' fitted function is the coefficient of \code{by} as it changes with the
#' covariate.
#'
#' A factor \code{by} is where a smooth's block can be SPARSE: each row sits
#' in the block of its own level and nowhere else, a density of \eqn{1/m}.
#' \code{sparse = TRUE} builds it that way rather than building the dense
#' matrix and compressing it -- measured at 2000 rows, \eqn{k = 10} and 200
#' levels, 0.35 MB against 28.93 MB, the numbers identical. \code{sparse =
#' NULL}, the default, settles it at build from the size of the block, the
#' dense form holding \eqn{n m k} cells against \eqn{n k} non-zeros; see
#' \code{\link{.resolve_sparse}} for the threshold and what it was measured
#' against. An explicit \code{TRUE} is refused without a factor \code{by},
#' where there would be nothing to build on: the basis is dense by
#' construction, the Demmler-Reinsch rotation making it so, and a numeric
#' \code{by} merely multiplies it.
#'
#' The block alone is sparse. The PENALTY of a factor \code{by} is the same
#' matrix repeated blockwise, and \pkg{penalties7} returns it dense -- 25.92
#' MB at those sizes, at a density of 0.0005 -- which is a property of that
#' package's contract rather than of this construction.
#' }
#'
#' @param x The covariate, an expression evaluated in the data.
#' @param by An optional factor or numeric variable; see Details.
#' @param k The basis dimension before reparametrization. Defaults to 10.
#' @param degree The spline degree. Defaults to 3, a cubic spline.
#' @param basis An optional \pkg{basis7} basis to use in place of the
#'   default B-spline; its range is taken as given.
#' @param linear Whether the linear effect is carried in the block, and
#'   left unpenalized there. Defaults to \code{TRUE}.
#' @param label A single non-empty string prefixed to the coefficient
#'   names. Defaults to a name built from the covariate.
#' @param lambda The smoothing parameter, held at the value given and
#'   ESTIMATED when left \code{NULL}, which is the default. An
#'   anisotropic tensor product carries one per margin, so a vector
#'   of that length, or a named one holding some of them.
#' @param sparse Whether the block is built as a \code{dgCMatrix}.
#'   \code{NULL}, the default, settles it at build from the size of the block.
#'   Only a FACTOR \code{by} admits it, each row sitting in the block of its
#'   own level; without one an explicit \code{TRUE} is refused rather than
#'   ignored. See Details.
#'
#' @return An object of class \code{\link{SmoothTerm}} (a specification;
#'   see \code{\link{term_build}}).
#'
#' @references
#' Demmler, A. and Reinsch, C. (1975). Oscillation matrices with spline
#' smoothing. \emph{Numerische Mathematik}, 24, 375--382.
#'
#' Bach, P. and Klein, N. (2024). Bayesian effect selection in additive
#' models with an application to time-to-event data.
#'
#' @examples
#' set.seed(1)
#' dd <- data.frame(x = sort(runif(80)))
#' dd$y <- sin(2 * pi * dd$x) + rnorm(80, sd = 0.2)
#' built <- term_build(s(x, k = 8), dd)
#' term_npar(built)
#' term_penalty(built)@params
#'
#' @seealso \code{\link{te}}, \code{\link{random}}, \code{\link{nl}}
#' @export
s <- function(x, by = NULL, k = 10, degree = 3, basis = NULL,
              linear = TRUE, label = NULL, lambda = NULL,
              sparse = NULL) {
  xe <- substitute(x)
  .smooth_spec(list(xe), substitute(by), k, degree,
               if (is.null(basis)) NULL else list(basis), linear, label,
               sprintf("s(%s)", deparse(xe)), lambda, sparse = sparse)
}

#' Penalized Smooth of Several Covariates
#'
#' @description
#' A tensor-product smooth: a \pkg{basis7} basis in each covariate,
#' combined by \code{\link[basis7]{tensor_basis}}, with a roughness
#' penalty on the product coefficients.
#'
#' @details
#' The block is the tensor basis evaluated at the covariates, and the
#' penalty is built from the marginal roughness penalties carried into the
#' product, \eqn{P_v = I \otimes \cdots \otimes P_v \otimes \cdots
#' \otimes I}, each penalizing curvature in one direction.
#'
#' With \code{anisotropic = TRUE}, the default, those components enter
#' \code{\link[penalties7]{additive_penalty}} and keep a smoothing
#' parameter each, so the surface may be rough in one direction and smooth
#' in another -- which is the usual reason for fitting a tensor smooth
#' rather than an isotropic one. With \code{anisotropic = FALSE} they are
#' summed first and one parameter governs the total, which costs one
#' hyperparameter instead of one per margin.
#'
#' The marginal bases are not reparametrized, so unlike \code{\link{s}}
#' the marginal linear effects are not separated out: the null space of
#' the tensor penalty contains them, and they are shrunk towards no
#' surface at all rather than towards a plane.
#'
#' \subsection{Centering}{
#' The tensor product of the marginal bases contains the constant, which
#' the penalty's null space contains as well, so beside an intercept the
#' block would be rank deficient by exactly one and the penalty would not
#' cover the deficiency. The block therefore carries the sum-to-zero
#' constraint over the observed covariates
#' (\code{\link[basis7]{constrain_basis}}): the term has one column fewer
#' than the product of its marginal dimensions, every column sums to zero
#' over the data it was built on, and the penalty follows by congruence
#' with its rank unchanged, the direction removed having been one of its
#' null directions. The transform is stored in the blueprint and reapplied
#' by \code{\link{term_predict}}, as the Demmler-Reinsch transform of
#' \code{\link{s}} is.
#'
#' The level of the surface is then the model's intercept, so a formula
#' that removes it (\code{y ~ te(x, z) - 1}) fits a surface constrained to
#' average zero over the data.
#' }
#'
#' @param ... The covariates, expressions evaluated in the data, at least
#'   two of them.
#' @param by An optional factor or numeric variable, as in \code{\link{s}}.
#' @param k The basis dimension per margin, recycled to the number of
#'   covariates. Defaults to 5.
#' @param degree The spline degree per margin, recycled. Defaults to 3.
#' @param bases An optional list of \pkg{basis7} bases, one per covariate,
#'   used in place of the default B-splines.
#' @param anisotropic Keep a smoothing parameter per margin? Defaults to
#'   \code{TRUE}.
#' @param label A single non-empty string prefixed to the coefficient
#'   names. Defaults to a name built from the covariates.
#' @param lambda The smoothing parameter, held at the value given and
#'   ESTIMATED when left \code{NULL}, which is the default. An
#'   anisotropic tensor product carries one per margin, so a vector
#'   of that length, or a named one holding some of them.
#' @param sparse Whether the block is built as a \code{dgCMatrix}.
#'   \code{NULL}, the default, settles it at build from the size of the block.
#'   Only a FACTOR \code{by} admits it, each row sitting in the block of its
#'   own level; without one an explicit \code{TRUE} is refused rather than
#'   ignored. See \code{\link{s}}.
#'
#' @return An object of class \code{\link{SmoothTerm}} (a specification;
#'   see \code{\link{term_build}}).
#'
#' @examples
#' set.seed(2)
#' dd <- data.frame(x = runif(120), z = runif(120))
#' dd$y <- dd$x * dd$z + rnorm(120, sd = 0.1)
#' built <- term_build(te(x, z, k = 4), dd)
#' term_npar(built)
#'
#' @references
#' Wood, S. N. (2006). Low-rank scale-invariant tensor product smooths for
#' generalized additive mixed models. \emph{Biometrics} 62, 1025-1036.
#'
#' Wood, S. N. (2017). \emph{Generalized Additive Models: An Introduction
#' with R}, 2nd edition. Chapman and Hall/CRC.
#'
#' @seealso \code{\link{s}}, \code{\link{random}}, \code{\link{nl}}
#' @export
te <- function(..., by = NULL, k = 5, degree = 3, bases = NULL,
               anisotropic = TRUE, label = NULL, lambda = NULL,
               sparse = NULL) {
  vars <- as.list(substitute(list(...)))[-1L]
  if (length(vars) < 2L) {
    stop("'te' needs at least two covariates; use s() for one.",
         call. = FALSE)
  }
  if (!is.logical(anisotropic) || length(anisotropic) != 1L ||
      is.na(anisotropic)) {
    stop("'anisotropic' must be TRUE or FALSE.", call. = FALSE)
  }
  # ANISOTROPIC is one smoothing parameter per margin, so that is how many
  # names there are to hold; isotropic is the single one of a quadratic
  # penalty. Either way the names are the penalty's own, which is what the
  # summary prints and what the fit keys its hyperparameters by.
  nms <- if (anisotropic) paste0("lambda", seq_along(vars)) else "lambda"
  sp <- .smooth_spec(vars, substitute(by), k, degree, bases, FALSE, label,
                     sprintf("te(%s)",
                             paste(vapply(vars, deparse, character(1)),
                                   collapse = ",")),
                     lambda, nms, sparse = sparse)
  sp@spec$anisotropic <- anisotropic
  sp
}

.smooth_spec <- function(vars, by, k, degree, bases, linear, label,
                         default_label, lambda = NULL, names = "lambda",
                         sparse = NULL) {
  nv <- length(vars)
  chk <- function(v, nm, lo) {
    if (!is.numeric(v) || anyNA(v) || any(v < lo) || any(v != round(v))) {
      stop(sprintf("'%s' must be a whole number of at least %d.", nm, lo),
           call. = FALSE)
    }
    as.integer(rep_len(v, nv))
  }
  k <- chk(k, "k", 3L)
  degree <- chk(degree, "degree", 1L)
  # a B-spline of degree m spans at least m + 1 functions, and finding that
  # out inside term_build() would report it several frames from the call
  if (is.null(bases) && any(k <= degree)) {
    stop(sprintf(paste("'k' must exceed 'degree': a spline of degree %d needs",
                       "at least %d basis functions."),
                 max(degree), max(degree) + 1L), call. = FALSE)
  }
  if (!is.null(bases)) {
    if (!is.list(bases) || length(bases) != nv ||
        !all(vapply(bases, function(b) S7::S7_inherits(b, basis7::basis),
                    logical(1)))) {
      stop("'bases' must be a list of one basis7 basis per covariate.",
           call. = FALSE)
    }
  }
  if (is.null(label)) label <- default_label
  if (!is.character(label) || length(label) != 1L || is.na(label) ||
      !nzchar(label)) {
    stop("'label' must be a single non-empty character string.",
         call. = FALSE)
  }
  if (!is.null(sparse) &&
      (!is.logical(sparse) || length(sparse) != 1L || is.na(sparse))) {
    stop(paste0("'sparse' must be TRUE, FALSE, or NULL to settle it from",
                " the design."), call. = FALSE)
  }
  # ASKED AT CONSTRUCTION, where the caller can see it: the sparsity of a
  # smooth comes from a factor `by`, whose indicators put each row in one
  # block. Without a `by` there is only the basis, which is dense by
  # construction, so the argument would have nothing to build on.
  if (isTRUE(sparse) && is.null(by)) {
    stop(paste0("'sparse' has nothing to build on here: a smooth's basis is",
                " dense by\n  construction. Sparsity comes from a FACTOR",
                " 'by', whose indicators put each\n  row in the block of its",
                " own level."), call. = FALSE)
  }
  SmoothTerm(label = label, vars = vars, by = by, sparse = sparse,
             spec = list(k = k, degree = degree, bases = bases,
                         linear = isTRUE(linear)),
             hyper = smooth_hyper(lambda, names, label),
             X = NULL, coef_names = character(0),
             blueprint = list(), penalty = NULL)
}

# the covariate values, one column per variable
.smooth_x <- function(vars, data) {
  out <- lapply(vars, function(e) as.numeric(eval(e, data, baseenv())))
  n <- unique(lengths(out))
  if (length(n) != 1L || n != nrow(data)) {
    stop("every covariate of a smooth must give one value per row.",
         call. = FALSE)
  }
  if (any(vapply(out, anyNA, logical(1)))) {
    stop("a smooth's covariates must not contain missing values.",
         call. = FALSE)
  }
  out
}

# One copy of the basis per level of a factor `by`: row i is non-zero only
# in the block of its own level, which is the shape .random_block() has and
# a density of 1/m. Built sparse where the term asks for it, rather than
# built dense and compressed.
.smooth_by_block <- function(Z, g, m, sparse = FALSE) {
  n <- nrow(Z)
  cols <- ncol(Z)
  gi <- as.integer(g)
  if (!isTRUE(sparse)) {
    G <- outer(g, levels(g), `==`) * 1
    Zb <- matrix(0, n, m * cols)
    for (l in seq_len(m)) {
      Zb[, (l - 1L) * cols + seq_len(cols)] <- G[, l] * Z
    }
    return(Zb)
  }
  if (anyNA(gi)) {
    stop(paste0("a 'by' with missing values has no block to place those rows",
                " in;\n  sparse = TRUE cannot build them."), call. = FALSE)
  }
  Matrix::sparseMatrix(
    i = rep(seq_len(n), each = cols),
    j = as.vector(t(outer((gi - 1L) * cols, seq_len(cols), `+`))),
    x = as.vector(t(as.matrix(Z))),
    dims = c(n, m * cols))
}

# the by variable: a factor gives one copy of the block per level, a numeric
# multiplies it
.smooth_by <- function(by, data, levels = NULL) {
  if (is.null(by)) return(NULL)
  v <- eval(by, data, baseenv())
  if (is.factor(v) || is.character(v) || is.logical(v)) {
    f <- if (is.null(levels)) factor(v) else factor(v, levels = levels)
    if (any(is.na(f) & !is.na(v))) {
      stop("a 'by' level absent at build time cannot be predicted.",
           call. = FALSE)
    }
    list(kind = "factor", value = f, levels = levels(f))
  } else {
    list(kind = "numeric", value = as.numeric(v))
  }
}

S7::method(term_build, SmoothTerm) <- function(term, data, ...) {
  xs <- .smooth_x(term@vars, data)
  sp <- term@spec
  nv <- length(xs)

  marg <- lapply(seq_len(nv), function(j) {
    if (!is.null(sp$bases)) return(sp$bases[[j]])
    r <- range(xs[[j]])
    pad <- diff(r) * 0.001 + .Machine$double.eps
    basis7::bspline_basis(lower = r[1] - pad, upper = r[2] + pad,
                          dimension = sp$k[j], degree = sp$degree[j])
  })

  if (nv == 1L) {
    b <- marg[[1L]]
    d <- basis7::dr_basis(b, xs[[1L]])
    Z <- basis7::basis_eval(d, xs[[1L]])
    nm <- paste0("z", seq_len(ncol(Z)))
    if (sp$linear) {
      ctr <- mean(xs[[1L]])
      scl <- stats::sd(xs[[1L]])
      if (!is.finite(scl) || scl == 0) scl <- 1
      Z <- cbind((xs[[1L]] - ctr) / scl, Z)
      nm <- c("lin", nm)
      pen_diag <- c(0, rep(1, ncol(Z) - 1L))
      lin <- list(center = ctr, scale = scl)
    } else {
      pen_diag <- rep(1, ncol(Z))
      lin <- NULL
    }
    P <- diag(pen_diag, length(pen_diag))
    core <- list(kind = "dr", dr = d, lin = lin)
  } else {
    tb <- basis7::tensor_basis(marg)
    xm <- do.call(cbind, xs)
    # the marginal roughness penalties carried into the product: the second
    # derivative Gram of each margin, the identity in the others
    dims <- vapply(marg, function(b) b@dimension, integer(1))
    comps <- lapply(seq_len(nv), function(j) {
      Pj <- basis7::basis_gram(marg[[j]], order = 2L)
      Pj <- Pj / max(1, max(abs(Pj)))
      blocks <- lapply(seq_len(nv), function(i)
        if (i == j) Pj else diag(dims[i]))
      # tensor_basis varies the FIRST margin fastest, so the Kronecker
      # product is taken in reverse order
      Reduce(kronecker, rev(blocks))
    })
    # The tensor product contains the constant, which the penalty's null space
    # also contains, so a model carrying an intercept would be rank deficient
    # by exactly one and nothing in the penalty would cover it. The sum-to-zero
    # constraint over the observed covariates removes that direction: every
    # column of the constrained block sums to zero, hence is orthogonal to an
    # intercept column. The penalty follows by congruence, and since the
    # direction removed lies in its null space the rank is unchanged.
    tb <- basis7::constrain_basis(tb, colSums(basis7::basis_eval(tb, xm)))
    tmat <- tb@transform
    comps <- lapply(comps, function(Pk) {
      M <- crossprod(tmat, Pk %*% tmat)
      (M + t(M)) / 2
    })
    Z <- basis7::basis_eval(tb, xm)
    nm <- paste0("z", seq_len(ncol(Z)))
    P <- if (isTRUE(sp$anisotropic)) comps else Reduce(`+`, comps)
    core <- list(kind = "tensor", basis = tb)
  }

  by_blocks <- 1L
  # settled here and carried to the blueprint: without a factor `by` there is
  # nothing sparse to build, so the answer is FALSE whatever was asked
  sp <- FALSE
  by <- .smooth_by(term@by, data)
  if (!is.null(by)) {
    if (identical(by$kind, "factor")) {
      m <- nlevels(by$value)
      nm <- as.character(t(outer(levels(by$value), nm,
                                 function(a, b) paste(a, b, sep = "."))))
      # the block is m copies of a basis of ncol(Z), one row per level, so
      # the dense form holds n * m * ncol(Z) cells against n * ncol(Z)
      # non-zeros -- the same product .resolve_sparse() reads for a formula
      sp <- .resolve_sparse(term@sparse, nrow(Z), m * ncol(Z))
      Z <- .smooth_by_block(Z, by$value, m, sp)
      # ONE COPY PER LEVEL, and the penalty says so rather than being handed
      # the assembled product. What a quadratic penalty needs from its matrix
      # -- the rank, the log pseudo-determinant, a basis of the null space --
      # all follow from one block, the eigenvalues of I (x) P being P's
      # repeated m times, so nothing of size (mk)^2 is ever decomposed.
      # Measured at m = 200 over a basis of ten, 565 times faster to build
      # and 0.03 MB stored against 32.
      #
      # The ANISOTROPIC branch still assembles: additive_penalty() reads one
      # eigendecomposition of the SUM of its components, which is not a
      # blockwise quantity, so the same shortcut does not apply to it.
      by_blocks <- m
      if (is.list(P)) {
        P <- lapply(P, function(Pk) kronecker(diag(m), Pk))
        by_blocks <- 1L
      }
    } else {
      # a numeric `by` MULTIPLIES the basis, so the block is as dense as the
      # basis is, and there is nothing for a sparse storage to hold back
      if (isTRUE(term@sparse)) {
        stop(paste0("'sparse' has nothing to build on here: a numeric 'by'",
                    " multiplies the basis,\n  so the block is as dense as",
                    " the basis. Sparsity comes from a FACTOR 'by', whose",
                    "\n  indicators put each row in one block."),
             call. = FALSE)
      }
      Z <- by$value * Z
    }
  }

  cn <- paste(term@label, nm, sep = ".")
  colnames(Z) <- cn
  rownames(Z) <- NULL
  term@X <- Z
  term@coef_names <- cn
  term@blueprint <- list(core = core, marg = marg, spec = sp,
                         vars = term@vars, by = term@by,
                         by_levels = if (!is.null(by) &&
                                         identical(by$kind, "factor"))
                           by$levels else NULL,
                         sparse = sp,
                         nblock = ncol(Z))
  term@penalty <- if (is.list(P)) penalties7::additive_penalty(P)
                  else penalties7::quadratic_penalty(P, blocks = by_blocks)
  term
}

S7::method(term_predict, SmoothTerm) <- function(term, newdata, ...) {
  .assert_built(term)
  bp <- term@blueprint
  xs <- .smooth_x(bp$vars, newdata)
  core <- bp$core

  if (identical(core$kind, "dr")) {
    Z <- basis7::basis_eval(core$dr, xs[[1L]])
    if (!is.null(core$lin)) {
      Z <- cbind((xs[[1L]] - core$lin$center) / core$lin$scale, Z)
    }
  } else {
    Z <- basis7::basis_eval(core$basis, do.call(cbind, xs))
  }

  by <- .smooth_by(bp$by, newdata, levels = bp$by_levels)
  if (!is.null(by)) {
    if (identical(by$kind, "factor")) {
      # the STORAGE is part of the blueprint: a prediction that densified
      # would spend at new data what the build was careful not to
      Z <- .smooth_by_block(Z, by$value, length(bp$by_levels),
                            isTRUE(bp$sparse))
    } else {
      Z <- by$value * Z
    }
  }
  colnames(Z) <- term@coef_names
  rownames(Z) <- NULL
  Z
}
