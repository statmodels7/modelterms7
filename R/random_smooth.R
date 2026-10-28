#' @include random.R smooth.R
NULL

#' @title S7 Class for a Smooth Random Effect
#' @name RandomSmoothTerm
#'
#' @description
#' `RandomSmoothTerm` is the subclass of [RandomTerm()] whose within-group
#' design is a smooth: one random curve or surface per level of the grouping
#' variable. [random()] constructs it when the left side of the bar is an [s()]
#' or a [te()] call, as in `random(~ s(x) | g)`.
#'
#' @details
#' # The two properties of its own
#'
#' `smooth` is the unbuilt [SmoothTerm()] read from the left side of the bar,
#' which carries the basis, the roughness penalty and the null space through
#' its \pkg{basis7} smoother. `intercept` records whether each group also has
#' a random level: `TRUE` for `~ s(x) | g` and `FALSE` for `~ 0 + s(x) | g`.
#'
#' Every other property, and every generic that [RandomTerm()] implements,
#' has the meaning documented there: the block is sparse and ordered group by
#' group, a level that was not present at build time is predicted at a zero
#' effect on request, and [term_within()] returns the design that one group's
#' effect multiplies.
#'
#' @inheritParams RandomTerm
#' @param smooth The unbuilt [SmoothTerm()] of the left side of the bar.
#' @param intercept A single logical: whether each group has a random level.
#'
#' @return An S7 object of class `RandomSmoothTerm`, inheriting from
#'   [RandomTerm()], [additive_term()] and [model_term()].
#'
#' @seealso [random()], the constructor, whose section on smooth random
#'   effects gives the model; [s()] and [te()] for the smooth itself.
#'
#' @examples
#' tm <- random(~ s(x, smoother = basis7::bspline_smooth(k = 6)) | g)
#' S7::S7_inherits(tm, RandomSmoothTerm)
#' S7::S7_inherits(tm, RandomTerm)
#' tm@intercept
#'
#' @export
RandomSmoothTerm <- S7::new_class(
  name = "RandomSmoothTerm",
  parent = RandomTerm,
  properties = list(
    smooth = S7::class_any,
    intercept = S7::class_logical
  )
)


#' The Smooth on the Left Side of a Bar
#'
#' @description
#' Reads the left side of a bar formula and returns the smooth written there,
#' with whether the groups also carry a random level, or `NULL` where the left
#' side holds no [s()] or [te()] call.
#'
#' @details
#' The left side is split on `+`. A smooth random effect holds exactly one
#' [s()] or [te()] call, optionally beside `1`, `0` or `-1`, which follow the
#' intercept convention of a model formula. Any other covariate beside the
#' smooth signals an error naming the spelling that writes it as a term of its
#' own.
#'
#' The call is evaluated with this package's [s()] and [te()] in place of
#' whatever the formula's environment binds to those names, so that a session
#' with \pkg{mgcv} attached builds the same term.
#'
#' The smooth may not carry a `by` variable, a `sparse` setting, held
#' hyperparameters or an `id`. The grouping is the variable after the bar,
#' and the hyperparameters belong to the random effect, where `hyper` and `id`
#' of [random()] name them.
#'
#' @param within The left side of the bar, a language object.
#' @param env The formula's environment.
#'
#' @return `NULL`, or a list with `term`, the unbuilt [SmoothTerm()], and
#'   `intercept`, a single logical.
#'
#' @keywords internal
.random_smooth_lhs <- function(within, env) {
  pieces <- list()
  flat <- function(e) {
    if (is.call(e) && identical(e[[1L]], as.name("+")) && length(e) == 3L) {
      flat(e[[2L]])
      flat(e[[3L]])
    } else if (is.call(e) && identical(e[[1L]], as.name("-")) &&
               length(e) == 3L) {
      flat(e[[2L]])
      pieces[[length(pieces) + 1L]] <<- call("-", e[[3L]])
    } else {
      pieces[[length(pieces) + 1L]] <<- e
    }
  }
  flat(within)
  is_smooth <- vapply(pieces, function(e) {
    is.call(e) && (identical(e[[1L]], as.name("s")) ||
                     identical(e[[1L]], as.name("te")) ||
                     identical(e[[1L]], quote(modelterms7::s)) ||
                     identical(e[[1L]], quote(modelterms7::te)))
  }, logical(1))
  if (!any(is_smooth)) return(NULL)
  if (sum(is_smooth) > 1L) {
    stop(paste0(
      "the left side of the bar holds more than one smooth. A smooth random\n",
      "  effect has one: write one random() term for each smooth."),
      call. = FALSE)
  }
  intercept <- TRUE
  for (e in pieces[!is_smooth]) {
    if (identical(e, 1) || identical(e, 1L)) next
    if (identical(e, 0) || identical(e, 0L) ||
        identical(e, quote(-1)) || identical(e, call("-", 1))) {
      intercept <- FALSE
      next
    }
    stop(sprintf(paste0(
      "the left side of the bar holds '%s' beside a smooth. A smooth random\n",
      "  effect holds only the smooth, with 1 or 0 for the intercept; write\n",
      "  another covariate as random(~ %s | g), a term of its own."),
      deparse1(e), deparse1(e)), call. = FALSE)
  }
  cl <- pieces[[which(is_smooth)]]
  fun <- if (identical(cl[[1L]], as.name("te")) ||
             identical(cl[[1L]], quote(modelterms7::te))) te else s
  cl[[1L]] <- fun
  sm <- eval(cl, env)
  if (!is.null(sm@by)) {
    stop(paste0(
      "the smooth on the left side of the bar carries a 'by' variable. The\n",
      "  grouping of a random effect is the variable after the bar, so the\n",
      "  smooth is written without one."), call. = FALSE)
  }
  if (!is.null(sm@sparse)) {
    stop(paste0(
      "the smooth on the left side of the bar carries 'sparse'. The block of\n",
      "  a random effect is sparse by construction."), call. = FALSE)
  }
  if (length(sm@hyper) || length(sm@ids)) {
    stop(paste0(
      "the smooth on the left side of the bar carries 'hyper' or 'id'. The\n",
      "  hyperparameters of a smooth random effect belong to the random\n",
      "  effect: name them in the 'hyper' and 'id' arguments of random()."),
      call. = FALSE)
  }
  list(term = sm, intercept = intercept)
}


#' The Component Matrices of a Built Smooth
#'
#' @description
#' Returns the roughness matrices of a built [SmoothTerm()] as a list, one per
#' smoothing parameter: one matrix for a single penalty, several for an
#' adaptive smoother or an anisotropic tensor product.
#'
#' @details
#' A smoother carrying a penalty factory has no roughness matrix, because the
#' factory declares its penalty over the penalized coordinates, and the
#' function signals an error for such a smoother.
#'
#' @param sb A built [SmoothTerm()] without a `by` variable.
#'
#' @return A list of symmetric numeric matrices of side [term_npar()].
#'
#' @keywords internal
.random_smooth_components <- function(sb) {
  pen <- sb@penalty
  if (S7::S7_inherits(pen, penalties7::QuadraticPenalty)) {
    return(list(as.matrix(pen@P)))
  }
  if (S7::S7_inherits(pen, penalties7::AdditivePenalty)) {
    return(lapply(pen@matrices, as.matrix))
  }
  stop(paste0(
    "the smooth on the left side of the bar carries a penalty of its own (a\n",
    "  penalty factory), which a smooth random effect does not support. Use\n",
    "  the smoother's roughness penalty."), call. = FALSE)
}


#' The Null Space and the Penalized Coordinates of a Smooth
#'
#' @description
#' Splits the coefficients of a smooth into the directions on which every
#' component of its penalty vanishes and the directions that the penalty
#' covers, and returns the linear map from the smooth's coefficients to those
#' two groups.
#'
#' @details
#' The null space of a sum of positive semidefinite matrices is the
#' intersection of the components' null spaces, so it is read once from the
#' components normalized to a unit largest entry and does not move with the
#' smoothing parameters.
#'
#' Where the smoother has already placed the null space in its leading columns
#' with a zero border in every component, as every \pkg{basis7}
#' smoother does under `null_space = "keep"`, the map keeps those columns and
#' their names. Otherwise, as for a tensor product, an orthonormal basis of
#' the null space and one of its complement are taken from the eigenvectors
#' of the normalized sum, and the columns are named `null1`, ... and `z1`, ....
#'
#' With a single component the penalized coordinates are then whitened: with
#' \eqn{P = R^\top R} the penalty on the penalized coordinates
#' \eqn{\gamma}, the coordinates \eqn{\delta = R\gamma} carry the identity
#' penalty, so the random effect over them is independent and Gaussian with
#' one standard deviation. A Demmler-Reinsch smoother already has the identity
#' there, and its coordinates are unchanged. With several components the
#' penalized coordinates are kept and each component is carried onto them.
#'
#' @param comps The component matrices, as [.random_smooth_components()]
#'   returns them.
#' @param names The smooth's coordinate names, without its label.
#' @param tol The relative eigenvalue tolerance for a zero.
#'
#' @return A list with `T`, the map (\eqn{p} rows, one column per coordinate
#'   of the random effect, null space first); `nnull`, the dimension of the
#'   null space; `npen`, the number of penalized coordinates; `comps`, the
#'   components on the penalized coordinates; and `names`.
#'
#' @keywords internal
.random_smooth_split <- function(comps, names, tol = 1e-10) {
  p <- nrow(comps[[1L]])
  stacked <- Reduce(`+`, lapply(comps, function(P) P / max(abs(P))))
  es <- eigen(stacked, symmetric = TRUE)
  keep <- es$values > tol * max(es$values)
  r <- sum(keep)
  k <- p - r
  border <- k == 0L || all(vapply(comps, function(P)
    all(P[, seq_len(k)] == 0), logical(1)))
  if (border) {
    N <- diag(1, p)[, seq_len(k), drop = FALSE]
    U <- diag(1, p)[, k + seq_len(r), drop = FALSE]
    nn <- names[seq_len(k)]
    pn <- names[k + seq_len(r)]
  } else {
    N <- es$vectors[, !keep, drop = FALSE]
    U <- es$vectors[, keep, drop = FALSE]
    nn <- if (k) paste0("null", seq_len(k)) else character(0)
    pn <- paste0("z", seq_len(r))
  }
  pc <- lapply(comps, function(P) {
    M <- crossprod(U, P %*% U)
    (M + t(M)) / 2
  })
  if (length(pc) == 1L) {
    R <- chol(pc[[1L]])
    U <- U %*% backsolve(R, diag(1, r))
    pc <- list(diag(1, r))
  }
  list(T = cbind(N, U), nnull = k, npen = r, comps = pc,
       names = c(nn, pn))
}


#' @title Build a Smooth Random Effect
#' @name term_build.RandomSmoothTerm
#'
#' @description
#' Builds the smooth on the left side of the bar, separates its null space
#' from its penalized coordinates, interacts the resulting within-group design
#' with the group indicators, and attaches one prior to the null space and one
#' to the penalized coordinates.
#'
#' @details
#' # The within-group design
#'
#' The design that the effect of one group multiplies is the constant (unless
#' the left side carries `0`), then the null space of the smooth's penalty
#' without the constant, which every smoother removes, then the penalized
#' coordinates, as
#' [.random_smooth_split()] defines them. The null space is whatever the
#' smoother's operator leaves unpenalized: a straight line for a second
#' derivative, a sine and a cosine for the harmonic operator, an exponential
#' for \eqn{D^2 + aD}, nothing for a cyclic smoother.
#'
#' The argument `null_space` of the smoother determines what the null space
#' contributes. Under `"keep"` its columns are there and carry the null-space
#' prior together with the constant; under `"shrink"` they are penalized with
#' the rest and belong to the penalized coordinates; under `"drop"` they are
#' absent.
#'
#' # The two priors
#'
#' The constant and the null-space columns, \eqn{d_0} of them, have the prior
#' that [random()] gives \eqn{d_0} within-group columns: a centered
#' multivariate Gaussian on an unstructured covariance, or one standard
#' deviation per column under `correlated = FALSE`, or `distrib` where it is
#' given. A full covariance makes the model invariant to the origin and the
#' scale of the covariate, and to the phase of a periodic one.
#'
#' The penalized coordinates of every group have the prior of the smooth's own
#' penalty. With one component that is an independent centered Gaussian with
#' one standard deviation, because the coordinates have been whitened, and its
#' hyperparameter is `sigma.smooth`. With several components it is the sum
#' \eqn{\sum_k \lambda_k (I_m \otimes S_k)}, with hyperparameters
#' `lambda1.smooth`, `lambda2.smooth` and so on. The groups share these
#' hyperparameters, so all the curves have the same roughness.
#'
#' @param term An unbuilt or built [RandomSmoothTerm()].
#' @param data A data frame carrying the grouping variable and the smooth's
#'   covariates.
#' @param ... Unused.
#'
#' @return The term with `X` (a `dgCMatrix`), `coef_names`, `build_info` and
#'   `penalty` filled.
#'
#' @seealso [random()], [term_build.RandomTerm()], [term_penalties()].
#'
#' @examples
#' set.seed(1)
#' dd <- data.frame(x = runif(40), g = factor(rep(letters[1:4], 10)))
#' b <- term_build(random(~ s(x, smoother = basis7::bspline_smooth(k = 5)) | g),
#'                 dd)
#' term_group(b)$names
#' vapply(term_penalties(b), function(e) paste(e$penalty@params, collapse = ","),
#'        character(1))
#'
#' @keywords internal
S7::method(term_build, RandomSmoothTerm) <- function(term, data, ...) {
  parts <- .random_parts(term@formula)
  g <- .random_group(parts$group, data)
  m <- nlevels(g)
  if (m < 2L) {
    stop("the grouping variable must have at least two levels.",
         call. = FALSE)
  }

  sb <- term_build(term@smooth, data)
  X <- as.matrix(term_matrix(sb))
  inner <- substring(sb@coef_names, nchar(sb@label) + 2L)
  sp <- .random_smooth_split(.random_smooth_components(sb), inner)
  W <- X %*% sp$T
  wn <- sp$names
  if (isTRUE(term@intercept)) {
    # every smoother removes the constant; a smooth that kept it would make
    # the random level aliased with its own coordinates
    res <- 1 - W %*% qr.solve(W, rep(1, nrow(W)))
    if (max(abs(res)) < 1e-8) {
      stop(paste0(
        "the smooth on the left side of the bar contains the constant, so a\n",
        "  random level beside it is aliased. Write ~ 0 + s(...) | g."),
        call. = FALSE)
    }
    W <- cbind(1, W)
    wn <- c("(Intercept)", wn)
  }
  dimnames(W) <- list(NULL, wn)
  d <- ncol(W)
  d0 <- as.integer(isTRUE(term@intercept)) + sp$nnull
  r <- sp$npen

  Z <- .random_block(g, W)
  cn <- .random_names(term@label, levels(g), wn)
  colnames(Z) <- cn

  base <- (seq_len(m) - 1L) * d
  entries <- list()
  if (d0 > 0L) {
    prior <- term@distrib
    if (is.null(prior)) {
      prior <- .random_default(d0, term@correlated)
    } else if (is.list(prior)) {
      if (length(prior) != d0) {
        stop(sprintf(paste0(
          "'distrib' is a list of %d and the random level and the smooth's\n",
          "  null space have %d columns; a list gives one distribution per\n",
          "  column."), length(prior), d0), call. = FALSE)
      }
      for (j in seq_len(d0)) {
        .random_check_prior(prior[[j]], d0, sprintf("'distrib[[%d]]'", j))
      }
    } else {
      .random_check_prior(prior, d0, "'distrib'")
    }
    pos0 <- as.vector(outer(seq_len(d0), base, `+`))
    entries <- lapply(.random_entries(prior, d0, m, wn[seq_len(d0)]),
                      function(en) {
                        en$index <- pos0[en$index]
                        en
                      })
  } else if (!is.null(term@distrib)) {
    stop(paste0(
      "'distrib' is the prior of the random level and of the smooth's null\n",
      "  space, and this term has neither: the left side carries 0 and the\n",
      "  smoother leaves no unpenalized direction."), call. = FALSE)
  }
  pos1 <- as.vector(outer(d0 + seq_len(r), base, `+`))
  pen1 <- if (length(sp$comps) == 1L) {
    penalties7::distrib_penalty(
      distributions7::fixed(distributions7::gaussian1_distrib(), mu = 0),
      n_coef = m * r)
  } else {
    penalties7::additive_penalty(sp$comps, blocks = m, readable = "sd")
  }
  entries[[length(entries) + 1L]] <- list(
    name = "smooth", index = pos1, penalty = pen1, fixed = list(),
    n_values = list(), values = list(), min_ratio = numeric(0),
    ids = character(0))
  entries <- .entry_hyper(entries, term@hyper, term@ids, term@label)

  term@X <- Z
  term@coef_names <- cn
  term@build_info <- list(
    kind = "smooth", gexpr = parts$group, glevels = levels(g),
    inner = sb, transform = sp$T, intercept = isTRUE(term@intercept),
    entries = entries, dim = d, nnull = d0, wnames = wn
  )
  term@hyper <- do.call(c, c(list(list()), lapply(entries, function(en) {
    if (!length(en$fixed) || !nzchar(en$name)) return(en$fixed)
    stats::setNames(en$fixed, paste0(names(en$fixed), ".", en$name))
  })))
  term@penalty <- NULL
  term
}
