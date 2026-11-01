#' @include penalized.R penalized_ctors.R
NULL

#' @title S7 Class for the Jeffreys Penalty
#' @name JeffreysPenalty
#'
#' @description
#' `JeffreysPenalty` is the subclass of [penalties7::penalty()] that a [jfp()]
#' term attaches to its block. It carries the power of the determinant and the
#' setting of the constant direction, and it has no hyperparameters.
#'
#' @details
#' The penalty is \eqn{-a\log|\tilde I(\beta)|}, where \eqn{\tilde I} is the
#' expected Fisher information of the block's coefficients, augmented by the
#' constant direction of the formula level where [jfp()] includes it. That
#' matrix depends on the family and on every predictor of the model, so the
#' penalty cannot be evaluated from the block's coefficients alone, and the
#' methods of [penalties7::penalty_value()], [penalties7::penalty_gradient()],
#' [penalties7::penalty_hessian()] and of the other generics that read the
#' coefficients signal an error. The fitting layer evaluates it on the model.
#'
#' The methods that do not need the model return their answer:
#' [penalties7::penalty_kinks()] is empty, [penalties7::is_proper()] and
#' [penalties7::is_quadratic()] are `FALSE`.
#'
#' @param power The power \eqn{a} of the determinant, a single positive
#'   number.
#' @param intercept `NULL`, `TRUE` or `FALSE`, the setting of the constant
#'   direction (see [jfp()]).
#' @param n_coef The number of coefficients in the block.
#' @param penalty_name,map,params,params_bounds,link_params,params_smooth
#'   Properties inherited from [penalties7::penalty()]. The penalty has no map
#'   and no hyperparameters, so `jeffreys_penalty()` sets the name
#'   `"jeffreys"`, `map = NULL`, an empty `params`, empty lists of bounds and
#'   links, and an empty logical `params_smooth`.
#'
#' @return `jeffreys_penalty()` returns an S7 object of class
#'   `JeffreysPenalty`, inheriting from [penalties7::penalty()].
#'
#' @seealso [jfp()], the term that attaches it.
#'
#' @examples
#' p <- jeffreys_penalty(3)
#' p@power
#' p@params
#' penalties7::penalty_kinks(p, list())
#' penalties7::is_proper(p)
#'
#' @export
JeffreysPenalty <- S7::new_class(
  name = "JeffreysPenalty",
  parent = penalties7::penalty,
  properties = list(
    power = S7::class_numeric,
    intercept = S7::class_any
  )
)

#' @rdname JeffreysPenalty
#' @export
jeffreys_penalty <- function(n_coef, power = 0.5, intercept = NULL) {
  .check_jfp_power(power)
  .check_jfp_intercept(intercept)
  JeffreysPenalty(penalty_name = "jeffreys", map = NULL,
                  n_coef = as.numeric(n_coef), params = character(0),
                  params_bounds = list(), link_params = list(),
                  params_smooth = logical(0), power = power,
                  intercept = intercept)
}

.check_jfp_power <- function(power) {
  if (!is.numeric(power) || length(power) != 1L || !is.finite(power) ||
      power <= 0) {
    stop("'power' must be a single positive number.", call. = FALSE)
  }
  invisible(NULL)
}

.check_jfp_intercept <- function(intercept) {
  if (!is.null(intercept) &&
      (!is.logical(intercept) || length(intercept) != 1L ||
       is.na(intercept))) {
    stop("'intercept' must be NULL, TRUE or FALSE.", call. = FALSE)
  }
  invisible(NULL)
}

.jfp_stop <- function() {
  stop(paste0("The Jeffreys penalty depends on the family and on every ",
              "predictor of the model,\n  so it is evaluated by the fitting ",
              "layer and not from the coefficients alone."),
       call. = FALSE)
}

# A method with the formals of the generic it is registered on, so that S7
# finds every argument it dispatches and passes.
.jfp_needs_model <- function(g) {
  f <- function(penalty, ...) .jfp_stop()
  formals(f) <- formals(utils::getFromNamespace(g, "penalties7"))
  f
}

# The penalties7 generics, declared as external generics so that the
# methods below are registered by S7::methods_register() when penalties7 is
# loaded.
.ext_penalty_value <- S7::new_external_generic("penalties7", "penalty_value", "penalty")
.ext_penalty_gradient <- S7::new_external_generic("penalties7", "penalty_gradient", "penalty")
.ext_penalty_hessian <- S7::new_external_generic("penalties7", "penalty_hessian", "penalty")
.ext_penalty_grad_theta <- S7::new_external_generic("penalties7", "penalty_grad_theta", "penalty")
.ext_penalty_hess_theta <- S7::new_external_generic("penalties7", "penalty_hess_theta", "penalty")
.ext_penalty_cross <- S7::new_external_generic("penalties7", "penalty_cross", "penalty")
.ext_penalty_dhessian <- S7::new_external_generic("penalties7", "penalty_dhessian", "penalty")
.ext_penalty_dcross <- S7::new_external_generic("penalties7", "penalty_dcross", "penalty")
.ext_penalty_dhessian_beta <- S7::new_external_generic("penalties7", "penalty_dhessian_beta", "penalty")
.ext_penalty_dhessian_beta_theta <- S7::new_external_generic("penalties7", "penalty_dhessian_beta_theta", "penalty")
.ext_penalty_draw <- S7::new_external_generic("penalties7", "penalty_draw", "penalty")
.ext_penalty_prox <- S7::new_external_generic("penalties7", "penalty_prox", "penalty")
.ext_penalty_prox_spec <- S7::new_external_generic("penalties7", "penalty_prox_spec", "penalty")
.ext_penalty_kinks <- S7::new_external_generic("penalties7", "penalty_kinks", "penalty")
.ext_is_proper <- S7::new_external_generic("penalties7", "is_proper", "penalty")
.ext_is_quadratic <- S7::new_external_generic("penalties7", "is_quadratic", "penalty")

S7::method(.ext_penalty_value, JeffreysPenalty) <-
  .jfp_needs_model("penalty_value")
S7::method(.ext_penalty_gradient, JeffreysPenalty) <-
  .jfp_needs_model("penalty_gradient")
S7::method(.ext_penalty_hessian, JeffreysPenalty) <-
  .jfp_needs_model("penalty_hessian")
S7::method(.ext_penalty_grad_theta, JeffreysPenalty) <-
  .jfp_needs_model("penalty_grad_theta")
S7::method(.ext_penalty_hess_theta, JeffreysPenalty) <-
  .jfp_needs_model("penalty_hess_theta")
S7::method(.ext_penalty_cross, JeffreysPenalty) <-
  .jfp_needs_model("penalty_cross")
S7::method(.ext_penalty_dhessian, JeffreysPenalty) <-
  .jfp_needs_model("penalty_dhessian")
S7::method(.ext_penalty_dcross, JeffreysPenalty) <-
  .jfp_needs_model("penalty_dcross")
S7::method(.ext_penalty_dhessian_beta, JeffreysPenalty) <-
  .jfp_needs_model("penalty_dhessian_beta")
S7::method(.ext_penalty_dhessian_beta_theta, JeffreysPenalty) <-
  .jfp_needs_model("penalty_dhessian_beta_theta")
S7::method(.ext_penalty_draw, JeffreysPenalty) <-
  .jfp_needs_model("penalty_draw")
S7::method(.ext_penalty_prox, JeffreysPenalty) <-
  .jfp_needs_model("penalty_prox")
S7::method(.ext_penalty_prox_spec, JeffreysPenalty) <-
  .jfp_needs_model("penalty_prox_spec")
S7::method(.ext_penalty_kinks, JeffreysPenalty) <-
  function(penalty, theta, ...) numeric(0)
S7::method(.ext_is_proper, JeffreysPenalty) <- function(penalty, ...) FALSE
S7::method(.ext_is_quadratic, JeffreysPenalty) <-
  function(penalty, ...) FALSE


#' @title S7 Class for the Jeffreys Term
#' @name JeffreysTerm
#'
#' @description
#' `JeffreysTerm` is the subclass of [PenalizedTerm()] that [jfp()]
#' constructs. It adds the power of the determinant and the setting of the
#' constant direction; the input handling, the storage and the prediction are
#' those of the penalized terms.
#'
#' @details
#' The build rejects a block whose columns are linearly dependent, by the
#' rank of the pivoted QR decomposition with the tolerance `1e-7` with which
#' \pkg{statmodels7} detects aliased coefficients: the determinant of the
#' block's information would be zero.
#'
#' @inheritParams PenalizedTerm
#' @param power The power of the determinant.
#' @param intercept The setting of the constant direction.
#'
#' @return An S7 object of class `JeffreysTerm`, inheriting from
#'   [PenalizedTerm()].
#'
#' @seealso [jfp()]
#'
#' @examples
#' t <- jfp(~ x)
#' S7::S7_inherits(t, JeffreysTerm)
#' S7::S7_inherits(t, PenalizedTerm)
#'
#' @export
JeffreysTerm <- S7::new_class(
  name = "JeffreysTerm",
  parent = PenalizedTerm,
  properties = list(
    power = S7::class_numeric,
    intercept = S7::class_any
  )
)


#' Jeffreys Prior on a Block of Coefficients
#'
#' @description
#' Specifies a block of coefficients under the prior
#' \eqn{|\tilde I(\beta)|^{a}}, where \eqn{\tilde I} is the expected Fisher
#' information of the block's coefficients. With \eqn{a = 1/2} this is the
#' Jeffreys prior on the block, and over every coefficient of a binomial
#' regression it is the penalty of Firth (1993).
#'
#' @details
#' The penalty added to the negative log-likelihood is
#' \deqn{\rho(\beta) = -a\log|\tilde I(\beta)|,\qquad
#'   \tilde I = \sum_t w_t\sum_{q,r}\mathcal I_{qr}(\eta_t)\,
#'   \tilde v_{q,t}^\top\tilde v_{r,t},}
#' with \eqn{\mathcal I_{qr}} the expected information of the family on the
#' link scale and \eqn{\tilde v_{q,t}} the derivative of the predictor of
#' equation \eqn{q} with respect to the block's coefficients, preceded by the
#' constant direction when it is included. It depends on the family and on
#' every predictor of the model, and the fitting layer evaluates it together
#' with its derivatives.
#'
#' The power is fixed. There is no hyperparameter, so nothing is estimated
#' for this term and it adds no row to an outer criterion.
#'
#' The constant direction of a formula level is the derivative of the
#' predictor along a constant shift of that level's linear predictor: the
#' column of ones for an equation. With `intercept = NULL` it is included
#' when the column of ones lies in the span of the level's unpenalized
#' columns and not in the span of the block, so the penalty is invariant to
#' centering a column of the block whenever the model is. `TRUE` includes it
#' unless the block already spans the column of ones, and `FALSE` never
#' includes it.
#'
#' The penalty moves an estimate as \eqn{a\,\tilde p} pseudo-observations
#' would, with \eqn{\tilde p} the size of \eqn{\tilde I}: with a Poisson log
#' link and an intercept in the model,
#' \eqn{\sum_t(\hat\mu_t - y_t) = a\tilde p}. On a dispersion or shape
#' parameter the effect is strong: the degrees of freedom of a Student t
#' fitted to Gaussian data stay finite.
#'
#' The matrix \eqn{\tilde I} is the expected information as the family gives
#' it, without modification, as in the definition of Jeffreys. The
#' information of a block's equation generally depends on the other
#' parameters of the distribution as well: for a gaussian the information on
#' the mean is divided by \eqn{\sigma^2}, and for a negative binomial it
#' depends on \eqn{\theta}. The prior therefore also acts on the coefficients
#' of those parameters, including the parameters that an outer criterion
#' such as `reml()` estimates. The shift is of order \eqn{1/n}. For a
#' gaussian with the block in the mean, the penalty equals a constant minus
#' \eqn{2a\tilde p\log\sigma}, so the REML estimate of \eqn{\sigma^2} is the
#' residual sum of squares divided by \eqn{n - p + 2a\tilde p} instead of
#' \eqn{n - p}, and \eqn{n\,\Delta\log\hat\sigma} tends to \eqn{-a\tilde p}.
#' In simulations over seven families with the block in the mean, the median
#' shift of a dispersion or shape parameter was 0.1 to 0.3 standard errors at
#' \eqn{n = 50} and 0.02 to 0.06 standard errors at \eqn{n = 1000}. A block
#' in the equation of a parameter other than the mean moves the coefficients
#' of the mean only when the information of that parameter depends on the
#' mean, as for the \eqn{\theta} of a negative binomial or the probability of
#' a zero-inflated family, and not for the \eqn{\sigma} of a gaussian.
#'
#' A formula is coded with its contrasts, as [linpar()] codes it, and the
#' intercept column is then removed, so `jfp(~ g)` over a factor of \eqn{k}
#' levels gives the \eqn{k - 1} contrasts and `jfp(~ 0 + g)` all \eqn{k}
#' indicators. This differs from [ridge()], whose formula gives every
#' indicator: the Jeffreys penalty depends on the predictor alone, so it does
#' not identify the direction in which a full set of indicators and the
#' equation's intercept trade a constant. A block whose columns are linearly
#' dependent is rejected when the term is built.
#'
#' @inheritParams penalized_terms
#' @param label A single non-empty character string prefixed to the
#'   coefficient names as `label.name`, `"jfp"` by default.
#' @param power The power \eqn{a} of the determinant, a single positive
#'   number; `0.5` by default.
#' @param intercept `NULL` (the default), `TRUE` or `FALSE`: whether the
#'   constant direction of the formula level enters the determinant, by the
#'   rule above.
#'
#' @return An unbuilt [JeffreysTerm()]. The build attaches a
#'   [JeffreysPenalty()] over the block's columns.
#'
#' @references
#' Firth, D. (1993). Bias reduction of maximum likelihood estimates.
#' *Biometrika* 80, 27--38.
#'
#' Kosmidis, I. and Firth, D. (2021). Jeffreys-prior penalty, finiteness and
#' shrinkage in binomial-response generalized linear models. *Biometrika*
#' 108, 71--82.
#'
#' @seealso [JeffreysPenalty()], [ridge()] for the same input handling.
#'
#' @examples
#' set.seed(3)
#' dd <- data.frame(x1 = rnorm(20), x2 = rnorm(20))
#' b <- term_build(jfp(~ x1 + x2), dd)
#' b@coef_names
#' b@penalty@power
#'
#' # one entry, with no hyperparameter
#' e <- term_penalties(b)
#' e[[1]]$penalty@params
#'
#' @export
jfp <- function(x, label = "jfp", power = 0.5, intercept = NULL,
                sparse = NULL, ...) {
  .check_jfp_power(power)
  .check_jfp_intercept(intercept)
  extra <- list(...)
  if (length(extra)) {
    stop(sprintf(paste0("'%s' has no argument '%s'.\n  It has no ",
                        "hyperparameters: the power is held by 'power'."),
                 label, names(extra)[[1L]]), call. = FALSE)
  }
  pt <- .penalized_spec(x, substitute(x), label, FALSE,
                        function(k, map = NULL) {
                          jeffreys_penalty(k, power = power,
                                           intercept = intercept)
                        },
                        list(), list(), sparse = sparse)
  do.call(JeffreysTerm, c(S7::props(pt),
                          list(power = power, intercept = intercept)),
          quote = TRUE)
}

S7::method(term_build, JeffreysTerm) <- function(term, data, ...) {
  b <- term_build(S7::super(term, PenalizedTerm), data, ...)
  if (!S7::S7_inherits(b, JeffreysTerm)) {
    b <- do.call(JeffreysTerm, c(S7::props(b),
                                 list(power = term@power,
                                      intercept = term@intercept)),
                 quote = TRUE)
  }
  # A FORMULA IS CODED WITH ITS CONTRASTS, as linpar() codes it, and the
  # intercept column is then removed. The Jeffreys penalty depends on the
  # predictor alone, so it is invariant along the direction in which a
  # full set of indicators and the equation's intercept trade a constant:
  # coded as ridge() codes it, jfp(~ g) was aliased beside every intercept.
  if (inherits(term@input, "formula") && !.intercept_only(term@input) &&
      attr(stats::terms(term@input), "intercept") == 1L) {
    mf <- stats::model.frame(term@input, data, na.action = stats::na.pass,
                             drop.unused.levels = FALSE)
    tt <- attr(mf, "terms")
    sp <- isTRUE(b@build_info$sparse)
    d <- .design_matrix(tt, mf, NULL, sp)
    X <- d$X[, colnames(d$X) != "(Intercept)", drop = FALSE]
    if (!ncol(X)) {
      stop(sprintf("the block of '%s' has no columns.", term@label),
           call. = FALSE)
    }
    cn <- paste(term@label, colnames(X), sep = ".")
    colnames(X) <- cn
    rownames(X) <- NULL
    b@X <- X
    b@coef_names <- cn
    b@build_info <- list(kind = "formula", terms = stats::delete.response(tt),
                         xlev = stats::.getXlevels(tt, mf),
                         contrasts = d$contrasts, sparse = sp,
                         contrast_coded = TRUE)
    b@penalty <- term@penalty_fun(ncol(X))
  }
  X <- b@X
  r <- qr(as.matrix(X), tol = 1e-7)$rank
  if (r < ncol(X)) {
    stop(sprintf(paste0("the columns of '%s' are linearly dependent (rank %d ",
                        "of %d), and the\n  determinant of their information ",
                        "would be zero."), term@label, r, ncol(X)),
         call. = FALSE)
  }
  b
}

S7::method(term_predict, JeffreysTerm) <- function(term, newdata, ...) {
  bp <- term@build_info
  if (!isTRUE(bp$contrast_coded)) {
    return(term_predict(S7::super(term, PenalizedTerm), newdata, ...))
  }
  .assert_built(term)
  mf <- stats::model.frame(bp$terms, newdata, na.action = stats::na.pass,
                           xlev = bp$xlev)
  X <- .design_matrix(bp$terms, mf, bp$contrasts, isTRUE(bp$sparse))$X
  X <- X[, colnames(X) != "(Intercept)", drop = FALSE]
  colnames(X) <- term@coef_names
  rownames(X) <- NULL
  X
}
