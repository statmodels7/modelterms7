#' @include term_classes.R generics.R
NULL

#' @title S7 Class for Nonlinear Parametric Terms
#' @name NlTerm
#'
#' @description
#' A subclass of \code{\link{additive_term}} for a parametric function
#' that is nonlinear in its own parameters. The design block is the
#' Jacobian of that function in the parameters, so the term is linear in
#' the sense the model layer needs while the function is not; the block
#' depends on where the parameters currently are and is recomputed by
#' \code{\link{term_refresh}}.
#'
#' @inheritParams additive_term
#' @param fn The function or formula defining the contribution.
#' @param nl_params The names of the nonlinear parameters.
#' @param links One link per parameter.
#' @param subformulas One optional formula per parameter.
#' @param deriv_mode How the derivatives are obtained.
#' @param spec The resolved construction settings.
#'
#' @return An object of class \code{NlTerm}.
#'
#' @seealso \code{\link{nl}}
#' @examples
#' S7::S7_inherits(nl(~ theta1 * exp(theta2 * x)), NlTerm)
#' @export
NlTerm <- S7::new_class(
  name = "NlTerm",
  parent = additive_term,
  properties = list(
    fn = S7::class_any,
    nl_params = S7::class_character,
    links = S7::class_list,
    subformulas = S7::class_list,
    deriv_mode = S7::class_character,
    spec = S7::class_list
  )
)

#' Nonlinear Parametric Term
#'
#' @description
#' A contribution \eqn{f(x; \theta)} that is nonlinear in its own
#' parameters, given either as a formula in the covariates and the
#' parameters or as a function of both. Each parameter may carry a link,
#' and, when the term is given as a formula, may itself be modeled with
#' covariates.
#'
#' @details
#' While \eqn{f} is differentiable the term is an ordinary additive one at
#' every point: the contribution is linearized as
#' \deqn{f(x;\theta(\beta)) \approx f(x;\theta(\beta_0))
#'   + J(\beta_0)\,(\beta - \beta_0), \qquad
#'   J = \frac{\partial f}{\partial \beta},}
#' so the design block is the Jacobian, and the only thing that
#' distinguishes the term from a linear one is that the block is refreshed
#' as the parameters move. \code{\link{term_refresh}} does that, and
#' \code{\link{term_value}} reports the contribution itself, which a
#' Gauss-Newton step needs beside the Jacobian.
#'
#' \subsection{Two ways to give the function, with different reach}{
#' A \strong{formula} such as \code{~ theta1 * exp(theta2 * x)} is read
#' symbolically: the names it uses that are not columns of the data are
#' the parameters, and the derivatives come from
#' \code{\link[stats]{deriv}} where that succeeds and from a central
#' difference where it does not. A \strong{function} \code{f(x, theta)},
#' vectorized in both, is treated as opaque: its derivatives are always
#' differenced, and its parameters must be named in \code{params}.
#'
#' The difference is not only in the derivatives. Modeling a parameter
#' with covariates means replacing \eqn{\theta_j} by
#' \eqn{g_j^{-1}(Z\gamma_j)} inside \eqn{f}, which requires knowing where
#' \eqn{\theta_j} enters; a formula says so and an opaque function does
#' not. \code{subformulas} is therefore available on the formula route
#' only, and is rejected on the other.
#' }
#'
#' \subsection{Links and submodels}{
#' \code{links} carries each parameter to an unconstrained scale, so a
#' rate constrained positive is estimated as its logarithm and the
#' optimizer never proposes a negative one. \code{subformulas} develops a
#' parameter as \eqn{\theta_j = g_j^{-1}(Z\gamma_j)} for a design
#' \eqn{Z} built from a one-sided formula, which gives a parameter that
#' varies by group or with a covariate; the coefficients are then the
#' \eqn{\gamma_j}, and the Jacobian carries the chain rule
#' \eqn{\partial f/\partial\theta_j \cdot (g_j^{-1})' \cdot Z}.
#' }
#'
#' @param fn A one-sided formula in the covariates and the parameters, or
#'   a function of \code{(x, theta)} vectorized in both.
#' @param params The parameter names. Required when \code{fn} is a
#'   function; inferred from the formula otherwise.
#' @param x The covariate expression handed to a function \code{fn},
#'   evaluated in the data. Unused for a formula.
#' @param links An optional named list of \pkg{linkfunctions7} links, one
#'   per parameter. Parameters without one carry the identity.
#' @param subformulas An optional named list of one-sided formulas, one
#'   per parameter to be modeled with covariates. Formula input only.
#' @param start An optional named list of starting values for the
#'   parameters, on the parameter scale. Defaults to the inverse link at
#'   zero.
#' @param label A single non-empty string prefixed to the coefficient
#'   names.
#'
#' @return An object of class \code{\link{NlTerm}} (a specification; see
#'   \code{\link{term_build}}).
#'
#' @examples
#' set.seed(1)
#' dd <- data.frame(x = seq(0, 3, length.out = 60))
#' dd$y <- 2 * exp(-1.3 * dd$x) + rnorm(60, sd = 0.05)
#'
#' # an exponential decay, with the rate held positive by a log link
#' spec <- nl(~ a * exp(-r * x),
#'            links = list(r = linkfunctions7::log_link()),
#'            start = list(a = 1, r = 1))
#' built <- term_build(spec, dd)
#' term_coef_names(built)
#' dim(term_matrix(built))
#'
#' @export
nl <- function(fn, params = NULL, x = NULL, links = NULL,
               subformulas = NULL, start = NULL, label = "nl") {
  xe <- substitute(x)
  is_formula <- inherits(fn, "formula")
  if (!is_formula && !is.function(fn)) {
    stop("'fn' must be a one-sided formula or a function of (x, theta).",
         call. = FALSE)
  }
  if (is_formula && length(fn) != 2L) {
    stop("a formula 'fn' must be one-sided, e.g. ~ a * exp(-r * x).",
         call. = FALSE)
  }
  if (!is_formula) {
    if (is.null(params) || !is.character(params) || !length(params)) {
      stop(paste("'params' must name the parameters when 'fn' is a",
                 "function: an opaque function does not say what they are."),
           call. = FALSE)
    }
    if (!is.null(subformulas)) {
      stop(paste("'subformulas' needs to know where a parameter enters the",
                 "model, which only a formula says; give 'fn' as a formula",
                 "to develop a parameter with covariates."), call. = FALSE)
    }
  }
  for (nmv in list(list(links, "links"), list(subformulas, "subformulas"),
                   list(start, "start"))) {
    v <- nmv[[1L]]
    if (!is.null(v) && (!is.list(v) || is.null(names(v)) ||
                        anyNA(names(v)) || !all(nzchar(names(v))))) {
      stop(sprintf("'%s' must be a named list.", nmv[[2L]]), call. = FALSE)
    }
  }
  if (!is.character(label) || length(label) != 1L || is.na(label) ||
      !nzchar(label)) {
    stop("'label' must be a single non-empty character string.",
         call. = FALSE)
  }
  NlTerm(label = label, fn = fn,
         nl_params = if (is_formula) character(0) else params,
         links = if (is.null(links)) list() else links,
         subformulas = if (is.null(subformulas)) list() else subformulas,
         deriv_mode = if (is_formula) "symbolic" else "numeric",
         spec = list(x = xe, start = start, is_formula = is_formula),
         X = NULL, coef_names = character(0),
         blueprint = list(), penalty = NULL)
}

# the parameters of a formula: every name it uses that the data do not
# supply and that is not a function being called
.nl_formula_params <- function(f, data) {
  e <- f[[2L]]
  nms <- all.vars(e)
  fns <- all.names(e)
  called <- setdiff(fns, nms)
  out <- setdiff(nms, c(names(data), called))
  if (!length(out)) {
    stop(paste("the formula uses no name the data do not supply, so it has",
               "no parameters to estimate."), call. = FALSE)
  }
  out
}

.nl_link <- function(links, p) {
  if (!is.null(links[[p]])) links[[p]] else linkfunctions7::identity_link()
}

# theta and d theta / d eta at the current coefficients, per parameter
.nl_theta <- function(bp, coef) {
  out <- list(value = list(), deriv = list())
  for (p in bp$params) {
    idx <- bp$index[[p]]
    Z <- bp$Z[[p]]
    eta <- if (is.null(Z)) coef[idx] else as.numeric(Z %*% coef[idx])
    lk <- bp$links[[p]]
    out$value[[p]] <- linkfunctions7::linkinv(lk, eta)
    out$deriv[[p]] <- linkfunctions7::dlinkinv(lk, eta)
  }
  out
}

# f and its derivatives in the parameters, at the current theta
.nl_eval <- function(bp, th) {
  env <- c(bp$data_vars, th)
  if (identical(bp$mode, "symbolic")) {
    res <- eval(bp$dexpr, env)
    g <- attr(res, "gradient")
    list(value = as.numeric(res) + 0 * bp$one,
         grad = stats::setNames(
           lapply(bp$params, function(p) g[, p] + 0 * bp$one), bp$params))
  } else if (identical(bp$mode, "formula_fd")) {
    val <- as.numeric(eval(bp$expr, env)) + 0 * bp$one
    list(value = val,
         grad = stats::setNames(lapply(bp$params, function(p) {
           h <- .nl_step(th[[p]])
           tp <- th; tm <- th
           tp[[p]] <- th[[p]] + h
           tm[[p]] <- th[[p]] - h
           (as.numeric(eval(bp$expr, c(bp$data_vars, tp))) -
              as.numeric(eval(bp$expr, c(bp$data_vars, tm)))) / (2 * h)
         }), bp$params))
  } else {
    val <- as.numeric(bp$fn(bp$xval, th)) + 0 * bp$one
    list(value = val,
         grad = stats::setNames(lapply(bp$params, function(p) {
           h <- .nl_step(th[[p]])
           tp <- th; tm <- th
           tp[[p]] <- th[[p]] + h
           tm[[p]] <- th[[p]] - h
           (as.numeric(bp$fn(bp$xval, tp)) -
              as.numeric(bp$fn(bp$xval, tm))) / (2 * h)
         }), bp$params))
  }
}

.nl_step <- function(v) .Machine$double.eps^(1 / 3) * pmax(1, abs(v))

# the Jacobian in the coefficients and the value, at one coefficient vector
.nl_jacobian <- function(bp, coef) {
  th <- .nl_theta(bp, coef)
  fv <- .nl_eval(bp, th$value)
  J <- matrix(0, bp$n, bp$ncoef)
  for (p in bp$params) {
    idx <- bp$index[[p]]
    w <- fv$grad[[p]] * th$deriv[[p]]
    Z <- bp$Z[[p]]
    J[, idx] <- if (is.null(Z)) w else w * Z
  }
  list(J = J, value = fv$value)
}

S7::method(term_build, NlTerm) <- function(term, data, ...) {
  n <- nrow(data)
  is_f <- term@spec$is_formula
  params <- if (is_f) .nl_formula_params(term@fn, data) else term@nl_params

  bad <- setdiff(names(term@subformulas), params)
  if (length(bad)) {
    stop(sprintf("'subformulas' names '%s', which is not a parameter.",
                 bad[1L]), call. = FALSE)
  }
  bad <- setdiff(names(term@links), params)
  if (length(bad)) {
    stop(sprintf("'links' names '%s', which is not a parameter.", bad[1L]),
         call. = FALSE)
  }

  links <- stats::setNames(lapply(params, .nl_link, links = term@links),
                           params)

  # the submodel design of each parameter, and where its coefficients sit
  Z <- stats::setNames(vector("list", length(params)), params)
  sub_bp <- stats::setNames(vector("list", length(params)), params)
  index <- stats::setNames(vector("list", length(params)), params)
  cn <- character(0)
  pos <- 0L
  for (p in params) {
    sf <- term@subformulas[[p]]
    if (is.null(sf)) {
      pos <- pos + 1L
      index[[p]] <- pos
      cn <- c(cn, p)
    } else {
      mf <- stats::model.frame(sf, data, na.action = stats::na.pass,
                               drop.unused.levels = FALSE)
      tt <- attr(mf, "terms")
      Zp <- stats::model.matrix(tt, mf)
      # the submodel needs the same blueprint discipline as any other
      # design: its levels and contrasts are recorded here so that
      # prediction reapplies them instead of re-deriving them
      sub_bp[[p]] <- list(terms = stats::delete.response(tt),
                          xlev = stats::.getXlevels(tt, mf),
                          contrasts = attr(Zp, "contrasts"))
      attr(Zp, "assign") <- NULL
      attr(Zp, "contrasts") <- NULL
      Z[[p]] <- Zp
      index[[p]] <- pos + seq_len(ncol(Zp))
      pos <- pos + ncol(Zp)
      cn <- c(cn, paste(p, colnames(Zp), sep = "."))
    }
  }

  # the derivative route: symbolic where deriv() manages the expression,
  # a central difference where it does not, and always a difference for a
  # function the term cannot look inside
  mode <- "numeric"
  dexpr <- NULL
  if (is_f) {
    dexpr <- tryCatch(stats::deriv(term@fn[[2L]], params),
                      error = function(e) NULL)
    mode <- if (is.null(dexpr)) "formula_fd" else "symbolic"
  }

  data_vars <- as.list(data)
  xval <- if (!is_f && !is.null(term@spec$x)) {
    eval(term@spec$x, data, baseenv())
  } else if (!is_f) {
    data
  } else NULL

  bp <- list(params = params, links = links, Z = Z, index = index,
             n = n, ncoef = pos, mode = mode, dexpr = dexpr,
             expr = if (is_f) term@fn[[2L]] else NULL,
             fn = if (is_f) NULL else term@fn, xval = xval,
             data_vars = data_vars, one = numeric(n),
             subformulas = term@subformulas, sub_bp = sub_bp,
             is_formula = is_f)

  # the coefficients the block is built at: the link of the starting value,
  # or zero, which is the inverse link's own natural point
  start <- term@spec$start
  coef0 <- numeric(pos)
  for (p in params) {
    if (!is.null(start[[p]])) {
      coef0[index[[p]][1L]] <- linkfunctions7::linkfun(links[[p]],
                                                       start[[p]])
    }
  }

  jv <- .nl_jacobian(bp, coef0)
  bp$coef <- coef0
  bp$value <- jv$value

  cn <- paste(term@label, cn, sep = ".")
  X <- jv$J
  colnames(X) <- cn
  term@X <- X
  term@coef_names <- cn
  term@nl_params <- params
  term@links <- links
  term@deriv_mode <- mode
  term@blueprint <- bp
  term
}

S7::method(term_predict, NlTerm) <- function(term, newdata, ...) {
  .assert_built(term)
  bp <- term@blueprint
  nb <- bp
  nb$n <- nrow(newdata)
  nb$one <- numeric(nrow(newdata))
  nb$data_vars <- as.list(newdata)
  for (p in bp$params) {
    if (!is.null(bp$Z[[p]])) {
      sb <- bp$sub_bp[[p]]
      mf <- stats::model.frame(sb$terms, newdata, na.action = stats::na.pass,
                               xlev = sb$xlev)
      Zp <- stats::model.matrix(sb$terms, mf, contrasts.arg = sb$contrasts)
      attr(Zp, "assign") <- NULL
      attr(Zp, "contrasts") <- NULL
      nb$Z[[p]] <- Zp
    }
  }
  if (!bp$is_formula) {
    nb$xval <- if (!is.null(term@spec$x)) {
      eval(term@spec$x, newdata, baseenv())
    } else newdata
  }
  X <- .nl_jacobian(nb, bp$coef)$J
  colnames(X) <- term@coef_names
  X
}

#' @title Refresh a Term at New Coefficients
#'
#' @description
#' Recomputes whatever a term's design block depends on when the
#' coefficients move. For every ordinary term this is the identity, the
#' block being a function of the data alone; for a nonlinear term the
#' block is the Jacobian of its contribution, which is a function of where
#' the parameters currently are.
#'
#' @param term A built term.
#' @param coef The current coefficients of the term's block.
#' @param ... Passed to methods.
#'
#' @return A built term, refreshed.
#'
#' @seealso \code{\link{nl}}, \code{\link{term_value}}
#'
#' @examples
#' dd <- data.frame(x = seq(0, 2, length.out = 20))
#' built <- term_build(nl(~ a * exp(-r * x), start = list(a = 1, r = 1)), dd)
#' r1 <- term_refresh(built, c(2, 0.5))
#' max(abs(term_matrix(r1) - term_matrix(built))) > 0
#'
#' @export
term_refresh <- S7::new_generic("term_refresh", "term",
  function(term, coef, ...) S7::S7_dispatch())

S7::method(term_refresh, model_term) <- function(term, coef, ...) term

S7::method(term_refresh, NlTerm) <- function(term, coef, ...) {
  .assert_built(term)
  bp <- term@blueprint
  coef <- as.numeric(coef)
  if (length(coef) != bp$ncoef) {
    stop(sprintf("'coef' must have length %d.", bp$ncoef), call. = FALSE)
  }
  jv <- .nl_jacobian(bp, coef)
  bp$coef <- coef
  bp$value <- jv$value
  X <- jv$J
  colnames(X) <- term@coef_names
  term@X <- X
  term@blueprint <- bp
  term
}

#' @title The Contribution of a Term at Its Current Coefficients
#'
#' @description
#' The values a term contributes to the linear predictor. For a linear
#' term this is the block times the coefficients and carries no
#' information the block does not; for a nonlinear one it is
#' \eqn{f(x;\theta)}, which the Jacobian alone does not give, and which a
#' Gauss-Newton step needs beside it.
#'
#' @param term A built term.
#' @param coef The coefficients. Optional for a nonlinear term, which
#'   carries the ones it was last refreshed at.
#' @param ... Passed to methods.
#'
#' @return A numeric vector, one value per observation.
#'
#' @seealso \code{\link{term_refresh}}
#'
#' @examples
#' dd <- data.frame(x = seq(0, 2, length.out = 20))
#' built <- term_build(nl(~ a * exp(-r * x), start = list(a = 2, r = 1)), dd)
#' head(term_value(built), 3)
#'
#' @export
term_value <- S7::new_generic("term_value", "term",
  function(term, coef = NULL, ...) S7::S7_dispatch())

S7::method(term_value, additive_term) <- function(term, coef = NULL, ...) {
  .assert_built(term)
  if (is.null(coef)) {
    stop("'coef' is required: a linear term contributes the block times it.",
         call. = FALSE)
  }
  as.numeric(term@X %*% coef)
}

S7::method(term_value, NlTerm) <- function(term, coef = NULL, ...) {
  .assert_built(term)
  if (is.null(coef)) return(term@blueprint$value)
  .nl_jacobian(term@blueprint, as.numeric(coef))$value
}

S7::method(print, NlTerm) <- function(x, ...) {
  if (term_is_built(x)) {
    cat(sprintf("<NlTerm> '%s' built: %d coefficient%s; %s derivatives\n",
                x@label, ncol(x@X), if (ncol(x@X) == 1L) "" else "s",
                x@deriv_mode))
    cat("  parameters: ", paste(x@nl_params, collapse = ", "), "\n", sep = "")
  } else {
    cat(sprintf("<NlTerm> '%s' (specification; call term_build() with data)\n",
                x@label))
  }
  invisible(x)
}
