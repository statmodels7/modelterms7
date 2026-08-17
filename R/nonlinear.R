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
#' optimizer never proposes a negative one. A subformula develops a
#' parameter as \eqn{\theta_j = g_j^{-1}(Z\gamma_j)} for a design
#' \eqn{Z} built from its right-hand side, which gives a parameter that
#' varies by group or with a covariate; the coefficients are then the
#' \eqn{\gamma_j}, and the Jacobian carries the chain rule
#' \eqn{\partial f/\partial\theta_j \cdot (g_j^{-1})' \cdot Z}. Because
#' the development acts on the unconstrained scale, the parameter stays
#' in its own set at every observation whatever the coefficients are: a
#' rate on the log link is positive for every row of \eqn{Z}.
#'
#' A subformula is written in \code{...} as a two-sided formula whose
#' left side names the parameter, \code{theta1 ~ z}, or programmatically
#' as \code{subformulas = list(theta1 = ~z)}; the two spellings are the
#' same and a parameter may carry only one. The right-hand side goes
#' through \code{\link{interpret_formula}}, so it takes any term of this
#' package: \code{theta1 ~ ridge(g)} is a population value (the
#' intercept) plus shrunken departures, \code{theta1 ~ s(z)} lets the
#' parameter move smoothly with a covariate, and
#' \code{theta1 ~ random(~1 | g)} is a random intercept on the
#' parameter's unconstrained scale. The penalties the sub-terms carry are
#' reported through \code{\link{term_penalties}} under the key
#' \code{parameter::subterm}, so a fitting layer estimates their
#' hyperparameters as it does any other term's. A structural term, and a
#' term whose block moves with its coefficients, are rejected: a
#' parameter's submodel must be a fixed design.
#' }
#'
#' \subsection{Penalizing a parameter}{
#' A penalty is asked for inside the subformula, where the sub-term that
#' carries it declares it: \code{theta1 ~ lasso(~z1 + z2)} selects which
#' covariates the parameter depends on, and
#' \code{theta1 ~ ridge(~g)} or \code{theta1 ~ random(~1 | g)} shrinks
#' per-group departures towards a population value (the intercept, which
#' stays unpenalized). Each sub-term brings its own hyperparameter, which
#' is what two parameters of a nonlinear function want, being on scales of
#' their own. Earlier releases carried \code{penalty} and \code{penalize}
#' arguments over the parameters' coefficients; both are gone, the
#' sub-terms covering the cases that matter with the hyperparameters in
#' the right place.
#' }
#'
#' \subsection{Supplying the derivatives}{
#' The derivatives of \eqn{f} in its own parameters drive everything: the
#' first is the design block, and a fitting layer reads the second and third
#' where the block moves with the coefficients. They come from whichever of
#' four routes is available, and the choice is made ONE ORDER AT A TIME:
#' a function given here, then symbolic differentiation of the expression,
#' then one stencil applied to the highest order that IS analytic.
#'
#' \code{gradient}, \code{hessian}, \code{deriv3} and \code{deriv4} are
#' independent, so the orders worth writing out by hand can be written and the
#' rest left alone. \strong{Writing the Hessian pays twice}: the third and
#' fourth orders are then one difference away from an exact second rather than
#' from the function. Measured on \eqn{f = a e^{-rx}} given as an opaque
#' function, so that nothing is symbolic, against the closed forms:
#'
#' \tabular{lrrrr}{
#'   \tab order 1 \tab 2 \tab 3 \tab 4 \cr
#'   nothing supplied \tab 5.2e-13 \tab 8.4e-03 \tab 4.62 \tab 1.96e+03 \cr
#'   gradient and hessian \tab 0 \tab 0 \tab 2.2e-12 \tab 8.9e-11
#' }
#'
#' Each function takes \code{(theta, data)} -- a named list of the parameters,
#' each of length \code{n} or 1, and the variables the formula names -- and
#' returns a named list of numeric vectors of length \code{n}. The names are
#' the components of that order, keyed as \pkg{distributions7} keys its own
#' derivative surfaces: the parameter names joined by \code{"_"}. Order two of
#' \code{c("a", "r")} is therefore \code{a_a}, \code{a_r}, \code{r_r}.
#'
#' The names are normalized here, so \code{r_a} and \code{a_r} are the same
#' component and the order in which they are returned does not matter; what is
#' checked is the SET. A name that is not a component of this term, a missing
#' component or a repeated one is an error at \code{\link{term_build}} rather
#' than a silent fall-through to the numerical route, an exact derivative that
#' is quietly not used being worse than none.
#'
#' The derivatives are in the parameters \eqn{\theta}, not in the
#' coefficients: the chain rule onto the coefficients -- each parameter's link,
#' and a subformula's design -- is the term's, which is the only thing that
#' knows them. \code{\link{nl_fderiv}} reads any order back.
#' }
#'
#' @param fn A one-sided formula in the covariates and the parameters, or
#'   a function of \code{(x, theta)} vectorized in both.
#' @param ... Two-sided formulas whose left side names a parameter, one
#'   per parameter to be modeled with covariates, e.g.
#'   \code{theta1 ~ s(z)}. Formula input only; the same as naming the
#'   right-hand side in \code{subformulas}.
#' @param params The parameter names. Required when \code{fn} is a
#'   function; inferred from the formula otherwise.
#' @param x The covariate expression handed to a function \code{fn},
#'   evaluated in the data. Unused for a formula.
#' @param links An optional named list of \pkg{linkfunctions7} links, one
#'   per parameter. Parameters without one carry the identity.
#' @param subformulas An optional named list of one-sided formulas, the
#'   programmatic spelling of the formulas of \code{...}. Formula input
#'   only.
#' @param start An optional named list of starting values for the
#'   parameters, on the parameter scale. Defaults to the inverse link at
#'   zero.
#' @param gradient,hessian,deriv3,deriv4 Optional functions
#'   \code{function(theta, data)} returning the derivatives of \code{fn} in
#'   its own parameters, each a NAMED LIST of numeric vectors. Supply as many
#'   as are worth writing out and leave the rest; see the section below.
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
#' # the amplitude developed by group: a population value plus shrunken
#' # departures, whose hyperparameter a fitting layer estimates
#' dd$g <- factor(rep(c("u", "v"), 30))
#' sub <- term_build(nl(~ a * exp(-r * x), a ~ ridge(~g),
#'                      start = list(a = 1, r = 1)), dd)
#' vapply(term_penalties(sub), function(e) e$name, character(1))
#'
#' # the derivatives written out: as many orders as are worth the work, the
#' # rest one difference away from the highest that is exact
#' g <- function(theta, data) {
#'   e <- exp(-theta$r * data$x)
#'   list(a = e, r = -theta$a * data$x * e)
#' }
#' h <- function(theta, data) {
#'   e <- exp(-theta$r * data$x)
#'   list(a_a = 0 * e, r_r = theta$a * data$x^2 * e, a_r = -data$x * e)
#' }
#' ex <- term_build(nl(~ a * exp(-r * x), gradient = g, hessian = h,
#'                     start = list(a = 1, r = 1)), dd)
#' names(nl_fderiv(ex, order = 3))
#'
#' @seealso \code{\link{s}}, \code{\link{te}}, \code{\link{random}}
#' @export
nl <- function(fn, ..., params = NULL, x = NULL, links = NULL,
               subformulas = NULL, start = NULL,
               gradient = NULL, hessian = NULL, deriv3 = NULL, deriv4 = NULL,
               label = "nl") {
  xe <- substitute(x)
  subformulas <- .nl_gather_subformulas(list(...), subformulas)
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
  supplied <- list(gradient = gradient, hessian = hessian,
                   deriv3 = deriv3, deriv4 = deriv4)
  for (nm in names(supplied)) {
    f <- supplied[[nm]]
    if (is.null(f)) next
    if (!is.function(f)) {
      stop(sprintf("'%s' must be a function of (theta, data).", nm),
           call. = FALSE)
    }
    if (length(formals(f)) < 2L) {
      stop(sprintf(paste("'%s' must take (theta, data): a named list of the",
                         "parameters and\n  the data the formula names."), nm),
           call. = FALSE)
    }
  }
  supplied <- supplied[!vapply(supplied, is.null, logical(1))]

  NlTerm(label = label, fn = fn,
         nl_params = if (is_formula) character(0) else params,
         links = if (is.null(links)) list() else links,
         subformulas = if (is.null(subformulas)) list() else subformulas,
         deriv_mode = if (is_formula) "symbolic" else "numeric",
         spec = list(x = xe, start = start, is_formula = is_formula,
                     supplied = supplied),
         X = NULL, coef_names = character(0),
         blueprint = list(), penalty = NULL)
}

# The two spellings of a parameter's subformula, merged: a two-sided formula
# in `...` whose left side names the parameter, and an entry of the
# `subformulas` list. A parameter may carry only one.
.nl_gather_subformulas <- function(dots, subformulas) {
  if (!length(dots)) return(subformulas)
  # a named entry is an argument the signature does not take -- the shape a
  # mistyped or removed argument arrives in -- and it is reported by name
  nms <- names(dots)
  if (!is.null(nms) && any(nzchar(nms))) {
    stop(sprintf(paste("unused argument '%s': '...' takes two-sided",
                       "formulas whose left side names a parameter."),
                 nms[nzchar(nms)][1L]), call. = FALSE)
  }
  if (!is.null(subformulas) &&
      (!is.list(subformulas) || is.null(names(subformulas)) ||
       anyNA(names(subformulas)) || !all(nzchar(names(subformulas))))) {
    stop("'subformulas' must be a named list.", call. = FALSE)
  }
  out <- if (is.null(subformulas)) list() else subformulas
  for (d in dots) {
    if (!inherits(d, "formula") || length(d) != 3L || !is.name(d[[2L]])) {
      stop(paste("arguments in '...' must be two-sided formulas whose left",
                 "side names a parameter, e.g. theta1 ~ s(z)."),
           call. = FALSE)
    }
    pn <- as.character(d[[2L]])
    if (!is.null(out[[pn]])) {
      stop(sprintf("parameter '%s' carries two subformulas.", pn),
           call. = FALSE)
    }
    out[[pn]] <- stats::as.formula(call("~", d[[3L]]), env = environment(d))
  }
  out
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
    out$deriv2[[p]] <- linkfunctions7::d2linkinv(lk, eta)
  }
  out
}

# f and its derivatives in the parameters, at the current theta
.nl_eval <- function(bp, th) {
  env <- c(bp$data_vars, th)
  # a gradient written out by the caller is the block, and is preferred to
  # every other route: it is what makes the Jacobian -- and so every step the
  # fit takes -- exact rather than a difference
  sup <- bp$supplied$gradient
  if (!is.null(sup)) {
    val <- if (bp$is_formula) as.numeric(eval(bp$expr, env)) else
      as.numeric(bp$fn(bp$xval, th))
    return(list(value = val + 0 * bp$one,
                grad = .nl_user_order(sup, bp$params, 1L, th, bp$data_vars,
                                      bp$one)))
  }
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
           .nl_fd(function(v) {
             tv <- th
             tv[[p]] <- v
             as.numeric(eval(bp$expr, c(bp$data_vars, tv)))
           }, th[[p]], bp$one)
         }), bp$params))
  } else {
    val <- as.numeric(bp$fn(bp$xval, th)) + 0 * bp$one
    list(value = val,
         grad = stats::setNames(lapply(bp$params, function(p) {
           .nl_fd(function(v) {
             tv <- th
             tv[[p]] <- v
             as.numeric(bp$fn(bp$xval, tv))
           }, th[[p]], bp$one)
         }), bp$params))
  }
}

# The derivative of the contribution in ONE parameter, by a stencil from
# numericals7 rather than a difference written out here: the nodes, the
# weights and the step are that package's, so a five-point rule is a change
# of one argument and the step is never a constant copied from somewhere
# else. The parameter is a scalar or one value per observation while the
# function returns one value per observation either way, which is why the
# stencil is assembled here instead of calling fd_derivative(), whose f maps
# a vector of points to the values at those points.
#
# The accuracy is four -- five nodes -- and not two. It costs two extra
# evaluations of f and buys an order of magnitude, which is the trade
# distributions7 measured for the skew t's degrees of freedom and the same
# one applies here: this derivative IS the design block, so what it costs is
# the accuracy of every step the fit takes.
.nl_fd <- function(f, v, one, accuracy = 4L) {
  s <- numericals7::fd_offsets(1L, accuracy)$central
  w <- numericals7::fd_weights(s, 1L)
  h <- numericals7::fd_step(v, 1L, accuracy)
  acc <- 0 * one
  for (k in seq_along(s)) {
    if (w[k] == 0) next
    acc <- acc + w[k] * f(v + s[k] * h)
  }
  acc / h
}

# --- the derivatives of f in its OWN parameters -------------------------
#
# The parameters are ordered as they APPEAR: all.vars() reads the expression
# left to right and setdiff() keeps that order, so `params` is the canonical
# order and every component name is built from it.
#
# A user supplies whichever orders are worth writing out and leaves the rest;
# what is not supplied comes from the symbolic route where the expression can
# be differentiated and otherwise from ONE stencil applied to the highest
# order that IS analytic -- including a supplied one, which is what makes
# writing the Hessian pay twice over.

# The canonical name of a component, built from the parameter order and never
# parsed back. The multiset of indices is sorted, so d2f/da dr and d2f/dr da
# are one component under one name.
.nl_key <- function(params, idx) paste(params[sort(idx)], collapse = "_")

# Every index multiset of a given order, in the enumeration distributions7
# uses: sorted, and lexicographic in the sorted form.
.nl_index_list <- function(np, order) {
  if (order == 1L) return(lapply(seq_len(np), identity))
  prev <- .nl_index_list(np, order - 1L)
  out <- list()
  for (ix in prev) {
    for (j in seq(from = ix[length(ix)], to = np)) out[[length(out) + 1L]] <- c(ix, j)
  }
  out
}

.nl_keys <- function(params, order) {
  vapply(.nl_index_list(length(params), order),
         function(ix) .nl_key(params, ix), character(1))
}

# Which spellings of a component name are accepted, and what each one means.
#
# A user writes "a_r" or "r_a" for the same mixed second derivative, so the
# names are normalized here -- POST HOC, by us, as the canonical order says.
# The map is BUILT from every permutation rather than obtained by splitting
# the user's string: a parameter whose own name contains an underscore makes
# a name ambiguous to read back, which is the trap this package records for
# Hessian component names, and building cannot fall into it.
.nl_key_map <- function(params, order) {
  keys <- character(0)
  canon <- character(0)
  for (ix in .nl_index_list(length(params), order)) {
    ck <- .nl_key(params, ix)
    for (pm in .nl_permutations(ix)) {
      keys <- c(keys, paste(params[pm], collapse = "_"))
      canon <- c(canon, ck)
    }
  }
  dup <- duplicated(keys)
  if (any(dup)) {
    bad <- unique(keys[dup])
    # two different components spelling the same string: only reachable when a
    # parameter's own name contains an underscore, and it is reported rather
    # than resolved by a guess
    if (any(vapply(bad, function(b) length(unique(canon[keys == b])) > 1L,
                   logical(1)))) {
      stop(sprintf(paste("the parameter names make the component name '%s'",
                         "ambiguous at order %d.\n  Rename a parameter: a name",
                         "containing '_' cannot be told apart from a\n",
                         " combination of others."), bad[1L], order),
           call. = FALSE)
    }
  }
  stats::setNames(canon[!dup], keys[!dup])
}

.nl_permutations <- function(x) {
  n <- length(x)
  if (n <= 1L) return(list(x))
  out <- list()
  for (i in seq_len(n)) {
    for (rest in .nl_permutations(x[-i])) out[[length(out) + 1L]] <- c(x[i], rest)
  }
  out
}

# One order from the user, with the names normalized and checked. A wrong name
# is an ERROR and not a silent fall-through to the numerical route: a supplied
# derivative that is quietly not used is an exact quantity thrown away.
.nl_user_order <- function(f, params, order, th, data_vars, one) {
  out <- f(th, data_vars)
  if (!is.list(out) || is.null(names(out)) || anyNA(names(out)) ||
      !all(nzchar(names(out)))) {
    stop(sprintf(paste("the '%s' function must return a NAMED list, one entry",
                       "per component."), .nl_arg_name(order)), call. = FALSE)
  }
  map <- .nl_key_map(params, order)
  want <- .nl_keys(params, order)
  unknown <- setdiff(names(out), names(map))
  if (length(unknown)) {
    stop(sprintf(paste("'%s' returned the component '%s', which is not one of",
                       "this term's.\n  Expected: %s."),
                 .nl_arg_name(order), unknown[1L],
                 paste(want, collapse = ", ")), call. = FALSE)
  }
  names(out) <- unname(map[names(out)])
  if (anyDuplicated(names(out))) {
    stop(sprintf("'%s' returned the same component twice ('%s').",
                 .nl_arg_name(order), names(out)[anyDuplicated(names(out))]),
         call. = FALSE)
  }
  missing <- setdiff(want, names(out))
  if (length(missing)) {
    stop(sprintf(paste("'%s' did not return the component '%s'.\n  All %d are",
                       "required: %s."), .nl_arg_name(order), missing[1L],
                 length(want), paste(want, collapse = ", ")), call. = FALSE)
  }
  lapply(out[want], function(v) as.numeric(v) + 0 * one)
}

.nl_arg_name <- function(order) {
  c("gradient", "hessian", "deriv3", "deriv4")[order]
}

# The symbolic route at any order: D() applied once per index, which composes
# to the whole multiset. deriv() is kept for order one, where it computes the
# gradient in a single compiled pass.
.nl_symbolic_exprs <- function(expr, params, order) {
  out <- list()
  for (ix in .nl_index_list(length(params), order)) {
    e <- expr
    ok <- TRUE
    for (j in ix) {
      e <- tryCatch(stats::D(e, params[j]), error = function(err) NULL)
      if (is.null(e)) { ok <- FALSE; break }
    }
    if (!ok) return(NULL)
    out[[.nl_key(params, ix)]] <- e
  }
  out
}

# One order of derivatives of f in its own parameters, by whichever route is
# available, highest first. Returns a named list keyed by .nl_keys().
.nl_forder <- function(bp, th, order) {
  params <- bp$params
  one <- bp$one
  # order 0 is the value itself, which is what a stencil starts from when the
  # term has no analytic derivative at all -- an opaque function. One stencil
  # of order k on the value is a single difference of an analytic quantity;
  # differencing the numerical FIRST derivative would be the nesting the
  # package forbids everywhere.
  if (order == 0L) return(stats::setNames(list(.nl_eval(bp, th)$value), ""))
  sup <- bp$supplied[[.nl_arg_name(order)]]
  if (!is.null(sup)) {
    return(.nl_user_order(sup, params, order, th, bp$data_vars, one))
  }
  # order one has a route of its own already: the supplied function above, the
  # compiled gradient of deriv(), or one stencil on the value
  if (order == 1L) return(.nl_eval(bp, th)$grad)
  ex <- bp$sym_exprs[[as.character(order)]]
  if (!is.null(ex)) {
    env <- c(bp$data_vars, th)
    return(lapply(ex, function(e) as.numeric(eval(e, env)) + 0 * one))
  }
  # the numerical route: ONE stencil per distinct parameter still to be
  # differentiated in, applied to the highest ANALYTIC order there is. Two
  # stencils in DIFFERENT parameters compose and are not the nested
  # differencing the package forbids; a repeated parameter takes a single
  # stencil of its own multiplicity instead of a chain.
  m <- min(bp$analytic_order, order - 1L)
  base <- function(tv) .nl_forder(bp, tv, m)
  out <- list()
  for (ix in .nl_index_list(length(params), order)) {
    key <- .nl_key(params, ix)
    # the analytic part is the first m indices of the sorted multiset; the
    # rest is differenced
    s <- sort(ix)
    keep <- if (m > 0L) s[seq_len(m)] else integer(0)
    rest <- if (m > 0L) s[-seq_len(m)] else s
    bkey <- .nl_key(params, keep)
    out[[key]] <- .nl_stencil(base, bkey, params[rest], th, one)
  }
  out
}

# f differentiated in each parameter of `vars`, with the multiplicity of a
# repeated one taken by a single stencil of that order.
.nl_stencil <- function(base, key, vars, th, one) {
  if (!length(vars)) return(base(th)[[key]])
  p <- vars[1L]
  r <- sum(vars == p)
  rest <- vars[vars != p]
  acc <- 4L
  s <- numericals7::fd_offsets(r, acc)$central
  w <- numericals7::fd_weights(s, r)
  h <- numericals7::fd_step(th[[p]], r, acc)
  val <- 0 * one
  for (k in seq_along(s)) {
    if (w[k] == 0) next
    tv <- th
    tv[[p]] <- th[[p]] + s[k] * h
    val <- val + w[k] * .nl_stencil(base, key, rest, tv, one)
  }
  val / h^r
}


#' The Derivatives of a Nonlinear Term in Its Own Parameters
#'
#' @description
#' \eqn{\partial^k f/\partial\theta^k} at the given coefficients, one named
#' component per index multiset, keyed as \pkg{distributions7} keys its own
#' derivative surfaces.
#'
#' @details
#' The route is chosen per ORDER, highest first: a function supplied to
#' \code{\link{nl}} is used where there is one, then the symbolic route where
#' the expression can be differentiated, and otherwise one stencil applied to
#' the highest order that is analytic -- which includes a supplied one. That is
#' what makes writing out a Hessian pay twice: the third and fourth orders are
#' then one difference away from an exact second rather than from the function.
#'
#' The components are in the term's OWN parameters, not in its coefficients.
#' The chain rule onto the coefficients -- the links, and a subformula's design
#' -- belongs to the term, which is the only thing that knows them.
#'
#' @param term A built \code{\link{NlTerm}}.
#' @param coef The coefficients, or \code{NULL} for the ones the term carries.
#' @param order 1, 2, 3 or 4.
#'
#' @return A named list of numeric vectors, one per component of that order.
#'
#' @examples
#' dd <- data.frame(x = seq(0.2, 3, length.out = 20))
#' dd$y <- 2 * exp(-1.3 * dd$x)
#' b <- term_build(nl(~ a * exp(-r * x), start = list(a = 2, r = 1.3)), dd)
#' names(nl_fderiv(b, order = 2))
#'
#' @seealso \code{\link{nl}}
#'
#' @export
nl_fderiv <- function(term, coef = NULL, order = 1L) {
  if (!S7::S7_inherits(term, NlTerm)) {
    stop("'term' must be a built nl() term.", call. = FALSE)
  }
  order <- as.integer(order)
  if (length(order) != 1L || is.na(order) || order < 1L || order > 4L) {
    stop("'order' must be 1, 2, 3 or 4.", call. = FALSE)
  }
  bp <- term@blueprint
  if (!length(bp)) stop("the term is not built.", call. = FALSE)
  cf <- if (is.null(coef)) bp$coef else as.numeric(coef)
  .nl_forder(bp, .nl_theta(bp, cf)$value, order)
}


# The Jacobian in the coefficients and the value, at one coefficient vector.
# It is assembled block by block rather than written into a preallocated
# matrix: scaling a sub-design by a per-row weight keeps its sparsity
# pattern, so a sparse submodel (indicators over many groups) stays sparse
# through the block it becomes.
.nl_jacobian <- function(bp, coef) {
  th <- .nl_theta(bp, coef)
  fv <- .nl_eval(bp, th$value)
  blocks <- vector("list", length(bp$params))
  for (k in seq_along(bp$params)) {
    p <- bp$params[k]
    w <- fv$grad[[p]] * th$deriv[[p]]
    Z <- bp$Z[[p]]
    blocks[[k]] <- if (is.null(Z)) matrix(w, ncol = 1L) else w * Z
  }
  list(J = .nl_bind(blocks), value = fv$value)
}

# blocks side by side, plain when every one is plain, Matrix when any is
.nl_bind <- function(blocks) {
  if (length(blocks) == 1L) return(blocks[[1L]])
  if (!any(vapply(blocks, function(m) inherits(m, "Matrix"), logical(1)))) {
    return(do.call(cbind, blocks))
  }
  Reduce(Matrix::cbind2, blocks)
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

  # the submodel design of each parameter, and where its coefficients sit.
  # A subformula goes through the interpreter, so it takes any term of the
  # package; the built sub-terms are kept, since prediction reapplies each
  # one's own blueprint, and the penalties they carry are collected here
  # with the key parameter::subterm
  Z <- stats::setNames(vector("list", length(params)), params)
  subs <- stats::setNames(vector("list", length(params)), params)
  index <- stats::setNames(vector("list", length(params)), params)
  sub_pens <- list()
  cn <- character(0)
  pos <- 0L
  for (p in params) {
    sf <- term@subformulas[[p]]
    if (is.null(sf)) {
      pos <- pos + 1L
      index[[p]] <- pos
      cn <- c(cn, p)
    } else {
      sub <- .nl_submodel(p, sf, data)
      subs[[p]] <- sub$terms
      Z[[p]] <- sub$Z
      index[[p]] <- pos + seq_len(ncol(sub$Z))
      pos <- pos + ncol(sub$Z)
      cn <- c(cn, paste(p, sub$coef_names, sep = "."))
      for (e in sub$penalties) {
        sub_pens[[length(sub_pens) + 1L]] <- list(
          name = paste0(p, "::", e$name),
          index = index[[p]][e$index],
          penalty = e$penalty, fixed = e$fixed, n_values = e$n_values,
          values = e$values, min_ratio = e$min_ratio, search = e$search)
      }
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

  # the orders a user wrote out, checked against this term's own parameters
  # now that they are known: a name is an error at BUILD rather than a silent
  # fall-through at the first evaluation
  supplied <- term@spec$supplied
  if (is.null(supplied)) supplied <- list()
  if (length(supplied)) {
    th0 <- stats::setNames(lapply(params, function(p) numeric(n) + 0.37),
                           params)
    for (nm in names(supplied)) {
      k <- match(nm, c("gradient", "hessian", "deriv3", "deriv4"))
      .nl_user_order(supplied[[nm]], params, k, th0, data_vars, numeric(n))
    }
  }
  # the symbolic expressions, differentiated once at build and evaluated at
  # every point thereafter; orders 3 and 4 are not formed unless the term is
  # asked for them, D() on a large expression being the expensive part
  sym <- list()
  if (is_f) {
    for (k in 2:4) {
      if (!is.null(supplied[[.nl_arg_name(k)]])) next
      e <- .nl_symbolic_exprs(term@fn[[2L]], params, k)
      if (is.null(e)) break
      sym[[as.character(k)]] <- e
    }
  }
  # the highest order that is analytic, which is where a stencil starts from
  ana <- 0L
  for (k in 1:4) {
    ok <- !is.null(supplied[[.nl_arg_name(k)]]) ||
      (k == 1L && identical(mode, "symbolic")) ||
      (k >= 2L && !is.null(sym[[as.character(k)]]))
    if (ok) ana <- k else break
  }

  bp <- list(params = params, links = links, Z = Z, index = index,
             n = n, ncoef = pos, mode = mode, dexpr = dexpr,
             expr = if (is_f) term@fn[[2L]] else NULL,
             fn = if (is_f) NULL else term@fn, xval = xval,
             data_vars = data_vars, one = numeric(n),
             subformulas = term@subformulas, subs = subs,
             is_formula = is_f,
             supplied = supplied, sym_exprs = sym,
             analytic_order = max(1L, ana))

  # the penalties are the sub-terms' own, collected above with the key
  # parameter::subterm: a penalty is asked for inside the subformula,
  # where the sub-term that carries it brings its own hyperparameter
  bp$penalties <- sub_pens

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

# One parameter's submodel: the subformula through the interpreter, every
# sub-term built against the data, the blocks side by side, and the
# penalties the sub-terms declare with their indices in the bound design.
# What is kept is the BUILT terms, because prediction reapplies each one's
# own blueprint (levels, knots, constants) rather than rebuilding it.
.nl_submodel <- function(p, sf, data) {
  ir <- interpret_formula(sf, data)
  for (lb in names(ir$terms)) {
    .nl_reject_subterm(ir$terms[[lb]], p, lb)
  }
  subs <- lapply(ir$terms, term_build, data = data)
  mats <- lapply(subs, term_matrix)
  pens <- list()
  off <- 0L
  for (lb in names(subs)) {
    for (e in term_penalties(subs[[lb]])) {
      key <- if (is.null(e$name) || !nzchar(e$name)) lb
        else paste0(lb, "::", e$name)
      pens[[length(pens) + 1L]] <- list(name = key, index = off + e$index,
                                        penalty = e$penalty, fixed = e$fixed,
                                        n_values = e$n_values,
                                        values = e$values,
                                        min_ratio = e$min_ratio,
                                        search = e$search)
    }
    off <- off + ncol(mats[[lb]])
  }
  list(terms = subs, Z = .nl_bind(mats),
       coef_names = unlist(lapply(subs, term_coef_names), use.names = FALSE),
       penalties = pens)
}

# What a parameter's submodel must be: a fixed design. A structural term
# reports no block at all, and a term whose block moves with its own
# coefficients would put an iteration inside the Jacobian of another one.
.nl_reject_subterm <- function(tm, p, lb) {
  if (S7::S7_inherits(tm, structural_term)) {
    stop(sprintf(paste("the subformula of '%s' carries '%s', a structural",
                       "term; a parameter's submodel must be a fixed",
                       "design."), p, lb), call. = FALSE)
  }
  m <- S7::method(term_refresh, S7::S7_class(tm))
  cls <- attr(m, "signature")[[1L]]
  if (!identical(attr(cls, "name"), "model_term")) {
    stop(sprintf(paste("the subformula of '%s' carries '%s', whose block",
                       "moves with its coefficients; a parameter's submodel",
                       "must be a fixed design."), p, lb), call. = FALSE)
  }
  invisible(NULL)
}

#' @title Penalties of a Nonlinear Term
#' @name term_penalties.NlTerm
#' @description
#' The penalties the sub-terms of the subformulas declare, each under the
#' key \code{parameter::subterm} and covering that sub-term's coefficients
#' in the term's own numbering. The list is empty when there are none, and
#' for a specification, whose parameters a formula does not name until the
#' data say which of them the data supply.
#' @param term A built \code{\link{NlTerm}}.
#' @param ... Unused.
#' @return A list of entries, as \code{\link{term_penalties}} documents.
#' @keywords internal
S7::method(term_penalties, NlTerm) <- function(term, ...) {
  pens <- term@blueprint$penalties
  if (is.null(pens)) list() else pens
}

# the blueprint of the same term read on other rows: the submodel designs
# are REAPPLIED through each sub-term's own blueprint, never rebuilt
.nl_blueprint_at <- function(term, newdata) {
  bp <- term@blueprint
  nb <- bp
  nb$n <- nrow(newdata)
  nb$one <- numeric(nrow(newdata))
  nb$data_vars <- as.list(newdata)
  for (p in bp$params) {
    if (!is.null(bp$Z[[p]])) {
      nb$Z[[p]] <- .nl_bind(lapply(bp$subs[[p]], term_predict,
                                   newdata = newdata))
    }
  }
  if (!bp$is_formula) {
    nb$xval <- if (!is.null(term@spec$x)) {
      eval(term@spec$x, newdata, baseenv())
    } else newdata
  }
  nb
}

S7::method(term_predict, NlTerm) <- function(term, newdata, ...) {
  .assert_built(term)
  X <- .nl_jacobian(.nl_blueprint_at(term, newdata), term@blueprint$coef)$J
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
#' @details
#' A term whose contribution \eqn{f(x; \theta)} is nonlinear in its
#' parameters is carried by the linearization
#'
#' \deqn{f(x; \theta + h) \approx f(x; \theta) + J(\theta) h,
#'   \qquad J(\theta)_{ij}
#'     = \frac{\partial f(x_i; \theta)}{\partial \theta_j},}
#'
#' so the design block at the current \eqn{\theta} is \eqn{J(\theta)} and
#' the coefficient it multiplies is the increment \eqn{h}. Refreshing at
#' the new \eqn{\theta} and solving the linear problem again is the
#' Gauss-Newton iteration; \code{\link{term_value}} reports
#' \eqn{f(x; \theta)} itself, which is the other half a step needs. For a
#' break-point term the same shape holds with a different block: the
#' Jacobian for \code{\link{seg}}, and the frozen-weight columns of
#' \code{\link{jump}}, from which the break-point is read rather than
#' incremented.
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


#' How a Term's Block Moves With Its Coefficients
#'
#' @description
#' \eqn{\sum_{i,j} A_{ij}\,\partial X_{ij}/\partial\beta_c}, one entry per
#' coefficient of the term: the block's own derivative contracted against a
#' weight the caller supplies.
#'
#' @details
#' A term whose block is a fixed design does not move, and the base method
#' returns zeros. A term registering \code{\link{term_refresh}} does move, and
#' a consumer that differentiates anything built from \eqn{X} needs this: a
#' marginal criterion reads \eqn{\log|K|} with
#' \eqn{K = -\sum_i w_i \ell_{ab,i} X_a X_b'}, so
#' \eqn{\partial K/\partial\beta} gains everything coming from
#' \eqn{\partial X/\partial\beta}.
#'
#' \strong{The contraction and not the derivative}, because
#' \eqn{\partial X/\partial\beta} is \eqn{n \times m \times m} and nothing
#' needs it whole. For \code{\link{nl}} it is closed form at \eqn{O(nm)}:
#' writing \eqn{X_{i,c} = w_{p}(i) Z_p[i,c]} with
#' \eqn{w_p = (\partial f/\partial\theta_p)\,h_p'},
#' \deqn{\frac{\partial X_{i,c_1}}{\partial\beta_{c_2}} =
#'   Z_{p_1}[i,c_1]Z_{p_2}[i,c_2]\Big(f_{p_1p_2}h_{p_1}'h_{p_2}'
#'   + \delta_{p_1p_2} f_{p_1} h_{p_1}''\Big),}
#' so with \eqn{s_{p}(i) = \sum_{c\in p} A_{ic}Z_{p}[i,c]} the answer is
#' \eqn{Z_{p_2}'\sum_{p_1} q_{p_1p_2}s_{p_1}}, one crossprod per parameter.
#' The chain rule onto the coefficients -- each parameter's link and its
#' subformula's design -- stays here, which is the only place that knows them.
#'
#' ⚠️ The base method is right for a FIXED block and is what
#' \code{\link{seg}}, \code{\link{jump}} and \code{\link{jseg}} inherit, whose
#' blocks do move. Theirs is not a difference a caller could take either: a
#' break-point column is a step function in its break-point, so the quotient
#' diverges as the step shrinks -- measured at h, h/4 and h/16, 3.6e4, 1.4e5
#' and 5.8e5, against \code{\link{nl}}'s 0.6038 throughout. They need their own
#' method, written from the closed forms, which is not this release.
#'
#' @param term A built term.
#' @param coef The coefficients, or \code{NULL} for the ones it carries.
#' @param A A numeric matrix, one row per observation and one column per
#'   coefficient of the term.
#' @param ... Passed to methods.
#'
#' @return A numeric vector as long as the term's coefficients.
#'
#' @examples
#' dd <- data.frame(x = seq(0.2, 3, length.out = 20))
#' dd$y <- 2 * exp(-1.3 * dd$x)
#' b <- term_build(nl(~ a * exp(-r * x), start = list(a = 2, r = 1.3)), dd)
#' term_block_contract(b, A = matrix(1, 20, 2))
#'
#' @seealso \code{\link{term_refresh}}, \code{\link{nl_fderiv}}
#'
#' @export
term_block_contract <- S7::new_generic(
  "term_block_contract", "term",
  function(term, coef = NULL, A, ...) S7::S7_dispatch())

S7::method(term_block_contract, model_term) <- function(term, coef = NULL, A,
                                                        ...) {
  numeric(length(term_coef_names(term)))
}

S7::method(term_block_contract, NlTerm) <- function(term, coef = NULL, A, ...) {
  bp <- term@blueprint
  if (!length(bp)) stop("the term is not built.", call. = FALSE)
  cf <- if (is.null(coef)) bp$coef else as.numeric(coef)
  A <- as.matrix(A)
  if (nrow(A) != bp$n || ncol(A) != bp$ncoef) {
    stop(sprintf("'A' must be %d by %d, the observations by the term's coefficients.",
                 bp$n, bp$ncoef), call. = FALSE)
  }
  params <- bp$params
  th <- .nl_theta(bp, cf)
  f1 <- .nl_forder(bp, th$value, 1L)
  f2 <- .nl_forder(bp, th$value, 2L)
  # s_p, one column per parameter: the weight already summed over that
  # parameter's own sub-design
  s <- lapply(params, function(p) {
    idx <- bp$index[[p]]
    Z <- bp$Z[[p]]
    Ap <- A[, idx, drop = FALSE]
    if (is.null(Z)) as.numeric(Ap[, 1L]) else
      as.numeric(rowSums(Ap * as.matrix(Z)))
  })
  names(s) <- params
  out <- numeric(bp$ncoef)
  for (p2 in params) {
    acc <- 0 * bp$one
    for (p1 in params) {
      q <- as.numeric(f2[[.nl_key(params, c(match(p1, params),
                                            match(p2, params)))]]) *
        th$deriv[[p1]] * th$deriv[[p2]]
      if (identical(p1, p2)) {
        q <- q + as.numeric(f1[[p1]]) * th$deriv2[[p1]]
      }
      acc <- acc + q * s[[p1]]
    }
    idx <- bp$index[[p2]]
    Z <- bp$Z[[p2]]
    out[idx] <- if (is.null(Z)) sum(acc) else
      as.numeric(crossprod(as.matrix(Z), acc))
  }
  out
}

#' @title Has a Term's Own Iteration Settled?
#'
#' @description
#' \code{TRUE} when a term whose block is refreshed has nothing further to
#' say about where its own parameters are. A term whose block does not move
#' answers \code{TRUE}, having no iteration of its own.
#'
#' @details
#' It exists because a score cannot always answer the question. Where the
#' block is the Jacobian of the contribution -- \code{\link{nl}},
#' \code{\link{seg}} -- the gradient of the model's objective is the block
#' times the derivative of the log-likelihood in the predictor, and its
#' vanishing is the test. Where the block is a working LINEARIZATION with a
#' frozen weight -- \code{\link{jump}}, \code{\link{jseg}} -- it is not: the
#' profile objective of a discontinuous term is a step function in the
#' break-point, so it has no gradient to vanish, and the quantity the
#' iteration actually drives to zero is the movement of the break-point.
#' That movement is what \code{\link{seg_converged}} reads, and a fitting
#' layer asks for it here without knowing which construction it holds.
#'
#' @param term A built term.
#' @param ... Passed to methods.
#'
#' @return A single logical.
#'
#' @examples
#' dd <- data.frame(x = seq(0, 2, length.out = 20))
#' term_converged(term_build(linpar(~x), dd))
#'
#' @seealso \code{\link{term_refresh}}, \code{\link{seg_converged}}
#' @export
term_converged <- S7::new_generic("term_converged", "term",
  function(term, ...) S7::S7_dispatch())

S7::method(term_converged, model_term) <- function(term, ...) TRUE

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
#' @details
#' \code{newdata} asks for the same contribution on other rows, and is what
#' a predictor needs where the block is a Jacobian: there
#' \code{term_predict()} times the coefficients is the linearization and
#' not the contribution, and the two differ by whatever the linearization
#' drops. For \code{\link{seg}} that difference is a step at the
#' break-point in a construction that is continuous. Rows arriving here are
#' treated as \code{\link{term_predict}} treats them, through the levels
#' and constants the blueprint recorded, never rebuilt.
#'
#' @param term A built term.
#' @param coef The coefficients. Optional for a nonlinear term, which
#'   carries the ones it was last refreshed at.
#' @param newdata An optional data frame; the contribution is returned on
#'   its rows instead of on the ones the term was built from.
#' @param ... Passed to methods.
#'
#' @return A numeric vector, one value per observation.
#'
#' @seealso \code{\link{term_refresh}}, \code{\link{term_predict}}
#'
#' @examples
#' dd <- data.frame(x = seq(0, 2, length.out = 20))
#' built <- term_build(nl(~ a * exp(-r * x), start = list(a = 2, r = 1)), dd)
#' head(term_value(built), 3)
#' head(term_value(built, newdata = data.frame(x = c(0, 1, 2))), 3)
#'
#' @export
term_value <- S7::new_generic("term_value", "term",
  function(term, coef = NULL, newdata = NULL, ...) S7::S7_dispatch())

S7::method(term_value, additive_term) <- function(term, coef = NULL,
                                                  newdata = NULL, ...) {
  .assert_built(term)
  if (is.null(coef)) {
    stop("'coef' is required: a linear term contributes the block times it.",
         call. = FALSE)
  }
  X <- if (is.null(newdata)) term@X else term_predict(term, newdata)
  as.numeric(X %*% coef)
}

S7::method(term_value, NlTerm) <- function(term, coef = NULL, newdata = NULL,
                                           ...) {
  .assert_built(term)
  bp <- if (is.null(newdata)) term@blueprint else
    .nl_blueprint_at(term, newdata)
  if (is.null(coef)) {
    if (is.null(newdata)) return(bp$value)
    coef <- term@blueprint$coef
  }
  .nl_jacobian(bp, as.numeric(coef))$value
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

#' @title Where a Nonlinear Term's Coefficients Begin
#' @name term_coef_start.NlTerm
#' @description
#' The coefficients \code{\link{term_build}} built the block at: the
#' starting value of each of the term's own parameters carried through its
#' link, from \code{start} where the caller gave one. The block is the
#' Jacobian at those values, so it is not the same block at any other
#' point, and a fitting layer that started at zero would linearize where
#' the term was never meant to be evaluated.
#' @param term A built \code{\link{NlTerm}}.
#' @param ... Unused.
#' @return A numeric vector, one value per column of the block.
#' @keywords internal
S7::method(term_coef_start, NlTerm) <- function(term, ...) {
  .assert_built(term)
  term@blueprint$coef
}
