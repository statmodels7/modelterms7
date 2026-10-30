#' @include gas.R
NULL

#' A Joint Score-Driven Filter
#'
#' @description
#' The S7 class of several [gas()] terms that sit in the equations of
#' different distribution parameters and are filtered as one recursion. Build
#' one with [gas_joint()].
#'
#' @details
#' A score-driven level in the equation of one parameter is driven by the
#' score of that parameter, which is evaluated at the current values of all
#' the others. Where two parameters both carry a level, each score depends on
#' the other level, so the two recursions are one: with \eqn{G} the filtered
#' parameters,
#' \deqn{f_{k,t} = \omega_k + \sum_{i=1}^{p_k} \kappa_{k,i}\, u_{k,t-i}
#'   + \sum_{l \in C_k} \sum_{i=1}^{p_{kl}} \kappa_{k,i,l}\, u_{l,t-i}
#'   + \sum_{j=1}^{q_k} \xi_{k,j} f_{k,t-j}, \qquad k \in G,}
#' where \eqn{u_t} is the vector of driving scores of the parameters in
#' \eqn{G}, read at the predictor of every equation at time \eqn{t},
#' \eqn{C_k} is the set of the parameters to which the `p_cross` of the term
#' of \eqn{k} gives a positive count, and \eqn{p_{kl}} is that count. The
#' persistence of each level is its own, so the autoregressive part is
#' stationary when each equation's is.
#'
#' The driving vector is the caller's: a fitting layer supplies
#' \eqn{u_t = \mathcal{I}_{GG,t}^{-d}\,\partial\ell_t/\partial\eta_{G,t}} with
#' the expected information of the filtered parameters and its derivative in
#' the predictors.
#'
#' @inheritParams model_term
#' @param terms A named list of built [GasTerm()] objects, one per filtered
#'   parameter, named by the parameter.
#' @param params The names of the filtered parameters, in the order of
#'   `terms`.
#' @param offsets The position before the first parameter of each term in the
#'   joint parameter vector.
#' @param partners A list with one integer vector per term: the positions in
#'   `params` of the parameters whose scores the term's level loads besides
#'   its own.
#' @param keys The names under which a fitting layer knew each term before
#'   joining them, named by the parameter, or an empty vector.
#' @param build_info The shared ordering: `order` (the rows of each group, in
#'   time order) and `n`. The property is inherited from [structural_term()].
#'
#' @return An S7 object of class `GasJoint`, inheriting from
#'   [structural_term()], so that a fitting layer treats the joint filter as
#'   one structural term whose parameters are those of its terms, each
#'   prefixed by the distribution parameter of its term.
#'
#' @seealso [gas_joint()], the constructor; [gas()] for one term.
#'
#' @examples
#' set.seed(1)
#' dd <- data.frame(t = 1:40, y = rnorm(40))
#' tm <- list(mu = term_build(gas(p = 1, q = 1, time = t), dd),
#'            sigma = term_build(gas(p = 1, q = 1, time = t), dd))
#' jt <- gas_joint(tm)
#' S7::S7_inherits(jt, GasJoint)
#' jt@params
#'
#' @export
GasJoint <- S7::new_class(
  name = "GasJoint",
  parent = structural_term,
  properties = list(
    terms = S7::class_list,
    params = S7::class_character,
    offsets = S7::class_integer,
    partners = S7::class_list,
    keys = S7::class_character
  )
)

#' Filter Several Score-Driven Terms as One Recursion
#'
#' @description
#' Checks that a set of built [gas()] terms can be filtered together and
#' returns them as a [GasJoint()] object, with the cross loadings that each
#' term's `p_cross` declares added to its parameters.
#'
#' @details
#' The terms must walk the same series: their `by` and `time` must give the
#' same groups and the same order, because each score is read at the
#' predictors of every equation at the same time. A term may develop its
#' level, its loadings, its persistence or its cross loadings with a
#' subformula.
#'
#' The terms must share one `scaling`, because the driving vector is scaled
#' by a power of the expected information of all the filtered parameters
#' together. With more than one term the power is 0, 1/2 or 1: for these the
#' derivatives of the matrix power are exact and short (the identity, the
#' symmetric square root through a Sylvester equation per order, and the
#' inverse).
#'
#' The `p_cross` of a term names the parameters whose scores its level loads
#' besides its own, with the number of lags of each. A parameter that it
#' names with a positive count must be one of the joined terms or one of
#' `drivers`, and a term cannot name its own parameter, whose score it loads
#' through its own loadings. The cross loadings of a term follow the
#' parameters in the order of the joined terms and then of `drivers`, lag by
#' lag.
#'
#' A parameter in `drivers` has no level of its own: its predictor is static,
#' and its score enters the driving vector so that the levels of the terms
#' can load it. It is carried as a [GasDriver()] member after the terms, with
#' no parameters and a level of zero, which is a score-driven term whose
#' loadings and persistence are zero. With a `scaling` other than zero the
#' driving vector is scaled by the expected information of the filtered
#' parameters and the drivers together.
#'
#' @param terms A named list of built [GasTerm()] objects, named by the
#'   distribution parameter in whose equation each one sits.
#' @param drivers The names of the distribution parameters without a
#'   [gas()] term whose scores the terms load, or an empty vector. A fitting
#'   layer, which knows the family, checks that they are parameters of it.
#'
#' @return A [GasJoint()] object. Its `terms` are the input terms with
#'   `cross_with` filled, so their [term_params()] include the cross
#'   loadings, followed by one [GasDriver()] per entry of `drivers`; its
#'   `params` are the names of both.
#'
#' @seealso [GasJoint()], [gas()]
#'
#' @examples
#' set.seed(1)
#' dd <- data.frame(t = 1:40, y = rnorm(40))
#' tm <- list(mu = term_build(gas(p = 1, q = 1, time = t,
#'                                p_cross = c(sigma = 2)), dd),
#'            sigma = term_build(gas(p = 1, q = 1, time = t), dd),
#'            nu = term_build(gas(p = 1, q = 1, time = t), dd))
#' jt <- gas_joint(tm)
#' lapply(jt@terms, term_params)
#'
#' # The level of mu loads the score of sigma, whose predictor is static.
#' jd <- gas_joint(list(mu = term_build(gas(p = 1, q = 1, time = t,
#'                                          p_cross = c(sigma = 1)), dd)),
#'                 drivers = "sigma")
#' jd@params
#' term_params(jd)
#'
#' @export
gas_joint <- function(terms, drivers = character(0)) {
  if (!is.list(terms) || !length(terms) || is.null(names(terms)) ||
      anyNA(names(terms)) || !all(nzchar(names(terms))) ||
      anyDuplicated(names(terms))) {
    stop("'terms' must be a list of gas() terms named by distinct parameters.",
         call. = FALSE)
  }
  params <- names(terms)
  for (k in params) {
    tm <- terms[[k]]
    if (!S7::S7_inherits(tm, GasTerm)) {
      stop(sprintf("the term of '%s' is not a gas() term.", k), call. = FALSE)
    }
    if (!length(tm@build_info)) {
      stop(sprintf("the gas() term of '%s' has not been built.", k),
           call. = FALSE)
    }
  }
  if ((length(terms) > 1L || length(drivers)) &&
      any(vapply(terms, .gas_static, logical(1))) &&
      any(vapply(terms, function(tm) length(tm@build_info$sub) > 0L,
                 logical(1)))) {
    stop(paste("init = \"static\" is not yet available for a joint filter",
               "whose terms carry subformulas."), call. = FALSE)
  }
  bp1 <- terms[[1L]]@build_info
  for (k in params[-1L]) {
    bp <- terms[[k]]@build_info
    if (!identical(bp$n, bp1$n) || !identical(bp$order, bp1$order)) {
      stop(sprintf(paste0(
        "the gas() terms of '%s' and '%s' walk different series:\n",
        "  their 'by' and 'time' must give the same groups in the same ",
        "order."), params[1L], k), call. = FALSE)
    }
  }
  d <- vapply(terms, function(tm) tm@scaling, numeric(1))
  if (length(unique(d)) > 1L) {
    stop(sprintf(paste0(
      "the gas() terms have different 'scaling' (%s):\n",
      "  the scores of the filtered parameters are scaled together, so ",
      "every term must have the same 'scaling'."),
      paste(sprintf("%s = %s", params, format(d)), collapse = ", ")),
      call. = FALSE)
  }
  if (!is.character(drivers) || anyNA(drivers) || !all(nzchar(drivers)) ||
      anyDuplicated(drivers)) {
    stop("'drivers' must be a vector of distinct parameter names.",
         call. = FALSE)
  }
  both <- intersect(drivers, params)
  if (length(both)) {
    stop(sprintf(paste0(
      "'%s' carries a gas() term and is named in 'drivers'; a driver has ",
      "no level of its own."), both[1L]), call. = FALSE)
  }
  if (length(terms) + length(drivers) > 1L && !(d[[1L]] %in% c(0, 0.5, 1))) {
    stop(sprintf(paste0(
      "'scaling = %s' is available for a single gas() term only; with ",
      "several\n  filtered parameters 'scaling' must be 0, 1/2 or 1."),
      format(d[[1L]])), call. = FALSE)
  }
  nt <- length(params)
  all_par <- c(params, drivers)
  partners <- vector("list", length(all_par))
  for (k in seq_len(nt)) {
    tm <- terms[[k]]
    pc <- tm@p_cross
    if (params[k] %in% names(pc)) {
      stop(sprintf(paste0(
        "the gas() term of '%s' names its own parameter in 'p_cross';\n",
        "  the level loads its own score through %s."), params[k],
        if (tm@p == 1L) "kappa1" else sprintf("kappa1, ..., kappa%d", tm@p)),
        call. = FALSE)
    }
    on <- names(pc)[pc > 0L]
    miss <- setdiff(on, all_par)
    if (length(miss)) {
      stop(sprintf(paste0(
        "the gas() term of '%s' loads the score of '%s' (p_cross), but '%s' ",
        "has no gas() term\n  and is not among 'drivers'."),
        params[k], miss[1L], miss[1L]), call. = FALSE)
    }
    # the partners in the order of the joined terms and then of the drivers,
    # which orders the cross loadings of the term
    partners[[k]] <- which(all_par %in% on)
    tm@cross_with <- all_par[partners[[k]]]
    terms[[k]] <- tm
  }
  for (k in seq_along(drivers)) {
    partners[[nt + k]] <- integer(0)
    terms[[drivers[k]]] <- .gas_driver(terms[[1L]], d[[1L]])
  }
  np <- vapply(terms, function(tm) length(term_params(tm)), integer(1))
  GasJoint(label = "gas_joint", terms = terms, params = all_par,
           offsets = as.integer(cumsum(np) - np), partners = partners,
           build_info = list(order = bp1$order, n = bp1$n))
}

# a driver beside the term `like`: its series, its scaling, no parameters
.gas_driver <- function(like, scaling) {
  bp <- like@build_info
  bp$sub <- NULL
  GasDriver(label = "driver", p = 0L, q = 0L, by = like@by, time = like@time,
            links = list(), submodels = list(), scaling = as.double(scaling),
            p_cross = stats::setNames(integer(0), character(0)),
            cross_with = character(0), build_info = bp)
}

# whether a member of a joint filter is a driver
.gas_is_driver <- function(tm) S7::S7_inherits(tm, GasDriver)

# the members of a joint filter that are terms, not drivers
.gas_joint_terms_only <- function(term) {
  Filter(function(tm) !.gas_is_driver(tm), term@terms)
}

#' A Static Parameter Whose Score Drives a Joint Filter
#'
#' @description
#' The S7 class of a member of a [GasJoint()] that carries no level: the
#' predictor of its parameter is static, and its score enters the driving
#' vector so that the levels of the other members can load it through their
#' `p_cross`. [gas_joint()] creates it from its argument `drivers`.
#'
#' @details
#' A driver is a score-driven term with no parameters, whose level is zero
#' at every observation: [term_params()], [term_links()], [term_start()],
#' [term_penalties()] and [term_components()] are empty, and the recursion
#' gives it no level and no derivative.
#'
#' @inheritParams GasTerm
#'
#' @return An S7 object of class `GasDriver`, inheriting from [GasTerm()].
#'
#' @seealso [gas_joint()]
#'
#' @examples
#' set.seed(1)
#' dd <- data.frame(t = 1:40, y = rnorm(40))
#' jd <- gas_joint(list(mu = term_build(gas(p = 1, q = 1, time = t,
#'                                          p_cross = c(sigma = 1)), dd)),
#'                 drivers = "sigma")
#' S7::S7_inherits(jd@terms$sigma, GasDriver)
#' term_params(jd@terms$sigma)
#'
#' @export
GasDriver <- S7::new_class(name = "GasDriver", parent = GasTerm)

#' @title The Empty Parameters of a Driver
#' @name term_params.GasDriver
#' @description
#' A [GasDriver()] has no parameters: its parameters, links, start,
#' penalties, components, level parameter and readable quantities are empty.
#' @param term A [GasDriver()].
#' @param zeta Unused.
#' @param ... Unused.
#' @return An empty vector or list of the type the generic returns.
#' @keywords internal
S7::method(term_params, GasDriver) <- function(term, ...) character(0)
S7::method(term_links, GasDriver) <- function(term, ...) list()
S7::method(term_start, GasDriver) <- function(term, ...) {
  stats::setNames(numeric(0), character(0))
}
S7::method(term_penalties, GasDriver) <- function(term, ...) list()
S7::method(term_components, GasDriver) <- function(term, ...) list()
S7::method(term_level_param, GasDriver) <- function(term, ...) character(0)
S7::method(term_readable, GasDriver) <- function(term, zeta, ...) {
  list(name = character(0), value = numeric(0),
       jacobian = matrix(0, 0, 0), scale = list())
}

#' @title The Parameters of a Joint Score-Driven Filter
#' @name term_params.GasJoint
#' @description
#' Returns the parameters of every term of the filter, in the order of the
#' terms, each prefixed by the distribution parameter of its term, as in
#' `mu.omega`.
#' @param term A [GasJoint()].
#' @param ... Unused.
#' @return A character vector.
#' @keywords internal
S7::method(term_params, GasJoint) <- function(term, ...) {
  .gas_joint_prefixed(term, function(tm) term_params(tm))
}

# a per-member character vector prefixed by the member's parameter, over the
# members that have parameters (a driver has none, and paste0() would turn
# its empty vector into "sigma.")
.gas_joint_prefixed <- function(term, f) {
  out <- character(0)
  for (k in term@params) {
    tm <- term@terms[[k]]
    if (.gas_is_driver(tm)) next
    v <- f(tm)
    if (length(v)) out <- c(out, paste0(k, ".", v))
  }
  out
}

# whether any term of the filter develops a parameter with a subformula,
# which makes the recursion's coefficients vary by observation
.gas_joint_has_sub <- function(term) {
  any(vapply(term@terms, function(tm)
    length(tm@build_info) > 0L && !is.null(tm@build_info$sub), logical(1)))
}

# The parameters as a list with one entry per filtered parameter. A caller
# may hand them as one flat list named as term_params() of the joint filter
# names them (mu.omega, ...), which is how a fitting layer keeps them.
.gas_joint_psi <- function(term, psi) {
  # a driver has no parameters, so a caller need not name it
  drv <- .gas_joint_drivers(term)
  if (all(setdiff(term@params, drv) %in% names(psi))) {
    for (k in drv) psi[[k]] <- list()
    return(psi)
  }
  flat <- unlist(psi)
  out <- list()
  for (k in term@params) {
    own <- term_params(term@terms[[k]])
    out[[k]] <- if (length(own))
      as.list(stats::setNames(flat[paste0(k, ".", own)], own)) else list()
  }
  out
}

# The coefficients of the joint recursion from the terms' parameters: per
# term, the chart of .gas_coefs() and the cross loadings, with the columns
# of the joint vector at which each coefficient's derivative sits. The cross
# loadings are a list of pairs in the order of .gas_cross_params(): pair h
# loads the score of the parameter at position xpar[h] of term@params at lag
# xlag[h], with the value cross[h] and its column col_cross[h].
.gas_joint_coefs <- function(term, psi) {
  psi <- .gas_joint_psi(term, psi)
  out <- vector("list", length(term@params))
  for (k in seq_along(term@params)) {
    tm <- term@terms[[term@params[k]]]
    if (.gas_is_driver(tm)) {
      out[[k]] <- .gas_driver_coefs()
      next
    }
    nm <- term_params(tm)
    v <- psi[[term@params[k]]]
    v <- unlist(v[nm])
    if (length(v) != length(nm) || anyNA(v)) {
      stop(sprintf("'psi$%s' must supply %s.", term@params[k],
                   paste(nm, collapse = ", ")), call. = FALSE)
    }
    p <- tm@p
    q <- tm@q
    base <- .gas_base_params(p, q)
    cf <- .gas_coefs(v[base], p, q)
    off <- term@offsets[k]
    nb <- length(base)
    pk <- term@partners[[k]]
    xp <- .gas_cross_pairs(tm)
    nx <- length(xp$name)
    out[[k]] <- list(p = p, q = q, omega = cf$omega, a = cf$a, b = cf$b,
                     db = cf$db, f0 = cf$f0, df0 = cf$df0,
                     own = off + seq_len(nb), col_omega = off + 1L,
                     col_a = off + cf$i_a, partners = as.integer(pk),
                     xpar = match(xp$partner, term@params),
                     xlag = xp$lag,
                     cross = as.numeric(v[nb + seq_len(nx)]),
                     col_cross = as.integer(off + nb + seq_len(nx)))
  }
  out
}

# the start of each member of a joint filter: `direct` and `d0`, the first
# level of each group under init = "static" (c minus the static predictor of
# the member's equation at the group's first row), from the static
# predictors `eta` (n x K)
.gas_joint_start <- function(term, cf, eta) {
  K <- length(term@params)
  for (k in seq_len(K)) {
    tm <- term@terms[[term@params[k]]]
    st <- !.gas_is_driver(tm) && .gas_static(tm)
    cf[[k]]$direct <- st
    cf[[k]]$d0 <- numeric(0)
    if (st) {
      if (is.null(eta)) {
        stop(paste("the static start of a joint filter is read from the",
                   "static predictors, and 'eta' was not given."),
             call. = FALSE)
      }
      e <- matrix(as.numeric(eta), term@build_info$n, K)
      cf[[k]]$d0 <- .gas_static_start(tm, e[, k])
    }
  }
  cf
}

# the coefficients of a driver: a level of zero with no lags and no columns
# (col_omega = 0 says that the constant has no derivative)
.gas_driver_coefs <- function() {
  list(p = 0L, q = 0L, omega = 0, a = numeric(0), b = numeric(0),
       db = matrix(0, 0, 0), f0 = 0, df0 = numeric(0), own = integer(0),
       col_omega = 0L, col_a = integer(0), partners = integer(0),
       xpar = integer(0), xlag = integer(0), cross = numeric(0),
       col_cross = integer(0))
}

#' @title Filter Several Score-Driven Terms Jointly
#' @name term_filter.GasJoint
#' @description
#' Runs the joint recursion over each group in time order and returns the
#' filtered predictors together with their exact derivatives in the joint
#' parameter vector, propagated alongside the state.
#' @details
#' With \eqn{F_{k,t}} the derivative of the level of parameter \eqn{k} and
#' \eqn{C_t = \partial u_t/\partial\eta_{G,t}^\top} the derivative of the
#' driving vector, the derivative of the driving score of parameter \eqn{l}
#' is \eqn{\sum_k C_{t,lk} F_{k,t}}, because the static predictors do not
#' depend on the terms' parameters. The recursion of \eqn{F} is the
#' recursion of the levels differentiated term by term.
#' @param term A [GasJoint()].
#' @param eta The static part of the predictors of the filtered parameters, a
#'   matrix with one column per entry of `term@params`.
#' @param y The response, unused directly.
#' @param score A function `function(e, i)` of the vector `e` of the filtered
#'   predictors at row `i`, returning the driving vector \eqn{u_i}, one value
#'   per filtered parameter.
#' @param curvature A function `function(e, i)` returning the matrix
#'   \eqn{\partial u_i/\partial e^\top}, whose entry \eqn{(l, k)} is the
#'   derivative of the driving score of parameter \eqn{l} in the predictor of
#'   parameter \eqn{k}.
#' @param psi A named list with one entry per filtered parameter, each the
#'   term's parameters on the parameter scale, named as [term_params()] of
#'   the term; or one flat list named as [term_params()] of the joint filter.
#' @param ... Unused.
#' @param fast `NULL`, or the context of the compiled route, a list with
#'   `family` (the name under which the C registry of \pkg{distributions7}
#'   knows the distribution), `constants` (its constants, each of length 1
#'   or one per observation), `links` (the class name of each parameter's
#'   link), `link_par` (a list of each link's own parameters), `bounds` (a
#'   list of each link's bounds), `y`, `eta` (the static predictors of every
#'   distribution parameter, one column each), `filtered` (the positions of
#'   the filtered parameters among them) and `scaling`. Where the joint
#'   registry of \pkg{distributions7} and the registry of
#'   \pkg{linkfunctions7} cover the family and the links, the recursion runs
#'   compiled and reads the family through them instead of the callbacks;
#'   where they do not, the context is inert and the callbacks run in R. A
#'   scaled filter also needs the family's expected information and its
#'   derivative in the registry, which some families do not have (their
#'   expected information is a quadrature); their scaled filters run in R.
#' @param threads How many threads the compiled recursion may use, over
#'   groups, when the distribution's entries are safe on a worker thread and
#'   there are at least eight groups. Each group writes only to its own rows,
#'   so the result does not depend on the count.
#' @return A list with `eta` (a matrix like the input), `jacobian` (a named
#'   list with one matrix per filtered parameter, one row per observation and
#'   one column per entry of [term_params()] of the joint filter) and `curv`
#'   (an array of dimension `c(n, K, K)` holding the curvature read at each
#'   row).
#' @keywords internal
S7::method(term_filter, GasJoint) <- function(term, eta, y, score, curvature,
                                              psi, ..., fast = NULL,
                                              threads = 1L) {
  # coefficients that vary by observation take the route of
  # gas_joint_sub.R, with the same result
  if (.gas_joint_has_sub(term)) {
    return(.gas_joint_sub_filter(term, eta, y, score, curvature, psi,
                                 fast = fast, threads = threads))
  }
  bp <- term@build_info
  K <- length(term@params)
  n <- bp$n
  eta <- matrix(as.numeric(eta), n, K)
  cf <- .gas_joint_start(term, .gas_joint_coefs(term, psi), eta)
  M <- length(term_params(term))
  nm <- term_params(term)
  if (!is.null(fast)) {
    out <- gas_joint_filter_cpp(eta, cf, M, bp$order, fast,
                                as.integer(threads))
    if (!is.null(out)) {
      jac <- lapply(out$jacobian, function(J) { colnames(J) <- nm; J })
      names(jac) <- term@params
      eta_out <- out$eta
      colnames(eta_out) <- term@params
      return(list(eta = eta_out, jacobian = jac, curv = out$curv))
    }
  }
  eta_out <- eta
  jac <- lapply(seq_len(K), function(k) matrix(0, n, M))
  cv <- array(0, c(n, K, K))
  for (gi in seq_along(bp$order)) {
    rows <- bp$order[[gi]]
    T_ <- length(rows)
    f <- matrix(0, K, T_)
    u <- matrix(0, K, T_)
    Fd <- array(0, c(K, T_, M))
    Sd <- array(0, c(K, T_, M))
    for (t in seq_len(T_)) {
      for (k in seq_len(K)) {
        c_k <- cf[[k]]
        # a static start: the first level is given and does not depend on
        # the parameters, and the levels before it equal it
        if (isTRUE(c_k$direct) && t == 1L) {
          f[k, t] <- c_k$d0[[gi]]
          Fd[k, t, ] <- 0
          next
        }
        ft <- c_k$omega
        Ft <- numeric(M)
        if (c_k$col_omega > 0L) Ft[c_k$col_omega] <- 1
        for (i in seq_len(c_k$p)) {
          if (t - i < 1L) next
          ft <- ft + c_k$a[[i]] * u[k, t - i]
          Ft <- Ft + c_k$a[[i]] * Sd[k, t - i, ]
          Ft[c_k$col_a[i]] <- Ft[c_k$col_a[i]] + u[k, t - i]
        }
        for (h in seq_along(c_k$xlag)) {
          i <- c_k$xlag[h]
          if (t - i < 1L) next
          l <- c_k$xpar[h]
          ft <- ft + c_k$cross[h] * u[l, t - i]
          Ft <- Ft + c_k$cross[h] * Sd[l, t - i, ]
          cc <- c_k$col_cross[h]
          Ft[cc] <- Ft[cc] + u[l, t - i]
        }
        for (j in seq_len(c_k$q)) {
          if (t - j >= 1L) {
            f_l <- f[k, t - j]
            F_l <- Fd[k, t - j, ]
          } else if (isTRUE(c_k$direct)) {
            f_l <- c_k$d0[[gi]]
            F_l <- numeric(M)
          } else {
            f_l <- c_k$f0
            F_l <- numeric(M)
            F_l[c_k$own] <- c_k$df0
          }
          ft <- ft + c_k$b[[j]] * f_l
          Ft <- Ft + c_k$b[[j]] * F_l
          Ft[c_k$own] <- Ft[c_k$own] + c_k$db[j, ] * f_l
        }
        f[k, t] <- ft
        Fd[k, t, ] <- Ft
      }
      row <- rows[t]
      e_t <- eta[row, ] + f[, t]
      u[, t] <- score(e_t, row)
      C <- matrix(curvature(e_t, row), K, K)
      cv[row, , ] <- C
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

#' @title Filter Several Score-Driven Terms Backwards
#' @name term_adjoint.GasJoint
#' @description
#' Runs the joint recursion of [term_filter()] in reverse and returns the
#' derivative of a caller's objective with respect to the static predictors
#' of the filtered parameters and with respect to the driving scores.
#' @details
#' With \eqn{\bar f_{k,t}} the derivative of the objective in the level and
#' \eqn{\bar u_{l,t}} in the driving score, the objective reaches the filtered
#' predictor \eqn{e_{k,t}} directly through `g` and through every driving
#' score read there, \eqn{\bar e_{k,t} = g_{k,t} + \sum_l \bar u_{l,t}
#' C_{t,lk}}; the level passes this on to its own earlier values through the
#' persistence and to the earlier driving scores through the loadings. The
#' pass reads the curvature \eqn{C_t} that the forward pass stored and runs
#' compiled.
#' @param term A [GasJoint()].
#' @param eta,y,score,curvature,psi As for [term_filter()].
#' @param g The direct derivative of the objective in the filtered
#'   predictors, a matrix like `eta`.
#' @param ... Unused.
#' @param stored_curv `NULL`, or the `curv` that [term_filter()] returned at
#'   the same parameters and static predictors; where it is `NULL`, the
#'   forward pass is run to obtain it. The name is not a prefix of
#'   `curvature`, so that it cannot be matched partially to that argument.
#' @param fast,threads Passed to that forward pass.
#' @return A list with `deta` and `dscore`, matrices like `eta`.
#' @keywords internal
S7::method(term_adjoint, GasJoint) <- function(term, eta, y, score, curvature,
                                               psi, g, ..., stored_curv = NULL,
                                               fast = NULL, threads = 1L) {
  bp <- term@build_info
  K <- length(term@params)
  n <- bp$n
  g <- matrix(as.numeric(g), n, K)
  cv <- if (is.null(stored_curv)) {
    term_filter(term, eta, y, score, curvature, psi, fast = fast,
                threads = threads)$curv
  } else {
    stored_curv
  }
  if (length(dim(cv)) != 3L || any(dim(cv) != c(n, K, K))) {
    stop(sprintf("'stored_curv' must be an array of dimension c(%d, %d, %d).",
                 n, K, K), call. = FALSE)
  }
  if (.gas_joint_has_sub(term)) {
    out <- .gas_joint_sub_adjoint(term, cv, g, psi)
    colnames(out$deta) <- colnames(out$dscore) <- term@params
    return(out)
  }
  cf <- .gas_joint_start(term, .gas_joint_coefs(term, psi), eta)
  out <- gas_joint_adjoint_cpp(as.numeric(cv), g, cf, bp$order)
  colnames(out$deta) <- colnames(out$dscore) <- term@params
  if (is.null(out$dstart)) out$dstart <- NULL
  else colnames(out$dstart) <- term@params
  out
}

# the per-term pieces of a joint vector keyed as term_params() of the joint
# filter: the entries of term k without their prefix
.gas_joint_split <- function(term, x) {
  nm <- term_params(term)
  out <- list()
  for (k in term@params) {
    own <- term_params(term@terms[[k]])
    out[[k]] <- if (length(own)) stats::setNames(x[paste0(k, ".", own)], own)
                else stats::setNames(x[0], character(0))
  }
  out
}

# the names of the drivers of a joint filter, in their order
.gas_joint_drivers <- function(term) {
  term@params[vapply(term@terms[term@params], .gas_is_driver, logical(1))]
}

#' @title The Charts of a Joint Score-Driven Filter
#' @name term_links.GasJoint
#' @description
#' Returns the links of every term of the filter, named as [term_params()]
#' of the joint filter names the parameters.
#' @param term A [GasJoint()].
#' @param ... Unused.
#' @return A named list of \pkg{linkfunctions7} links.
#' @keywords internal
S7::method(term_links, GasJoint) <- function(term, ...) {
  out <- list()
  for (k in term@params) {
    if (.gas_is_driver(term@terms[[k]])) next
    lk <- term_links(term@terms[[k]])
    names(lk) <- paste0(k, ".", names(lk))
    out <- c(out, lk)
  }
  out
}

#' @title Where a Joint Score-Driven Filter's Parameters Start
#' @name term_start.GasJoint
#' @description
#' Returns the start of every term of the filter, from [term_start()] of the
#' term, named as [term_params()] of the joint filter names the parameters.
#' The cross loadings start at zero, which is the filter without them.
#' @param term A [GasJoint()].
#' @param ... Unused.
#' @return A named numeric vector on the unconstrained scale.
#' @keywords internal
S7::method(term_start, GasJoint) <- function(term, ...) {
  out <- numeric(0)
  for (k in term@params) {
    if (.gas_is_driver(term@terms[[k]])) next
    z <- term_start(term@terms[[k]])
    names(z) <- paste0(k, ".", names(z))
    out <- c(out, z)
  }
  out
}

#' @title The Starting Points of a Joint Score-Driven Filter
#' @name term_starts.GasJoint
#' @description
#' Returns the first start, [term_start()]'s, and then, for each term of the
#' filter in turn, one start for each of that term's further starts
#' ([term_starts()] of the term), with the other terms at their first start.
#' The number of starts is therefore one plus the sum, over the terms, of
#' their starts beyond the first.
#' @param term A [GasJoint()].
#' @param ... Unused.
#' @param target Unused.
#' @return A non-empty list of named numeric vectors on the unconstrained
#'   scale, named as [term_params()].
#' @keywords internal
S7::method(term_starts, GasJoint) <- function(term, ..., target = NULL) {
  z0 <- term_start(term)
  out <- list(z0)
  for (k in term@params) {
    tk <- term@terms[[k]]
    if (.gas_is_driver(tk)) next
    sk <- term_starts(tk)
    for (s in sk[-1L]) {
      z <- z0
      z[paste0(k, ".", names(s))] <- s
      out[[length(out) + 1L]] <- z
    }
  }
  out
}

#' @title The Levels of a Joint Score-Driven Filter
#' @name term_level_param.GasJoint
#' @description
#' Returns the level of every term of the filter, such as `mu.omega` and
#' `sigma.omega`. Each adds a constant to the predictor of its own
#' parameter, so a fitting layer compares each with the design of that
#' parameter's equation.
#' @param term A [GasJoint()].
#' @param ... Unused.
#' @return A character vector with one entry per filtered parameter.
#' @keywords internal
S7::method(term_level_param, GasJoint) <- function(term, ...) {
  .gas_joint_prefixed(term, function(tm) term_level_param(tm))
}

#' @title The Reported Quantities of a Joint Score-Driven Filter
#' @name term_readable.GasJoint
#' @description
#' Returns the quantities of every term of the filter, from
#' [term_readable()] of the term, with the names prefixed by the term's
#' distribution parameter and the Jacobian placed in the term's own columns.
#' @param term A [GasJoint()].
#' @param zeta The parameters on the unconstrained scale, named as
#'   [term_params()] of the joint filter.
#' @param ... Unused.
#' @return A list with `name`, `value`, `jacobian` and `scale`, as
#'   [term_readable()] returns them.
#' @keywords internal
S7::method(term_readable, GasJoint) <- function(term, zeta, ...) {
  nm <- term_params(term)
  z <- unlist(zeta[nm])
  parts <- .gas_joint_split(term, z)
  name <- character(0)
  value <- numeric(0)
  scale <- list()
  rows <- list()
  for (k in term@params) {
    if (.gas_is_driver(term@terms[[k]])) next
    r <- term_readable(term@terms[[k]], as.list(parts[[k]]))
    J <- matrix(0, length(r$name), length(nm), dimnames = list(
      paste0(k, ".", r$name), nm))
    J[, paste0(k, ".", colnames(r$jacobian))] <- r$jacobian
    rows[[k]] <- J
    name <- c(name, paste0(k, ".", r$name))
    value <- c(value, r$value)
    sc <- r$scale
    names(sc) <- paste0(k, ".", names(sc))
    scale <- c(scale, sc)
  }
  list(name = name, value = value,
       jacobian = do.call(rbind, unname(rows)), scale = scale)
}

#' @title Penalties of a Joint Score-Driven Filter
#' @name term_penalties.GasJoint
#' @description
#' Returns the penalties that the subformulas of every term of the filter
#' declare, from [term_penalties()] of the term, each under the key
#' `parameter.key` and with its indices moved to the term's position in
#' [term_params()] of the joint filter. The scalar parameters and the cross
#' loadings are unpenalized, so the list is empty when no term carries a
#' subformula.
#' @param term A built [GasJoint()].
#' @param ... Unused.
#' @return A list of entries, as [term_penalties()] documents.
#' @keywords internal
S7::method(term_penalties, GasJoint) <- function(term, ...) {
  out <- list()
  for (k in seq_along(term@params)) {
    for (e in term_penalties(term@terms[[term@params[k]]])) {
      out[[length(out) + 1L]] <- .entry_relabel(
        e, paste0(term@params[k], ".", e$name), term@offsets[k] + e$index)
    }
  }
  out
}

#' @title How a Joint Score-Driven Filter's Parameters Divide
#' @name term_components.GasJoint
#' @description
#' Returns the entries of [term_components()] of every term of the filter,
#' named and indexed in [term_params()] of the joint filter, followed by one
#' entry per cross loading of the term.
#' @param term A built [GasJoint()].
#' @param ... Unused.
#' @return A named list of entries, as [term_components()] documents.
#' @keywords internal
S7::method(term_components, GasJoint) <- function(term, ...) {
  out <- list()
  for (k in seq_along(term@params)) {
    tm <- term@terms[[term@params[k]]]
    if (.gas_is_driver(tm)) next
    pre <- paste0(term@params[k], ".")
    off <- term@offsets[k]
    cps <- term_components(tm)
    for (cp in cps) {
      cp$name <- paste0(pre, cp$name)
      cp$index <- off + cp$index
      cp$sub_index <- lapply(cp$sub_index, function(i) off + i)
      out[[cp$name]] <- cp
    }
    own <- term_params(tm)
    covered <- unlist(lapply(cps, function(cp) cp$index))
    for (i in setdiff(seq_along(own), covered)) {
      nm <- paste0(pre, own[i])
      out[[nm]] <- list(name = nm, index = off + i, subs = list(),
                        sub_index = list())
    }
  }
  out
}

#' @title Build a Joint Score-Driven Filter Against Data
#' @name term_build.GasJoint
#' @description
#' Builds every term of the filter against `data` and joins them again with
#' [gas_joint()], which checks once more that they walk the same series.
#' @param term A [GasJoint()].
#' @param data A data frame.
#' @param ... Passed to [term_build()] of each term.
#' @return A built [GasJoint()].
#' @keywords internal
S7::method(term_build, GasJoint) <- function(term, data, ...) {
  # the drivers are made again by gas_joint() beside the rebuilt terms
  tms <- lapply(.gas_joint_terms_only(term), function(tm)
    term_build(tm, data, ...))
  out <- gas_joint(tms, drivers = .gas_joint_drivers(term))
  out@label <- term@label
  out@keys <- term@keys
  out
}

#' @title Print a Joint Score-Driven Filter
#' @name print.GasJoint
#' @description
#' Prints the filtered parameters with the orders of each term, the scores
#' that each level loads besides its own with their lags, and the scaling.
#' @param x A [GasJoint()].
#' @param ... Unused.
#' @return `x`, invisibly.
#' @keywords internal
S7::method(print, GasJoint) <- function(x, ...) {
  drv <- .gas_joint_drivers(x)
  cat(sprintf("<GasJoint> '%s': %d filtered parameters, %d group(s)\n",
              x@label, length(x@params) - length(drv),
              length(x@build_info$order)))
  for (k in x@params) {
    tm <- x@terms[[k]]
    if (.gas_is_driver(tm)) {
      cat(sprintf("  %s: static, its score drives the filter\n", k))
      next
    }
    cat(sprintf("  %s: p = %d, q = %d%s\n", k, tm@p, tm@q,
                if (length(tm@cross_with)) sprintf(", loads %s",
                  paste(vapply(tm@cross_with, function(l) {
                    L <- tm@p_cross[[l]]
                    if (L == 1L) sprintf("%s (lag 1)", l)
                    else sprintf("%s (lags 1-%d)", l, L)
                  }, ""), collapse = ", ")) else ""))
  }
  d <- x@terms[[1L]]@scaling
  if (d != 0) cat(sprintf("  scaling: information^%s\n", format(-d)))
  invisible(x)
}
