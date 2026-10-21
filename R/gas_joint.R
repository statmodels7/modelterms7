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
#' \deqn{f_{k,t} = \omega_k + \sum_{i=1}^{p_k} \Big(\kappa_{k,i}\, u_{k,t-i}
#'   + \sum_{l \in C_k} \kappa_{k,i,l}\, u_{l,t-i}\Big)
#'   + \sum_{j=1}^{q_k} \xi_{k,j} f_{k,t-j}, \qquad k \in G,}
#' where \eqn{u_t} is the vector of driving scores of the parameters in
#' \eqn{G}, read at the predictor of every equation at time \eqn{t}, and
#' \eqn{C_k} is the set of the other parameters of the group of `cross` that
#' contains \eqn{k} (empty when the term carries no label). The persistence
#' of each level is its own, so the autoregressive part is stationary when
#' each equation's is.
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
#' @param blueprint The shared ordering: `order` (the rows of each group, in
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
#' returns them as a [GasJoint()] object, with the loadings of `cross` added
#' to the terms that share a label.
#'
#' @details
#' The terms must walk the same series: their `by` and `time` must give the
#' same groups and the same order, because each score is read at the
#' predictors of every equation at the same time. A term whose parameters
#' carry a subformula is not accepted.
#'
#' The terms must share one `scaling`, because the driving vector is scaled
#' by a power of the expected information of all the filtered parameters
#' together. With more than one term the power is 0, 1/2 or 1: for these the
#' derivatives of the matrix power are exact and short (the identity, the
#' symmetric square root through a Sylvester equation per order, and the
#' inverse).
#'
#' The terms that carry the same label in `cross` form a group, in which each
#' level loads the scores of the other parameters of the group as well as its
#' own. A label carried by one term only is rejected, because it describes no
#' group.
#'
#' @param terms A named list of built [GasTerm()] objects, named by the
#'   distribution parameter in whose equation each one sits.
#'
#' @return A [GasJoint()] object. Its `terms` are the input terms with
#'   `cross_with` filled, so their [term_params()] include the cross
#'   loadings.
#'
#' @seealso [GasJoint()], [gas()]
#'
#' @examples
#' set.seed(1)
#' dd <- data.frame(t = 1:40, y = rnorm(40))
#' tm <- list(mu = term_build(gas(p = 1, q = 1, time = t, cross = "a"), dd),
#'            sigma = term_build(gas(p = 1, q = 1, time = t, cross = "a"), dd),
#'            nu = term_build(gas(p = 1, q = 1, time = t), dd))
#' jt <- gas_joint(tm)
#' lapply(jt@terms, term_params)
#'
#' @export
gas_joint <- function(terms) {
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
    if (!length(tm@blueprint)) {
      stop(sprintf("the gas() term of '%s' has not been built.", k),
           call. = FALSE)
    }
    if (!is.null(tm@blueprint$sub)) {
      stop(sprintf(paste0(
        "the gas() term of '%s' develops a parameter with a subformula,\n",
        "  which is not available when several parameters are filtered ",
        "together."), k), call. = FALSE)
    }
  }
  bp1 <- terms[[1L]]@blueprint
  for (k in params[-1L]) {
    bp <- terms[[k]]@blueprint
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
  if (length(terms) > 1L && !(d[[1L]] %in% c(0, 0.5, 1))) {
    stop(sprintf(paste0(
      "'scaling = %s' is available for a single gas() term only; with ",
      "several\n  filtered parameters 'scaling' must be 0, 1/2 or 1."),
      format(d[[1L]])), call. = FALSE)
  }
  labels <- lapply(terms, function(tm) tm@cross)
  lab <- unlist(labels)
  for (a in unique(lab)) {
    if (sum(lab == a) < 2L) {
      stop(sprintf(paste0(
        "the label cross = \"%s\" is carried by one gas() term only;\n",
        "  a label groups the terms whose levels load each other's ",
        "scores."), a), call. = FALSE)
    }
  }
  partners <- vector("list", length(params))
  for (k in seq_along(params)) {
    tm <- terms[[k]]
    if (length(tm@cross)) {
      same <- vapply(labels, function(l) length(l) && l == tm@cross,
                     logical(1))
      partners[[k]] <- setdiff(which(same), k)
    } else {
      partners[[k]] <- integer(0)
    }
    tm@cross_with <- params[partners[[k]]]
    terms[[k]] <- tm
  }
  np <- vapply(terms, function(tm) length(term_params(tm)), integer(1))
  GasJoint(label = "gas_joint", terms = terms, params = params,
           offsets = as.integer(cumsum(np) - np), partners = partners,
           blueprint = list(order = bp1$order, n = bp1$n))
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
  unlist(lapply(term@params, function(k)
    paste0(k, ".", term_params(term@terms[[k]]))), use.names = FALSE)
}

# The parameters as a list with one entry per filtered parameter. A caller
# may hand them as one flat list named as term_params() of the joint filter
# names them (mu.omega, ...), which is how a fitting layer keeps them.
.gas_joint_psi <- function(term, psi) {
  if (all(term@params %in% names(psi))) return(psi)
  flat <- unlist(psi)
  out <- list()
  for (k in term@params) {
    own <- term_params(term@terms[[k]])
    out[[k]] <- as.list(stats::setNames(flat[paste0(k, ".", own)], own))
  }
  out
}

# The coefficients of the joint recursion from the terms' parameters: per
# term, the chart of .gas_coefs() and the cross loadings, with the columns
# of the joint vector at which each coefficient's derivative sits.
.gas_joint_coefs <- function(term, psi) {
  psi <- .gas_joint_psi(term, psi)
  out <- vector("list", length(term@params))
  for (k in seq_along(term@params)) {
    tm <- term@terms[[term@params[k]]]
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
    cl <- matrix(0, length(pk), p)
    ccol <- matrix(0L, length(pk), p)
    for (j in seq_along(pk)) {
      for (i in seq_len(p)) {
        ccol[j, i] <- off + nb + (j - 1L) * p + i
        cl[j, i] <- v[[nb + (j - 1L) * p + i]]
      }
    }
    out[[k]] <- list(p = p, q = q, omega = cf$omega, a = cf$a, b = cf$b,
                     db = cf$db, f0 = cf$f0, df0 = cf$df0,
                     own = off + seq_len(nb), col_omega = off + 1L,
                     col_a = off + cf$i_a, partners = pk, cross = cl,
                     col_cross = ccol)
  }
  out
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
#'   where they do not, the context is inert and the callbacks run in R.
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
  bp <- term@blueprint
  K <- length(term@params)
  n <- bp$n
  eta <- matrix(as.numeric(eta), n, K)
  cf <- .gas_joint_coefs(term, psi)
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
  for (rows in bp$order) {
    T_ <- length(rows)
    f <- matrix(0, K, T_)
    u <- matrix(0, K, T_)
    Fd <- array(0, c(K, T_, M))
    Sd <- array(0, c(K, T_, M))
    for (t in seq_len(T_)) {
      for (k in seq_len(K)) {
        c_k <- cf[[k]]
        ft <- c_k$omega
        Ft <- numeric(M)
        Ft[c_k$col_omega] <- 1
        for (i in seq_len(c_k$p)) {
          if (t - i < 1L) next
          ft <- ft + c_k$a[[i]] * u[k, t - i]
          Ft <- Ft + c_k$a[[i]] * Sd[k, t - i, ]
          Ft[c_k$col_a[i]] <- Ft[c_k$col_a[i]] + u[k, t - i]
          for (j in seq_along(c_k$partners)) {
            l <- c_k$partners[j]
            ft <- ft + c_k$cross[j, i] * u[l, t - i]
            Ft <- Ft + c_k$cross[j, i] * Sd[l, t - i, ]
            cc <- c_k$col_cross[j, i]
            Ft[cc] <- Ft[cc] + u[l, t - i]
          }
        }
        for (j in seq_len(c_k$q)) {
          if (t - j >= 1L) {
            f_l <- f[k, t - j]
            F_l <- Fd[k, t - j, ]
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
  bp <- term@blueprint
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
  out <- gas_joint_adjoint_cpp(as.numeric(cv), g, .gas_joint_coefs(term, psi),
                               bp$order)
  colnames(out$deta) <- colnames(out$dscore) <- term@params
  out
}

# the per-term pieces of a joint vector keyed as term_params() of the joint
# filter: the entries of term k without their prefix
.gas_joint_split <- function(term, x) {
  nm <- term_params(term)
  out <- list()
  for (k in term@params) {
    own <- term_params(term@terms[[k]])
    out[[k]] <- stats::setNames(x[paste0(k, ".", own)], own)
  }
  out
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
    z <- term_start(term@terms[[k]])
    names(z) <- paste0(k, ".", names(z))
    out <- c(out, z)
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
  unlist(lapply(term@params, function(k)
    paste0(k, ".", term_level_param(term@terms[[k]]))), use.names = FALSE)
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
#' Returns an empty list: the terms of a joint filter carry no subformula,
#' and their scalar parameters are unpenalized.
#' @param term A [GasJoint()].
#' @param ... Unused.
#' @return An empty list.
#' @keywords internal
S7::method(term_penalties, GasJoint) <- function(term, ...) list()

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
  tms <- lapply(term@terms, function(tm) term_build(tm, data, ...))
  out <- gas_joint(tms)
  out@label <- term@label
  out@keys <- term@keys
  out
}

#' @title Print a Joint Score-Driven Filter
#' @name print.GasJoint
#' @description
#' Prints the filtered parameters with the orders of each term, the groups
#' of cross loadings, and the scaling.
#' @param x A [GasJoint()].
#' @param ... Unused.
#' @return `x`, invisibly.
#' @keywords internal
S7::method(print, GasJoint) <- function(x, ...) {
  cat(sprintf("<GasJoint> '%s': %d filtered parameters, %d group(s)\n",
              x@label, length(x@params), length(x@blueprint$order)))
  for (k in x@params) {
    tm <- x@terms[[k]]
    cat(sprintf("  %s: p = %d, q = %d%s\n", k, tm@p, tm@q,
                if (length(tm@cross_with)) sprintf(", loads %s",
                  paste(tm@cross_with, collapse = ", ")) else ""))
  }
  d <- x@terms[[1L]]@scaling
  if (d != 0) cat(sprintf("  scaling: information^%s\n", format(-d)))
  invisible(x)
}
