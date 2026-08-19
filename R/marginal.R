#' @include term_classes.R generics.R structural.R regime.R segmented.R
NULL

#' @title The Levels of a Likelihood-Shaped Structural Term
#'
#' @description
#' The shifts a term of the likelihood shape adds to its equation's
#' predictor, one per mixture component, in the order the columns of
#' \code{\link{term_posterior}} carry the components.
#'
#' @details
#' By Fisher's identity the derivative of a likelihood mixed over latent
#' states, in any predictor the model carries, is the posterior-weighted
#' derivative of the ordinary one, each component read at the predictor
#' shifted by its own level. \code{\link{term_posterior}} supplies the
#' weights; this supplies the levels, so a fitting layer assembles the
#' identity without reading the term's internals. For \code{\link{regime}}
#' the levels are the ordered regime means, one number per component; for a
#' marginal break-point term of the step kind they are the sums of the
#' changes of level over the active break-points, one number per side
#' pattern.
#'
#' A component's shift may vary by observation -- the quadrature nodes of a
#' marginal \code{\link{seg}} or \code{\link{jseg}} term shift each
#' observation by its own hinge value -- and the method then returns a
#' matrix with one row per observation and one column per component, whose
#' columns a caller reads in place of the constant levels.
#'
#' @param term A built structural term of the likelihood shape.
#' @param psi The term's parameters, named as \code{\link{term_params}}.
#' @param ... Passed to methods.
#'
#' @return A numeric vector with one level per component, or a matrix with
#'   one row per observation and one column per component.
#'
#' @examples
#' term_levels(regime(2), list(level1 = 0, gap2 = 3,
#'                             alr1.1 = 2, alr2.1 = -2))
#'
#' @seealso \code{\link{term_posterior}}, \code{\link{term_loglik}}
#' @export
term_levels <- S7::new_generic("term_levels", "term",
  function(term, psi, ...) S7::S7_dispatch())

S7::method(term_levels, structural_term) <- function(term, psi, ...) {
  stop(sprintf("the term class '%s' does not implement term_levels().",
               attr(S7::S7_class(term), "name")), call. = FALSE)
}

# Registered here rather than in regime.R because the generic is: this file
# loads after that one, so the class exists when the method is written.
#' @title Levels of a Regime Term
#' @name term_levels.RegimeTerm
#' @description
#' The regime means, ordered by construction: the first level and the
#' cumulative sums of the positive gaps.
#' @param term A \code{\link{RegimeTerm}}.
#' @param psi The term's parameters.
#' @param ... Unused.
#' @return A numeric vector of length \code{k}.
#' @keywords internal
S7::method(term_levels, RegimeTerm) <- function(term, psi, ...) {
  v <- unlist(psi[term_params(term)])
  k <- term@k
  unname(cumsum(c(v[["level1"]],
                  if (k > 1L) v[paste0("gap", seq.int(2L, k))]
                  else numeric(0))))
}


#' @title The Posterior of a Structural Term's Latent Variable
#'
#' @description
#' A summary of the latent variable a structural term integrates over,
#' given the whole sample: for a marginal break-point term, the posterior
#' mean and standard deviation of each group's break-points.
#'
#' @details
#' \code{\link{term_posterior}} answers the fitting layer's question, the
#' component weights Fisher's identity needs at every observation. This one
#' answers the reader's: where each group's latent positions sit once the
#' data have been seen. For the marginal break-point term the two come from
#' the same decomposition; the mean and variance within an interval are
#' those of the prior truncated to it, and under quadrature the moments of
#' the node posterior.
#'
#' @param term A built structural term.
#' @param eta The static predictor of the equation the term sits in.
#' @param y The response.
#' @param logdens The log-density as a function of a predictor value and a
#'   row index, as \code{\link{term_loglik}} takes it.
#' @param psi The term's parameters, named as \code{\link{term_params}}.
#' @param ... Passed to methods.
#'
#' @return A data frame with one row per group and break-point:
#'   \code{group}, \code{psi} (which break-point), \code{mean} and
#'   \code{sd}.
#'
#' @examples
#' set.seed(1)
#' dd <- data.frame(id = rep(1:3, each = 8), x = rep(1:8, 3))
#' dd$y <- rnorm(24, 2 * (dd$x >= 4.5), 0.4)
#' tm <- term_build(jump(x, psi ~ random(~1 | id), marginal = TRUE), dd)
#' term_latent(tm, rep(0, 24), dd$y,
#'             logdens = function(e, i) dnorm(dd$y[i], e, 0.4, log = TRUE),
#'             psi = list(m1 = 4.5, tau1 = 0.5, delta1 = 2))
#'
#' @seealso \code{\link{term_posterior}}, \code{\link{jump}}
#' @export
term_latent <- S7::new_generic("term_latent", "term",
  function(term, eta, y, logdens, psi, ...) S7::S7_dispatch())

S7::method(term_latent, structural_term) <- function(term, eta, y, logdens,
                                                     psi, ...) {
  stop(sprintf("the term class '%s' does not implement term_latent().",
               attr(S7::S7_class(term), "name")), call. = FALSE)
}


#' @title S7 Class for Marginal Break-Point Terms
#' @name MarginalBreakTerm
#'
#' @description
#' A subclass of \code{\link{structural_term}} for break-points that vary
#' by group as latent variables integrated out of the likelihood.
#' Constructed by \code{\link{jump}}, \code{\link{seg}} or
#' \code{\link{jseg}} with \code{marginal = TRUE}.
#'
#' @inheritParams model_term
#' @param kind Which of the three constructions.
#' @param var The covariate expression.
#' @param npsi The number of break-points.
#' @param linear Whether the term carries the linear effect as its own
#'   parameter (\code{seg} and \code{jseg}).
#' @param group The grouping expression, from the break-point's
#'   \code{random()} subformula.
#' @param prior The latent's distribution: \code{NULL} for the gaussian, or
#'   a \pkg{distributions7} object from \code{random(distrib = )}.
#' @param spec The resolved construction settings.
#' @param blueprint The resolved grouping and interval structure.
#'
#' @return An object of class \code{MarginalBreakTerm}.
#'
#' @seealso \code{\link{jump}}
#' @examples
#' S7::S7_inherits(jump(x, psi ~ random(~1 | id), marginal = TRUE),
#'                 MarginalBreakTerm)
#' @export
MarginalBreakTerm <- S7::new_class(
  name = "MarginalBreakTerm",
  parent = structural_term,
  properties = list(
    kind = S7::class_character,
    var = S7::class_any,
    npsi = S7::class_integer,
    linear = S7::class_logical,
    group = S7::class_any,
    prior = S7::class_any,
    spec = S7::class_list,
    blueprint = S7::class_list
  )
)

# The constructor branch of the three break-point terms under
# marginal = TRUE. The prior is read off the break-point's subformula,
# which must be a single random(~1 | g) call: the latent IS the random
# effect, and under the marginal its distribution is part of the
# likelihood rather than a penalty.
.marg_spec <- function(kind, var, npsi, psi, by, dots, smoothed, c0_set,
                       n_boot_set, label, linear = FALSE) {
  if (!is.null(smoothed)) {
    message(paste("'smoothed' is ignored under marginal = TRUE: the marginal",
                  "likelihood is smooth in every parameter already, and a",
                  "mollifier would only add its bias."))
  }
  if (c0_set) {
    message(paste("'c0' is the scaling schedule of the working",
                  "parametrization, and the marginal construction has none:",
                  "it is ignored."))
  }
  if (n_boot_set) {
    message(paste("'n_boot' is the restart count of the working",
                  "construction; the marginal term starts from the exact",
                  "profile and declares none."))
  }
  if (!is.null(by)) {
    stop(paste("'by' does not combine with marginal = TRUE: the grouping",
               "comes from the break-point's subformula, psi ~ random(~1 |",
               "g), and the term's own coefficients are scalars."),
         call. = FALSE)
  }
  if (!is.numeric(npsi) || length(npsi) != 1L || is.na(npsi) ||
      npsi != round(npsi) || npsi < 1) {
    stop("'npsi' must be a single integer of at least 1.", call. = FALSE)
  }
  npsi <- as.integer(npsi)
  if (kind == "jump" && npsi > 3L) {
    stop(paste("marginal = TRUE covers up to three break-points: the exact",
               "sum runs over the product partition of the intervals, whose",
               "(n+1)^K cells grow past what one evaluation can afford."),
         call. = FALSE)
  }
  if (kind != "jump" && npsi > 1L) {
    stop(sprintf(paste("marginal = TRUE covers one break-point for a %s",
                       "term: its conditional is smooth in the positions, so",
                       "several latents need a product quadrature whose",
                       "component count no fitting layer can carry."), kind),
         call. = FALSE)
  }
  if (!is.null(psi) && (!is.numeric(psi) || length(psi) != npsi ||
                        anyNA(psi))) {
    stop("'psi' must give one starting position per break-point.",
         call. = FALSE)
  }
  if (!is.character(label) || length(label) != 1L || is.na(label) ||
      !nzchar(label)) {
    stop("'label' must be a single non-empty character string.",
         call. = FALSE)
  }

  need <- paste("marginal = TRUE requires the break-point to carry the",
                "latent, psi ~ random(~1 | g): the marginal likelihood",
                "integrates over it, and without one there is nothing to",
                "integrate.")
  if (!length(dots)) stop(need, call. = FALSE)
  nms <- names(dots)
  sub <- NULL
  for (i in seq_along(dots)) {
    d <- dots[[i]]
    if (!is.null(nms) && nzchar(nms[i])) {
      stop(sprintf("'%s' is not an argument of a break-point term.", nms[i]),
           call. = FALSE)
    }
    if (!inherits(d, "formula") || length(d) != 3L || !is.name(d[[2L]])) {
      stop(paste("arguments in '...' must be two-sided formulas developing",
                 "a coefficient of the term, e.g. psi ~ random(~1 | g)."),
           call. = FALSE)
    }
    lhs <- as.character(d[[2L]])
    if (!lhs %in% c("psi", paste0("psi", seq_len(npsi)))) {
      stop(sprintf(paste("'%s' cannot be developed under marginal = TRUE:",
                         "the term's own coefficients are scalars, estimated",
                         "by maximum likelihood with the prior's",
                         "parameters."), lhs), call. = FALSE)
    }
    if (!identical(lhs, "psi") && npsi > 1L) {
      stop(paste("under marginal = TRUE the latent development is shared:",
                 "name the stem, psi ~ random(~1 | g), so that every",
                 "break-point carries it. A fixed break-point beside a",
                 "random one would put the step back into an ordinary",
                 "parameter."), call. = FALSE)
    }
    if (!is.null(sub)) stop("'psi' is developed twice.", call. = FALSE)
    sub <- d
  }
  if (is.null(sub)) stop(need, call. = FALSE)

  rt <- tryCatch(eval(sub[[3L]], environment(sub)), error = function(e) NULL)
  if (is.null(rt) || !S7::S7_inherits(rt, RandomTerm)) {
    stop(paste("under marginal = TRUE the break-point's subformula must be a",
               "single random(...) call: the latent IS the random effect,",
               "and its distribution is the prior the likelihood integrates",
               "against."), call. = FALSE)
  }
  e <- rt@formula[[2L]]
  if (!(is.numeric(e[[2L]]) && length(e[[2L]]) == 1L && e[[2L]] == 1)) {
    stop(paste("the latent break-point is one value per group: the random()",
               "of a marginal term takes an intercept-only formula,",
               "~ 1 | g."), call. = FALSE)
  }
  prior <- rt@distrib
  if (!is.null(prior)) {
    if (kind != "jump" || npsi > 1L) {
      stop(paste("a prior other than the gaussian is carried by the single",
                 "break-point of a jump term, whose interval masses are",
                 "differences of the prior's cdf; the quadrature",
                 "constructions read the prior's density and are built for",
                 "the gaussian."), call. = FALSE)
    }
    if (!S7::S7_inherits(prior, distributions7::continuous_distrib)) {
      stop("the latent's prior must be a continuous univariate distribution.",
           call. = FALSE)
    }
    ip <- tryCatch(prior@params_interpretation, error = function(e) NULL)
    loc <- names(ip)[unlist(ip) %in% c("location", "mean")]
    if (any(loc %in% prior@params)) {
      stop(sprintf(paste("the prior's location ('%s') must be fixed at zero",
                         "-- fixed(%s, %s = 0) -- since the population",
                         "position m1 carries it."),
                   loc[1L], "student_t1_distrib()", loc[1L]), call. = FALSE)
    }
  }
  if (length(rt@hyper)) {
    stop(paste("'hyper' has nothing to hold under marginal = TRUE: the prior",
               "is part of the likelihood and its parameters are ordinary",
               "ones estimated with the rest."), call. = FALSE)
  }

  MarginalBreakTerm(label = label, kind = kind, var = var, npsi = npsi,
                    linear = linear, group = e[[3L]], prior = prior,
                    spec = list(psi = psi), blueprint = list())
}

# The names of the prior's own parameters for one break-point: the
# gaussian's are the position and the scale, named by the break-point; an
# explicit prior contributes its own free names beside the position.
.marg_prior_names <- function(term, k) {
  if (is.null(term@prior)) c(paste0("m", k), paste0("tau", k))
  else c("m1", term@prior@params)
}

S7::method(term_params, MarginalBreakTerm) <- function(term, ...) {
  K <- term@npsi
  c(if (term@linear) "beta" else character(0),
    unlist(lapply(seq_len(K), function(k) .marg_prior_names(term, k))),
    if (term@kind %in% c("seg", "jseg")) paste0("gamma", seq_len(K))
    else character(0),
    if (term@kind %in% c("jump", "jseg")) paste0("delta", seq_len(K))
    else character(0))
}

S7::method(term_links, MarginalBreakTerm) <- function(term, ...) {
  nm <- term_params(term)
  out <- stats::setNames(vector("list", length(nm)), nm)
  for (p in nm) {
    out[[p]] <- if (grepl("^tau[0-9]+$", p)) {
      linkfunctions7::log_link()
    } else if (!is.null(term@prior) && p %in% term@prior@params) {
      term@prior@link_params[[p]]
    } else {
      linkfunctions7::identity_link()
    }
  }
  out
}

S7::method(term_build, MarginalBreakTerm) <- function(term, data, ...) {
  xv <- .seg_x(term@var, data)
  n <- length(xv)
  if (diff(range(xv)) <= 0) {
    stop("the covariate of a break-point term must vary.", call. = FALSE)
  }
  g <- eval(term@group, data, baseenv())
  if (length(g) != n) {
    stop(paste("the grouping of a marginal break-point term must give one",
               "value per row."), call. = FALSE)
  }
  if (anyNA(g)) {
    stop(paste("the grouping of a marginal break-point term must not contain",
               "missing values."), call. = FALSE)
  }
  f <- factor(g)
  rows <- split(seq_len(n), f)
  term@blueprint <- list(
    n = n, x = xv,
    # each group's rows sorted by the covariate: the intervals between them
    # are where the sum (or the quadrature) runs
    groups = lapply(rows, function(r) r[order(xv[r])]),
    labels = names(rows))
  term
}

# log(Phi(zu) - Phi(zl)), stable in either tail: the difference is taken on
# whichever side keeps both terms away from one, through the log-scale tail
# probabilities.
.marg_lmass <- function(zl, zu) {
  out <- numeric(length(zl))
  lo_inf <- !is.finite(zl) & zl < 0
  hi_inf <- !is.finite(zu) & zu > 0
  i <- lo_inf & hi_inf
  out[i] <- 0
  i <- lo_inf & !hi_inf
  out[i] <- stats::pnorm(zu[i], log.p = TRUE)
  i <- hi_inf & !lo_inf
  out[i] <- stats::pnorm(zl[i], lower.tail = FALSE, log.p = TRUE)
  i <- !lo_inf & !hi_inf
  if (any(i)) {
    a <- zl[i]
    b <- zu[i]
    lower <- (a + b) < 0
    la <- ifelse(lower, stats::pnorm(b, log.p = TRUE),
                 stats::pnorm(a, lower.tail = FALSE, log.p = TRUE))
    lb <- ifelse(lower, stats::pnorm(a, log.p = TRUE),
                 stats::pnorm(b, lower.tail = FALSE, log.p = TRUE))
    d <- lb - la
    out[i] <- la + ifelse(d > -log(2), log(-expm1(d)), log1p(-exp(d)))
  }
  out
}

.marg_lse <- function(a) {
  mx <- max(a)
  if (!is.finite(mx)) return(-Inf)
  mx + log(sum(exp(a - mx)))
}

# The interval decomposition of one group under the GAUSSIAN prior:
# boundaries between its sorted covariate values, the log masses and the
# derivatives of those masses in (m, tau) on the parameter scale. Interval
# j (0-based) is (x_(j), x_(j+1)] with x_(0) = -Inf and x_(n+1) = +Inf: a
# break-point in interval j leaves the observations of sorted rank above j
# on the shifted side. An interval whose mass underflows carries zero
# posterior weight, and its derivative rows are set to zero so that the
# products they enter stay finite.
.marg_intervals <- function(xs, m, tau, d2 = FALSE) {
  b <- c(-Inf, xs, Inf)
  nb <- length(b)
  zl <- (b[-nb] - m) / tau
  zu <- (b[-1L] - m) / tau
  lm <- .marg_lmass(zl, zu)
  mass <- exp(lm)
  ok <- mass > 0
  fd <- function(z) {
    phi <- ifelse(is.finite(z), stats::dnorm(z), 0)
    zf <- ifelse(is.finite(z), z, 0)
    list(phi = phi, zp = zf * phi, z2p = zf^2 * phi, z3p = zf^3 * phi)
  }
  L <- fd(zl)
  U <- fd(zu)
  # dF/dm = -phi(z)/tau and dF/dtau = -z phi(z)/tau at each boundary; a mass
  # derivative is the difference of the two boundary values
  dm <- ifelse(ok, ((L$phi - U$phi) / tau) / mass, 0)
  dt <- ifelse(ok, ((L$zp - U$zp) / tau) / mass, 0)
  out <- list(lm = lm, mass = mass, ok = ok, dm = dm, dt = dt,
              zl = zl, zu = zu, phi_l = L$phi, phi_u = U$phi,
              zp_l = L$zp, zp_u = U$zp)
  if (d2) {
    # second derivatives of F = Phi((b - m)/tau):
    #   F_mm = -z phi/tau^2, F_mtau = phi(1 - z^2)/tau^2,
    #   F_tautau = z phi (2 - z^2)/tau^2
    mm <- (L$zp - U$zp) / tau^2
    mt <- ((U$phi - U$z2p) - (L$phi - L$z2p)) / tau^2
    tt <- ((2 * U$zp - U$z3p) - (2 * L$zp - L$z3p)) / tau^2
    out$dmm <- ifelse(ok, mm / mass, 0) - dm * dm
    out$dmt <- ifelse(ok, mt / mass, 0) - dm * dt
    out$dtt <- ifelse(ok, tt / mass, 0) - dt * dt
  }
  out
}

# The same decomposition under an EXPLICIT prior (one break-point): masses
# are differences of the prior's cdf at the boundaries shifted by the
# position, and their derivatives come from the density (for the position)
# and from the cdf surface (for the prior's own parameters), the route the
# censored likelihoods use. Columns of dlm follow .marg_prior_names().
.marg_prior_intervals <- function(prior, xs, m, theta) {
  b <- c(-Inf, xs, Inf)
  nb <- length(b)
  l <- b[-nb]
  u <- b[-1L]
  th <- stats::setNames(as.list(theta[prior@params]), prior@params)
  cdf_at <- function(q) {
    out <- numeric(length(q))
    fin <- is.finite(q)
    if (any(fin)) {
      out[fin] <- as.numeric(distributions7::distrib_cdf(prior, q[fin] - m,
                                                         th))
    }
    out[!fin] <- ifelse(q[!fin] > 0, 1, 0)
    out
  }
  pdf_at <- function(q) {
    out <- numeric(length(q))
    fin <- is.finite(q)
    if (any(fin)) {
      out[fin] <- as.numeric(distributions7::distrib_pdf(prior, q[fin] - m,
                                                         th))
    }
    out
  }
  gcdf_at <- function(q) {
    out <- matrix(0, length(q), length(prior@params))
    fin <- is.finite(q)
    if (any(fin)) {
      g <- distributions7::distrib_grad_cdf(prior, q[fin] - m, th,
                                            lower.tail = TRUE, log = FALSE)
      for (j in seq_along(prior@params)) {
        out[fin, j] <- rep_len(as.numeric(g[[prior@params[j]]]), sum(fin))
      }
    }
    out
  }
  Fl <- cdf_at(l)
  Fu <- cdf_at(u)
  mass <- pmax(Fu - Fl, 0)
  ok <- mass > 0
  lm <- ifelse(ok, log(mass), -Inf)
  # d/dm F(b - m) = -f(b - m), so the position's column is the density
  # difference; the prior's own columns are the cdf-gradient differences
  dm_col <- ifelse(ok, (pdf_at(l) - pdf_at(u)) / mass, 0)
  Gd <- gcdf_at(u) - gcdf_at(l)
  Gd[!ok, ] <- 0
  Gd[ok, ] <- Gd[ok, , drop = FALSE] / mass[ok]
  dlm <- cbind(dm_col, Gd)
  colnames(dlm) <- c("m1", prior@params)
  list(lm = lm, mass = mass, ok = ok, dlm = dlm, l = l, u = u)
}

.marg_check_psi <- function(term, psi) {
  nm <- term_params(term)
  v <- unlist(psi[nm])
  if (length(v) != length(nm) || anyNA(v)) {
    stop(sprintf("'psi' must supply %s.", paste(nm, collapse = ", ")),
         call. = FALSE)
  }
  taus <- grep("^tau[0-9]+$", nm, value = TRUE)
  if (length(taus) && any(v[taus] <= 0)) {
    stop("every prior scale must be positive.", call. = FALSE)
  }
  v
}

.marg_built <- function(term) {
  bp <- term@blueprint
  if (!length(bp)) {
    stop("the term has not been built; call term_build(term, data) first.",
         call. = FALSE)
  }
  bp
}

# ---------------------------------------------------------------------------
# The step kind: cells. With K latent positions the conditional is constant
# on the product partition of each coordinate's intervals, so the marginal
# is an exact sum over (n+1)^K cells per group. The cells live in a K-array;
# one observation of sorted rank r is on the shifted side of break-point k
# exactly in the cells whose k-th index is at most r.

# sum_k vs[[k]][j_k] over the cell array
.marg_outer_sum <- function(vs) {
  A <- vs[[1L]]
  if (length(vs) > 1L) for (k in 2L:length(vs)) A <- outer(A, vs[[k]], "+")
  A
}

# the array of side patterns of one observation: entry sum_k 2^(k-1) *
# 1(j_k <= r), used to index the 2^K per-pattern values
.marg_pattern <- function(r, J, K) {
  .marg_outer_sum(lapply(seq_len(K), function(k)
    2^(k - 1L) * as.numeric(seq_len(J) <= r)))
}

# the indicator of one coordinate's side, broadcast over the cells
.marg_bit <- function(r, J, K, k) {
  vs <- rep(list(numeric(J)), K)
  vs[[k]] <- as.numeric(seq_len(J) <= r)
  .marg_outer_sum(vs)
}

# marginal sums of a cell array along one coordinate
.marg_margin <- function(w, k) {
  if (is.null(dim(w)) || length(dim(w)) == 1L) return(w)
  apply(w, k, sum)
}

# the per-coordinate interval structures and the prior-parameter columns
# they own, shared by every method of the step kind
.marg_jump_prior <- function(term, xs, v) {
  K <- term@npsi
  if (is.null(term@prior)) {
    ivs <- lapply(seq_len(K), function(k)
      .marg_intervals(xs, v[[paste0("m", k)]], v[[paste0("tau", k)]]))
    list(
      lm = lapply(ivs, `[[`, "lm"),
      dlm = lapply(seq_len(K), function(k) {
        M <- cbind(ivs[[k]]$dm, ivs[[k]]$dt)
        M[!ivs[[k]]$ok, ] <- 0
        colnames(M) <- c(paste0("m", k), paste0("tau", k))
        M
      }),
      ivs = ivs)
  } else {
    pv <- .marg_prior_intervals(term@prior, xs, v[["m1"]],
                                v[term@prior@params])
    list(lm = list(pv$lm), dlm = list(pv$dlm), ivs = list(pv))
  }
}

# the 2^K side patterns: their bit matrix and the shift each adds
.marg_bits <- function(K, deltas) {
  P <- 2^K
  bits <- matrix(0L, P, K)
  for (k in seq_len(K)) bits[, k] <- bitwAnd(seq_len(P) - 1L, 2^(k - 1L)) > 0
  list(bits = bits, shifts = as.numeric(bits %*% deltas))
}

#' @title Log-Likelihood of a Marginal Break-Point Term
#' @name term_loglik.MarginalBreakTerm
#' @description
#' The exact marginal likelihood of a break-point model with latent
#' positions per group, with the exact derivatives in the term's own
#' parameters on the parameter scale.
#' @details
#' For the step kind the conditional is constant on the product partition
#' of the intervals between a group's ordered observations, so the
#' integral over the latents is a finite sum: masses are differences of
#' the prior's cdf and the conditional is updated by one density ratio per
#' cell, one observation changing side with respect to one break-point.
#' For the continuous kinds the conditional is smooth within an interval
#' and the integral runs on a fixed Gauss-Kronrod panel per interval
#' (\code{\link[numericals7]{gauss_kronrod15}}), the interior nodes fixed
#' points of the data so that the derivatives in the prior's parameters
#' read the prior alone; the region below the data, where the hinge keeps
#' moving, is covered by panels that follow the prior's bulk, whose node
#' motion the derivatives carry, and the region above it, where the
#' conditional is constant, contributes its closed tail mass.
#'
#' The per-observation contributions are the one-step predictive densities
#' given the group's observations at smaller covariate values, which sum
#' to each group's marginal log-likelihood; producing them costs one pass
#' over the component weights per observation, so the decomposition is a
#' factor of the group's size dearer than the total alone, which the
#' per-group cell or node count already prices.
#' @param term A built \code{\link{MarginalBreakTerm}}.
#' @param eta The static predictor.
#' @param y The response, reaching the sum through the callbacks.
#' @param logdens,score The log-density and its derivative in the predictor.
#' @param psi The parameters, named as \code{\link{term_params}}.
#' @param ... Unused.
#' @return A list with \code{loglik} and \code{jacobian}, the latter on the
#'   parameter scale.
#' @keywords internal
S7::method(term_loglik, MarginalBreakTerm) <- function(term, eta, y, logdens,
                                                       score, psi, ...) {
  if (term@kind == "jump") {
    .marg_jump_loglik(term, eta, y, logdens, score, psi)
  } else {
    .marg_seg_loglik(term, eta, y, logdens, score, psi)
  }
}

.marg_jump_loglik <- function(term, eta, y, logdens, score, psi) {
  bp <- .marg_built(term)
  v <- .marg_check_psi(term, psi)
  nm <- term_params(term)
  K <- term@npsi
  n <- bp$n
  idx <- seq_len(n)
  pb <- .marg_bits(K, v[paste0("delta", seq_len(K))])
  P <- 2^K
  LD <- matrix(0, n, P)
  SC <- matrix(0, n, P)
  for (p in seq_len(P)) {
    LD[, p] <- as.numeric(logdens(eta + pb$shifts[p], idx))
    SC[, p] <- as.numeric(score(eta + pb$shifts[p], idx))
  }
  if (any(!is.finite(LD[, 1L])) && all(is.finite(y))) {
    # a non-finite log-density is the callback's to explain; the sum only
    # propagates it
  }

  ll <- numeric(n)
  jac <- matrix(0, n, length(nm), dimnames = list(NULL, nm))
  for (rs in bp$groups) {
    ng <- length(rs)
    J <- ng + 1L
    pr <- .marg_jump_prior(term, bp$x[rs], v)
    A <- .marg_outer_sum(pr$lm)
    tot <- .marg_lse(A)
    u <- exp(A - tot)
    Dk <- rep(list(A * 0), K)
    for (t in seq_len(ng)) {
      row <- rs[t]
      PA1 <- .marg_pattern(t, J, K) + 1
      lf <- LD[row, ][PA1]
      scv <- SC[row, ][PA1]
      if (K > 1L) {
        dim(lf) <- dim(A)
        dim(scv) <- dim(A)
      }
      A2 <- A + lf
      tot2 <- .marg_lse(A2)
      w <- exp(A2 - tot2)
      ll[row] <- tot2 - tot
      for (k in seq_len(K)) {
        bit <- .marg_bit(t, J, K, k)
        D2 <- Dk[[k]] + scv * bit
        jac[row, paste0("delta", k)] <- sum(w * D2) - sum(u * Dk[[k]])
        Dk[[k]] <- D2
      }
      for (k in seq_along(pr$dlm)) {
        wm <- .marg_margin(w, k) - .marg_margin(u, k)
        jac[row, colnames(pr$dlm[[k]])] <-
          jac[row, colnames(pr$dlm[[k]])] + as.numeric(wm %*% pr$dlm[[k]])
      }
      A <- A2
      tot <- tot2
      u <- w
    }
  }
  list(loglik = ll, jacobian = jac)
}

# ---------------------------------------------------------------------------
# The continuous kinds: nodes. The conditional is smooth within an interval,
# so each interval carries a fixed Gauss-Kronrod panel; the interval
# boundaries are the data, so interior nodes do not move with the prior's
# parameters. Below the data the hinge keeps varying, and the panels there
# follow the prior's bulk: their nodes move with (m, tau) and the
# derivatives carry that motion. Above the data every device is saturated
# at zero, the conditional is constant, and the tail contributes its
# closed mass through a node at +Inf.

# the node set of one group: positions, log quadrature weights (prior
# density included), and the geometry the derivatives need
.marg_nodes <- function(xs, m, tau) {
  gk <- numericals7::gauss_kronrod15()
  x1 <- xs[1L]
  xn <- xs[length(xs)]
  lo_pts <- numeric(0)
  hi_pts <- numeric(0)
  edge <- logical(0)
  # the lower region, from 8.5 prior sds below the bulk up to the first
  # observation; its panel edges are affine in the lower limit, which is
  # what the node-motion derivatives read
  a <- min(m, x1) - 8.5 * tau
  wid <- x1 - a
  ns <- min(8L, max(4L, ceiling(wid / (2.5 * tau))))
  ed <- seq(a, x1, length.out = ns + 1L)
  for (s in seq_len(ns)) {
    lo_pts <- c(lo_pts, ed[s])
    hi_pts <- c(hi_pts, ed[s + 1L])
    edge <- c(edge, TRUE)
  }
  # the interior intervals, split so that no panel spans more than a couple
  # of prior sds
  for (j in seq_len(length(xs) - 1L)) {
    wj <- xs[j + 1L] - xs[j]
    if (wj <= 0) next
    nsj <- min(6L, max(1L, ceiling(wj / (2.5 * tau))))
    edj <- seq(xs[j], xs[j + 1L], length.out = nsj + 1L)
    for (s in seq_len(nsj)) {
      lo_pts <- c(lo_pts, edj[s])
      hi_pts <- c(hi_pts, edj[s + 1L])
      edge <- c(edge, FALSE)
    }
  }
  nn <- length(gk$nodes)
  np <- length(lo_pts)
  h <- rep((hi_pts - lo_pts) / 2, each = nn)
  mid <- rep((hi_pts + lo_pts) / 2, each = nn)
  p <- mid + h * gk$nodes
  wk <- rep(gk$wk, np)
  z <- (p - m) / tau
  lw <- log(wk * h) + stats::dnorm(z, log = TRUE) - log(tau)
  is_edge <- rep(edge, each = nn)
  # the relative position of a node inside the lower region: its motion is
  # (1 - t) times the motion of the lower limit
  tfrac <- ifelse(is_edge, (p - a) / pmax(x1 - a, .Machine$double.eps), 0)
  # the closed upper tail: every device saturates at zero above the data,
  # so the conditional there is the one a node at +Inf evaluates
  zn <- (xn - m) / tau
  mr <- numericals7::mills_ratio(-zn)$r
  p <- c(p, Inf)
  lw <- c(lw, stats::pnorm(zn, lower.tail = FALSE, log.p = TRUE))
  z <- c(z, Inf)
  is_edge <- c(is_edge, FALSE)
  tfrac <- c(tfrac, 0)
  # d log(weight)/d(m, tau) at fixed node: the prior's own derivatives, the
  # tail's through the Mills ratio; the log panel width of the lower region
  # scales with the limit as well
  glw_m <- c(z[-length(p)] / tau, mr / tau)
  glw_t <- c((z[-length(p)]^2 - 1) / tau, zn * mr / tau)
  da_m <- as.numeric(m < x1)
  da_t <- -8.5
  glw_m <- glw_m + is_edge * (-da_m / pmax(x1 - a, .Machine$double.eps))
  glw_t <- glw_t + is_edge * (-da_t / pmax(x1 - a, .Machine$double.eps))
  dpsi_m <- ifelse(is_edge, (1 - tfrac) * da_m, 0)
  dpsi_t <- ifelse(is_edge, (1 - tfrac) * da_t, 0)
  # the node-motion part of the prior weight's derivative: d log phi/d psi
  dlphi_dpsi <- ifelse(is.finite(z), -z / tau, 0)
  list(p = p, lw = lw, z = z,
       glw_m = glw_m + dlphi_dpsi * dpsi_m,
       glw_t = glw_t + dlphi_dpsi * dpsi_t,
       dpsi_m = dpsi_m, dpsi_t = dpsi_t)
}

# the shift each node adds to each of the group's observations, and its
# derivatives in the term's own coefficients and in the position
.marg_seg_shift <- function(term, xg, nd, v) {
  hinge <- outer(xg, nd$p, function(x, p) pmax(x - ifelse(is.finite(p), p,
                                                          Inf), 0))
  step <- outer(xg, nd$p, function(x, p) as.numeric(is.finite(p) & x >= p))
  shift <- 0
  if (term@linear) shift <- shift + v[["beta"]] * xg
  shift <- shift + v[["gamma1"]] * hinge
  if (term@kind == "jseg") shift <- shift + v[["delta1"]] * step
  # d shift/d psi = -gamma 1(x > psi); the step's derivative is zero away
  # from the point
  dshift_dpsi <- -v[["gamma1"]] * (hinge > 0)
  list(shift = shift, hinge = hinge, step = step, dshift_dpsi = dshift_dpsi)
}

.marg_seg_loglik <- function(term, eta, y, logdens, score, psi) {
  bp <- .marg_built(term)
  v <- .marg_check_psi(term, psi)
  nm <- term_params(term)
  n <- bp$n
  ll <- numeric(n)
  jac <- matrix(0, n, length(nm), dimnames = list(NULL, nm))
  for (rs in bp$groups) {
    ng <- length(rs)
    nd <- .marg_nodes(bp$x[rs], v[["m1"]], v[["tau1"]])
    C <- length(nd$p)
    sh <- .marg_seg_shift(term, bp$x[rs], nd, v)
    ei <- rep(eta[rs], C) + as.numeric(sh$shift)
    ii <- rep(rs, C)
    LD <- matrix(as.numeric(logdens(ei, ii)), ng, C)
    SC <- matrix(as.numeric(score(ei, ii)), ng, C)

    A <- nd$lw
    tot <- .marg_lse(A)
    u <- exp(A - tot)
    acc <- list()
    own <- c(if (term@linear) "beta", "gamma1",
             if (term@kind == "jseg") "delta1")
    dmat <- list(beta = if (term@linear) matrix(bp$x[rs], ng, C) else NULL,
                 gamma1 = sh$hinge,
                 delta1 = if (term@kind == "jseg") sh$step else NULL)
    for (p in own) acc[[p]] <- numeric(C)
    accp <- numeric(C)   # the psi-motion accumulator, sum_t s * dshift/dpsi
    for (t in seq_len(ng)) {
      row <- rs[t]
      A2 <- A + LD[t, ]
      tot2 <- .marg_lse(A2)
      w <- exp(A2 - tot2)
      ll[row] <- tot2 - tot
      for (p in own) {
        a2 <- acc[[p]] + SC[t, ] * dmat[[p]][t, ]
        jac[row, p] <- sum(w * a2) - sum(u * acc[[p]])
        acc[[p]] <- a2
      }
      ap2 <- accp + SC[t, ] * sh$dshift_dpsi[t, ]
      jac[row, "m1"] <- sum(w * (nd$glw_m + ap2 * nd$dpsi_m)) -
        sum(u * (nd$glw_m + accp * nd$dpsi_m))
      jac[row, "tau1"] <- sum(w * (nd$glw_t + ap2 * nd$dpsi_t)) -
        sum(u * (nd$glw_t + accp * nd$dpsi_t))
      accp <- ap2
      A <- A2
      tot <- tot2
      u <- w
    }
  }
  list(loglik = ll, jacobian = jac)
}

#' @title Posterior Components of a Marginal Break-Point Term
#' @name term_posterior.MarginalBreakTerm
#' @description
#' The component weights Fisher's identity takes at every observation: for
#' the step kind, the posterior probability of each side pattern given the
#' whole sample, one column per pattern of active break-points; for the
#' continuous kinds, each group's posterior over its quadrature nodes,
#' repeated down the group's rows and zero-padded to the widest group.
#' @param term A built \code{\link{MarginalBreakTerm}}.
#' @param eta The static predictor.
#' @param y The response.
#' @param logdens The log-density.
#' @param psi The term's parameters.
#' @param ... Unused.
#' @return A numeric matrix, one row per observation, rows summing to one.
#' @keywords internal
S7::method(term_posterior, MarginalBreakTerm) <- function(term, eta, y,
                                                          logdens, psi, ...) {
  if (term@kind == "jump") {
    .marg_jump_posterior(term, eta, y, logdens, psi)
  } else {
    .marg_seg_posterior(term, eta, y, logdens, psi)$gamma
  }
}

# cumulative sums of a cell array along every coordinate, for box sums by
# inclusion-exclusion
.marg_prefix <- function(A) {
  d <- dim(A)
  if (is.null(d)) return(cumsum(A))
  K <- length(d)
  for (k in seq_len(K)) {
    A <- apply(A, setdiff(seq_len(K), k), cumsum)
    # apply puts the cumulated coordinate first; rotate it back into place
    A <- aperm(A, order(c(k, setdiff(seq_len(K), k))))
  }
  A
}

# the sum of a cell array over a box, from its prefix sums
.marg_boxsum <- function(S, lo, hi) {
  K <- length(lo)
  tot <- 0
  for (c0 in seq_len(2^K) - 1L) {
    ix <- integer(K)
    sgn <- 1
    ok <- TRUE
    for (k in seq_len(K)) {
      take_lo <- bitwAnd(c0, 2^(k - 1L)) > 0
      if (take_lo) {
        ix[k] <- lo[k] - 1L
        sgn <- -sgn
        if (ix[k] < 1L) {
          ok <- FALSE
          break
        }
      } else {
        ix[k] <- hi[k]
      }
    }
    if (!ok) next
    tot <- tot + sgn * S[matrix(ix, 1L)]
  }
  tot
}

# the final cell posterior of one group of the step kind
.marg_jump_cells <- function(term, xg, rows, v, LD) {
  ng <- length(rows)
  J <- ng + 1L
  K <- term@npsi
  pr <- .marg_jump_prior(term, xg, v)
  A <- .marg_outer_sum(pr$lm)
  for (t in seq_len(ng)) {
    lf <- LD[rows[t], ][.marg_pattern(t, J, K) + 1]
    if (K > 1L) dim(lf) <- dim(A)
    A <- A + lf
  }
  w <- exp(A - .marg_lse(A))
  list(w = w, pr = pr, J = J)
}

.marg_jump_posterior <- function(term, eta, y, logdens, psi) {
  bp <- .marg_built(term)
  v <- .marg_check_psi(term, psi)
  K <- term@npsi
  P <- 2^K
  n <- bp$n
  idx <- seq_len(n)
  pb <- .marg_bits(K, v[paste0("delta", seq_len(K))])
  LD <- matrix(0, n, P)
  for (p in seq_len(P)) {
    LD[, p] <- as.numeric(logdens(eta + pb$shifts[p], idx))
  }
  out <- matrix(0, n, P)
  for (rs in bp$groups) {
    ng <- length(rs)
    cl <- .marg_jump_cells(term, bp$x[rs], rs, v, LD)
    if (K == 1L) {
      cw <- cumsum(cl$w)
      for (t in seq_len(ng)) {
        out[rs[t], 2L] <- cw[t]
        out[rs[t], 1L] <- 1 - cw[t]
      }
    } else {
      S <- .marg_prefix(cl$w)
      J <- cl$J
      for (t in seq_len(ng)) {
        for (p in seq_len(P)) {
          lo <- ifelse(pb$bits[p, ], 1L, t + 1L)
          hi <- ifelse(pb$bits[p, ], t, J)
          out[rs[t], p] <- if (any(lo > hi)) 0 else
            .marg_boxsum(S, as.integer(lo), as.integer(hi))
        }
      }
    }
  }
  # rounding in the box sums; the rows are probabilities by construction
  out[out < 0] <- 0
  sw <- rowSums(out)
  out[sw > 0, ] <- out[sw > 0, , drop = FALSE] / sw[sw > 0]
  out
}

# the node posterior of the continuous kinds, with the shift matrix the
# levels method reads: both padded to the widest group
.marg_seg_posterior <- function(term, eta, y, logdens, psi) {
  bp <- .marg_built(term)
  v <- .marg_check_psi(term, psi)
  n <- bp$n
  states <- list()
  Cmax <- 0L
  for (g in seq_along(bp$groups)) {
    rs <- bp$groups[[g]]
    nd <- .marg_nodes(bp$x[rs], v[["m1"]], v[["tau1"]])
    sh <- .marg_seg_shift(term, bp$x[rs], nd, v)
    C <- length(nd$p)
    ei <- rep(eta[rs], C) + as.numeric(sh$shift)
    LD <- matrix(as.numeric(logdens(ei, rep(rs, C))), length(rs), C)
    la <- nd$lw + colSums(LD)
    w <- exp(la - .marg_lse(la))
    states[[g]] <- list(rs = rs, w = w, sh = sh, nd = nd)
    Cmax <- max(Cmax, C)
  }
  gamma <- matrix(0, n, Cmax)
  shift <- matrix(0, n, Cmax)
  for (st in states) {
    C <- length(st$w)
    gamma[st$rs, seq_len(C)] <- matrix(st$w, length(st$rs), C, byrow = TRUE)
    shift[st$rs, seq_len(C)] <- st$sh$shift
  }
  list(gamma = gamma, shift = shift, states = states)
}

#' @title Levels of a Marginal Break-Point Term
#' @name term_levels.MarginalBreakTerm
#' @description
#' For the step kind, the constant shift of each side pattern, the sums of
#' the changes of level over the active break-points. For the continuous
#' kinds the shift varies by observation -- each node's hinge value -- and
#' a matrix is returned, aligned with \code{\link{term_posterior}}'s
#' columns; it takes the callbacks because the node set is theirs to
#' rebuild.
#' @param term A built \code{\link{MarginalBreakTerm}}.
#' @param psi The term's parameters.
#' @param eta,y,logdens For the continuous kinds, the quantities the node
#'   set is built from; ignored by the step kind.
#' @param ... Unused.
#' @return A numeric vector, or a matrix with one row per observation.
#' @keywords internal
S7::method(term_levels, MarginalBreakTerm) <- function(term, psi, eta = NULL,
                                                       y = NULL,
                                                       logdens = NULL, ...) {
  v <- .marg_check_psi(term, psi)
  if (term@kind == "jump") {
    return(.marg_bits(term@npsi,
                      v[paste0("delta", seq_len(term@npsi))])$shifts)
  }
  if (is.null(eta) || is.null(logdens)) {
    stop(paste("the continuous kinds shift each observation by its node's",
               "hinge value: pass eta, y and logdens so the node set can be",
               "rebuilt."), call. = FALSE)
  }
  .marg_seg_posterior(term, eta, y, logdens, psi)$shift
}

#' @title Posterior Break-Points of a Marginal Term
#' @name term_latent.MarginalBreakTerm
#' @description
#' The posterior mean and standard deviation of each group's break-points.
#' For the step kind under the gaussian prior these are mixtures of
#' truncated-normal moments over the interval posterior, the edge
#' intervals read through \code{\link[numericals7]{mills_ratio}}; under an
#' explicit prior the truncated moments come from
#' \code{\link[distributions7]{truncated}} and
#' \code{\link[distributions7]{expectation}}, one interval at a time. For
#' the continuous kinds they are the moments of the node posterior, the
#' closed upper tail entering through its truncated-normal moments.
#' @param term A built \code{\link{MarginalBreakTerm}}.
#' @param eta The static predictor.
#' @param y The response.
#' @param logdens The log-density.
#' @param psi The term's parameters.
#' @param ... Unused.
#' @return A data frame with \code{group}, \code{psi}, \code{mean} and
#'   \code{sd}.
#' @keywords internal
S7::method(term_latent, MarginalBreakTerm) <- function(term, eta, y, logdens,
                                                       psi, ...) {
  bp <- .marg_built(term)
  v <- .marg_check_psi(term, psi)
  K <- term@npsi
  n <- bp$n
  idx <- seq_len(n)
  rows <- list()
  if (term@kind == "jump") {
    P <- 2^K
    pb <- .marg_bits(K, v[paste0("delta", seq_len(K))])
    LD <- matrix(0, n, P)
    for (p in seq_len(P)) {
      LD[, p] <- as.numeric(logdens(eta + pb$shifts[p], idx))
    }
    for (g in seq_along(bp$groups)) {
      rs <- bp$groups[[g]]
      cl <- .marg_jump_cells(term, bp$x[rs], rs, v, LD)
      for (k in seq_len(K)) {
        wm <- .marg_margin(cl$w, k)
        mo <- .marg_interval_moments(term, bp$x[rs], v, k, wm, cl$pr)
        rows[[length(rows) + 1L]] <- data.frame(
          group = bp$labels[g], psi = k, mean = mo$mean, sd = mo$sd,
          stringsAsFactors = FALSE)
      }
    }
  } else {
    ps <- .marg_seg_posterior(term, eta, y, logdens, psi)
    for (g in seq_along(ps$states)) {
      st <- ps$states[[g]]
      w <- st$w
      p <- st$nd$p
      C <- length(p)
      fin <- is.finite(p)
      m1 <- sum(w[fin] * p[fin])
      m2 <- sum(w[fin] * p[fin]^2)
      # the closed upper tail: the truncated normal above the last
      # observation
      wt <- w[C]
      if (wt > 0) {
        xn <- max(bp$x[st$rs])
        zn <- (xn - v[["m1"]]) / v[["tau1"]]
        lam <- numericals7::mills_ratio(-zn)$r
        mu_t <- v[["m1"]] + v[["tau1"]] * lam
        var_t <- v[["tau1"]]^2 * (1 + zn * lam - lam^2)
        m1 <- m1 + wt * mu_t
        m2 <- m2 + wt * (var_t + mu_t^2)
      }
      rows[[length(rows) + 1L]] <- data.frame(
        group = bp$labels[g], psi = 1L, mean = m1,
        sd = sqrt(max(m2 - m1^2, 0)), stringsAsFactors = FALSE)
    }
  }
  do.call(rbind, rows)
}

# the truncated moments of one break-point's prior over its intervals,
# mixed by the posterior marginal
.marg_interval_moments <- function(term, xs, v, k, wm, pr) {
  J <- length(wm)
  if (is.null(term@prior)) {
    m <- v[[paste0("m", k)]]
    tau <- v[[paste0("tau", k)]]
    iv <- pr$ivs[[k]]
    r1 <- ifelse(iv$ok, (iv$phi_l - iv$phi_u) / iv$mass, 0)
    r2 <- ifelse(iv$ok, (iv$zp_l - iv$zp_u) / iv$mass, 0)
    mr_u <- numericals7::mills_ratio(iv$zu[1L])$r
    mr_l <- numericals7::mills_ratio(-iv$zl[J])$r
    r1[1L] <- -mr_u
    r2[1L] <- -iv$zu[1L] * mr_u
    r1[J] <- mr_l
    r2[J] <- iv$zl[J] * mr_l
    mean_j <- m + tau * r1
    var_j <- tau^2 * (1 + r2 - r1^2)
    keep <- wm > 0
    mu <- sum(wm[keep] * mean_j[keep])
    m2 <- sum(wm[keep] * (var_j[keep] + mean_j[keep]^2))
    return(list(mean = mu, sd = sqrt(max(m2 - mu^2, 0))))
  }
  # an explicit prior: the moments of the prior truncated to each interval
  # that carries appreciable weight, through the batched engines. A moment
  # the engine cannot deliver is reported NA rather than approximated: a
  # heavy-tailed prior's edge intervals keep no mean below one degree of
  # freedom and no variance below two, and the quadrature's refusal is the
  # honest reading of that.
  iv <- pr$ivs[[k]]
  th <- stats::setNames(as.list(v[term@prior@params]), term@prior@params)
  mu <- 0
  m2 <- 0
  for (j in seq_len(J)) {
    if (wm[j] < 1e-9) next
    # an infinite endpoint is omitted rather than passed: truncating at the
    # support's own bound removes no mass and truncated() rejects it
    tr_args <- list(term@prior)
    if (is.finite(iv$l[j])) tr_args$lower <- iv$l[j] - v[["m1"]]
    if (is.finite(iv$u[j])) tr_args$upper <- iv$u[j] - v[["m1"]]
    tr <- do.call(distributions7::truncated, tr_args)
    e1 <- tryCatch(as.numeric(distributions7::expectation(
      tr, function(y, theta) y, th)), error = function(e) NA_real_)
    e2 <- tryCatch(as.numeric(distributions7::expectation(
      tr, function(y, theta) y^2, th)), error = function(e) NA_real_)
    mu <- mu + wm[j] * (v[["m1"]] + e1)
    m2 <- m2 + wm[j] * (e2 + 2 * v[["m1"]] * e1 + v[["m1"]]^2)
  }
  sc <- sum(wm[wm >= 1e-9])
  list(mean = mu / sc,
       sd = if (is.finite(m2) && is.finite(mu)) {
         sqrt(max(m2 / sc - (mu / sc)^2, 0))
       } else NA_real_)
}

#' @title Where a Marginal Break-Point Term Starts
#' @name term_start.MarginalBreakTerm
#' @description
#' Given a target (the response on the scale of the predictor, which the
#' fitting layer supplies), the exact profile: each group's best split by
#' least squares on its own observations under the term's own design, the
#' population positions and the prior scales read off those, and the
#' changes their pooled means. Without one, the covariate's quantiles and
#' a scale from its spread, with the changes at zero.
#' @param term A built \code{\link{MarginalBreakTerm}}.
#' @param ... Unused.
#' @param target The response on the predictor's scale, or \code{NULL}.
#' @return A named numeric vector on the unconstrained scale.
#' @keywords internal
S7::method(term_start, MarginalBreakTerm) <- function(term, ...,
                                                      target = NULL) {
  nm <- term_params(term)
  links <- term_links(term)
  out <- stats::setNames(numeric(length(nm)), nm)
  bp <- term@blueprint
  if (!length(bp)) return(out)
  xv <- bp$x
  K <- term@npsi
  qs <- stats::quantile(xv, c(0.05, 0.95), names = FALSE)
  m0 <- if (!is.null(term@spec$psi)) sort(term@spec$psi) else {
    as.numeric(stats::quantile(xv, seq_len(K) / (K + 1), names = FALSE))
  }
  tau0 <- rep((qs[2L] - qs[1L]) / 6, K)
  ch <- list(delta = rep(0, K), gamma = rep(0, K), beta = 0)
  if (!is.null(target) && length(target) == bp$n) {
    prof <- .marg_profile_start(term, bp, as.numeric(target))
    if (!is.null(prof)) {
      if (is.null(term@spec$psi)) m0 <- prof$m
      tau0 <- prof$tau
      ch <- prof$ch
    }
  }
  tau0 <- pmax(tau0, .marg_gap_floor(bp))
  for (k in seq_len(K)) {
    if (is.null(term@prior)) {
      out[[paste0("m", k)]] <- m0[k]
      out[[paste0("tau", k)]] <- log(tau0[k])
    }
  }
  if (!is.null(term@prior)) {
    out[["m1"]] <- m0[1L]
    # a prior parameter read as a scale starts at the profile's spread,
    # carried through its own chart; the others keep the chart's zero
    ip <- tryCatch(term@prior@params_interpretation, error = function(e) NULL)
    for (p in term@prior@params) {
      if (!is.null(ip) && !is.null(ip[[p]]) &&
          ip[[p]] %in% c("scale", "standard deviation")) {
        z <- tryCatch(linkfunctions7::linkfun(links[[p]], tau0[1L]),
                      error = function(e) NA_real_)
        if (is.finite(z)) out[[p]] <- z
      }
    }
  }
  if (term@linear) out[["beta"]] <- ch$beta
  for (k in seq_len(K)) {
    if (term@kind %in% c("jump", "jseg")) {
      out[[paste0("delta", k)]] <- ch$delta[k]
    }
    if (term@kind %in% c("seg", "jseg")) {
      out[[paste0("gamma", k)]] <- ch$gamma[k]
    }
  }
  out
}

# Half the median gap between a group's consecutive distinct covariate
# values, pooled over the groups: the finest scale at which the data can
# tell two positions apart, and the floor under a starting prior scale.
.marg_gap_floor <- function(bp) {
  gaps <- unlist(lapply(bp$groups, function(rs) {
    u <- sort(unique(bp$x[rs]))
    if (length(u) > 1L) diff(u) else numeric(0)
  }))
  if (!length(gaps)) gaps <- diff(range(bp$x))
  max(stats::median(gaps) / 2, sqrt(.Machine$double.eps))
}

# The exact least-squares profile of the target under the term's own
# design, in two stages. The POOLED stage estimates the population
# positions and the changes on all the groups at once, per-group
# intercepts absorbing whatever levels the target carries, the positions
# taken greedily off a quantile grid -- one group's own profile is not
# used for the coefficients, because a per-group fit with the full design
# on a noisy target overfits its own noise (measured on a Poisson panel:
# per-group linear effects with median -4.7 against a truth of 0.15, and
# the search then started in the wrong basin). The PER-GROUP stage moves
# one position at a time with everything else held at the pooled
# estimates, which is one free value per group and cannot overfit; the
# spread of those positions is the prior scale's start.
.marg_profile_start <- function(term, bp, target) {
  K <- term@npsi
  n <- bp$n
  xv <- bp$x
  keep <- is.finite(target)
  if (sum(keep) < K + 3L) return(NULL)
  G <- length(bp$groups)
  gvec <- integer(n)
  for (g in seq_len(G)) gvec[bp$groups[[g]]] <- g
  Gm <- matrix(0, n, G)
  Gm[cbind(seq_len(n), gvec)] <- 1

  cols_at <- function(cs) {
    X <- if (term@linear) cbind(Gm, xv) else Gm
    for (c0 in cs) {
      if (term@kind %in% c("seg", "jseg")) X <- cbind(X, pmax(xv - c0, 0))
      if (term@kind %in% c("jump", "jseg")) {
        X <- cbind(X, as.numeric(xv >= c0))
      }
    }
    X
  }
  cand <- unique(as.numeric(stats::quantile(
    xv[keep], seq(0.06, 0.94, length.out = 23L), names = FALSE)))
  pooled_rss <- function(cs) {
    X <- cols_at(cs)
    fit <- tryCatch(stats::lm.fit(X[keep, , drop = FALSE], target[keep]),
                    error = function(e) NULL)
    if (is.null(fit)) return(list(rss = Inf, cf = NULL))
    list(rss = sum(fit$residuals^2), cf = fit$coefficients)
  }

  # the columns of the final scoring fit, with PER-GROUP positions: shared
  # coefficients over each group's own hinge and step
  cols_per_group <- function(psis) {
    X <- if (term@linear) cbind(Gm, xv) else Gm
    pg <- psis[gvec, , drop = FALSE]
    for (k in seq_len(K)) {
      if (term@kind %in% c("seg", "jseg")) {
        X <- cbind(X, pmax(xv - pg[, k], 0))
      }
      if (term@kind %in% c("jump", "jseg")) {
        X <- cbind(X, as.numeric(xv >= pg[, k]))
      }
    }
    X
  }
  read_ch <- function(cf) {
    cf[!is.finite(cf)] <- 0
    pos <- G + as.integer(term@linear)
    beta0 <- if (term@linear) cf[G + 1L] else 0
    gamma0 <- numeric(K)
    delta0 <- numeric(K)
    for (k in seq_len(K)) {
      if (term@kind == "seg") gamma0[k] <- cf[pos + k]
      if (term@kind == "jump") delta0[k] <- cf[pos + k]
      if (term@kind == "jseg") {
        gamma0[k] <- cf[pos + 2L * k - 1L]
        delta0[k] <- cf[pos + 2L * k]
      }
    }
    list(beta = unname(beta0), gamma = gamma0, delta = delta0)
  }
  refine <- function(found, ch) {
    contrib <- function(cs) {
      v <- if (term@linear) ch$beta * xv else numeric(n)
      for (k in seq_along(cs)) {
        if (term@kind %in% c("seg", "jseg")) {
          v <- v + ch$gamma[k] * pmax(xv - cs[k], 0)
        }
        if (term@kind %in% c("jump", "jseg")) {
          v <- v + ch$delta[k] * as.numeric(xv >= cs[k])
        }
      }
      v
    }
    psis <- matrix(NA_real_, G, K)
    for (g in seq_len(G)) {
      rs <- bp$groups[[g]]
      okr <- rs[is.finite(target[rs])]
      if (length(okr) < 3L) next
      xs <- xv[okr]
      dist <- which(xs[-length(xs)] < xs[-1L])
      if (!length(dist)) next
      mids <- (xs[dist] + xs[dist + 1L]) / 2
      for (k in seq_len(K)) {
        best <- Inf
        bc <- found[k]
        for (c0 in mids) {
          cs <- found
          cs[k] <- c0
          r <- target[okr] - contrib(cs)[okr]
          rss <- sum((r - mean(r))^2)
          if (rss < best) {
            best <- rss
            bc <- c0
          }
        }
        psis[g, k] <- bc
      }
      psis[g, ] <- sort(psis[g, ])
    }
    psis
  }

  # The pooled single-position profile is nearly flat when the positions
  # vary by group -- their spread smears the break -- and its global
  # minimum is not reliably in the right basin (measured on a Poisson
  # panel: 103.6 at c = 8.3 against 104.8 at c = 4.4, the truth's basin
  # the runner-up). Every local minimum is therefore carried through the
  # per-group refinement, and the winner is the candidate whose FINAL fit,
  # shared coefficients over per-group positions, leaves the least behind:
  # at per-group positions the smearing is gone and the true basin wins.
  scan <- vapply(cand, function(c0) pooled_rss(c0)$rss, numeric(1))
  nc <- length(scan)
  locmin <- which(scan <= c(Inf, scan[-nc]) & scan <= c(scan[-1L], Inf) &
                    is.finite(scan))
  if (!length(locmin)) return(NULL)
  locmin <- locmin[order(scan[locmin])]
  locmin <- locmin[seq_len(min(4L, length(locmin)))]

  best <- NULL
  for (i0 in locmin) {
    found <- cand[i0]
    ok_seed <- TRUE
    if (K > 1L) {
      for (k in seq.int(2L, K)) {
        rssk <- vapply(setdiff(cand, found), function(c0)
          pooled_rss(c(found, c0))$rss, numeric(1))
        if (!any(is.finite(rssk))) {
          ok_seed <- FALSE
          break
        }
        found <- c(found, setdiff(cand, found)[which.min(rssk)])
      }
    }
    if (!ok_seed) next
    found <- sort(found)
    pf <- pooled_rss(found)
    if (is.null(pf$cf)) next
    ch <- read_ch(pf$cf)
    psis <- refine(found, ch)
    okg <- stats::complete.cases(psis)
    if (!any(okg)) next
    for (g in which(!okg)) psis[g, ] <- found
    Xf <- cols_per_group(psis)
    ff <- tryCatch(stats::lm.fit(Xf[keep, , drop = FALSE], target[keep]),
                   error = function(e) NULL)
    if (is.null(ff)) next
    rssf <- sum(ff$residuals^2)
    if (is.null(best) || rssf < best$rss) {
      best <- list(rss = rssf, psis = psis, ch = read_ch(ff$coefficients),
                   okg = okg)
    }
  }
  if (is.null(best)) return(NULL)
  psis <- best$psis[best$okg, , drop = FALSE]
  tau <- apply(psis, 2L, stats::mad)
  bad <- !is.finite(tau) | tau <= 0
  if (any(bad)) tau[bad] <- apply(psis, 2L, stats::sd, na.rm = TRUE)
  list(m = apply(psis, 2L, stats::median),
       tau = pmax(ifelse(is.finite(tau), tau, 0), 0),
       ch = best$ch)
}

# ---------------------------------------------------------------------------
# The observed Hessian. The step kind with one gaussian break-point carries
# the fully propagated one, the interval sum differentiated twice; every
# other configuration differences the analytic full gradient once, a single
# central stencil on the analytic order below, which is the licence the
# toolkit's non-closed derivatives already run on. The one-break-point
# gaussian route doubles as the control: the two must agree, and a test
# holds them to it.

# the analytic gradient of the marginal log-likelihood in ALL of a caller's
# unknowns, on the unconstrained scale. Only the level equation's static
# predictor is moved by the term, so it is the one handed over; the other
# equations enter through `grad`, read at the family's current parameters.
.marg_grad_full <- function(term, eta, y, logdens, grad, zeta, seed, cols,
                            level, w) {
  bp <- term@blueprint
  nm <- term_params(term)
  links <- term_links(term)
  v <- stats::setNames(vapply(nm, function(j)
    linkfunctions7::linkinv(links[[j]], zeta[[j]]), numeric(1)), nm)
  chain <- vapply(nm, function(j)
    linkfunctions7::dlinkinv(links[[j]], zeta[[j]]), numeric(1))
  psi <- as.list(v)
  n <- bp$n
  npar <- length(seed)
  mm <- ncol(seed[[1L]])
  gout <- numeric(mm)
  own_g <- stats::setNames(numeric(length(nm)), nm)
  idx <- seq_len(n)

  if (term@kind == "jump") {
    K <- term@npsi
    P <- 2^K
    pb <- .marg_bits(K, v[paste0("delta", seq_len(K))])
    LD <- matrix(0, n, P)
    G <- vector("list", P)
    for (p in seq_len(P)) {
      LD[, p] <- as.numeric(logdens(eta + pb$shifts[p], idx))
      G[[p]] <- as.matrix(grad(eta + pb$shifts[p], idx))
    }
    gamma <- .marg_jump_posterior(term, eta, y, logdens, psi)
    # the coefficient block by Fisher's identity, every predictor weighted
    # by the pattern posterior
    gv <- matrix(0, n, npar)
    for (p in seq_len(P)) gv <- gv + gamma[, p] * G[[p]]
    for (q in seq_len(npar)) {
      gout <- gout + as.numeric(crossprod(seed[[q]], w * gv[, q]))
    }
    # the changes of level: the level equation's score where the pattern
    # activates the break-point
    for (k in seq_len(K)) {
      s <- 0
      for (p in seq_len(P)) {
        if (pb$bits[p, k]) s <- s + sum(w * gamma[, p] * G[[p]][, level])
      }
      own_g[[paste0("delta", k)]] <- own_g[[paste0("delta", k)]] + s
    }
    # the prior's parameters: interval-posterior expectations of the mass
    # derivatives, per group
    for (g in seq_along(bp$groups)) {
      rs <- bp$groups[[g]]
      wi <- w[rs][1L]
      cl <- .marg_jump_cells(term, bp$x[rs], rs, v, LD)
      for (k in seq_along(cl$pr$dlm)) {
        wm <- .marg_margin(cl$w, k)
        contrib <- as.numeric(wm %*% cl$pr$dlm[[k]])
        pc <- colnames(cl$pr$dlm[[k]])
        own_g[pc] <- own_g[pc] + wi * contrib
      }
    }
  } else {
    ps <- .marg_seg_posterior(term, eta, y, logdens, psi)
    for (st in ps$states) {
      rs <- st$rs
      wi <- w[rs][1L]
      C <- length(st$w)
      ng <- length(rs)
      ei <- rep(eta[rs], C) + as.numeric(st$sh$shift)
      ii <- rep(rs, C)
      Gf <- as.matrix(grad(ei, ii))
      # the coefficient block: per-observation node-weighted scores
      for (q in seq_len(npar)) {
        Gq <- matrix(Gf[, q], ng, C)
        gvq <- as.numeric(Gq %*% st$w)
        gout <- gout +
          wi * as.numeric(crossprod(seed[[q]][rs, , drop = FALSE], gvq))
      }
      GL <- matrix(Gf[, level], ng, C)
      own <- list(beta = if (term@linear) matrix(bp$x[rs], ng, C),
                  gamma1 = st$sh$hinge,
                  delta1 = if (term@kind == "jseg") st$sh$step)
      for (p in names(own)) {
        if (is.null(own[[p]])) next
        own_g[[p]] <- own_g[[p]] + wi * sum(st$w * colSums(GL * own[[p]]))
      }
      cpsi <- colSums(GL * st$sh$dshift_dpsi)
      own_g[["m1"]] <- own_g[["m1"]] +
        wi * sum(st$w * (st$nd$glw_m + cpsi * st$nd$dpsi_m))
      own_g[["tau1"]] <- own_g[["tau1"]] +
        wi * sum(st$w * (st$nd$glw_t + cpsi * st$nd$dpsi_t))
    }
  }
  # the term's own columns onto the unconstrained scale
  gout[cols] <- gout[cols] + own_g * chain
  gout
}

#' @title Observed Hessian of a Marginal Break-Point Term
#' @name term_hessian.MarginalBreakTerm
#' @description
#' The exact observed Hessian of the marginal log-likelihood over a
#' caller's unknowns, the coefficients of every equation together with the
#' term's own parameters on the unconstrained scale.
#' @details
#' The step kind with one gaussian break-point differentiates the interval
#' sum twice, the second derivatives of the conditional collapsing into
#' per-observation Hessians weighted by each side's posterior probability
#' and the mass curvature closed in the normal density. Every other
#' configuration differences the analytic full gradient once -- a single
#' central stencil on the analytic order below, the licence the toolkit's
#' non-closed derivatives run on -- and the one-break-point gaussian route
#' is the control the tests hold it to.
#'
#' The marginal likelihood of a group does not factorize over its
#' observations, so an observation weight has a reading only when it is
#' constant within each group; anything else is rejected.
#' @param term A built \code{\link{MarginalBreakTerm}}.
#' @param eta The static predictor of the level equation.
#' @param y The response.
#' @param logdens,grad,hess The log-density and its first two derivatives
#'   in the predictors.
#' @param psi The term's parameters.
#' @param seed The derivative of each predictor in the caller's unknowns.
#' @param cols The columns the term's own parameters occupy.
#' @param level The distribution parameter the term shifts.
#' @param weights Observation weights, constant within each group.
#' @param ... Unused.
#' @return A list with \code{loglik}, \code{gradient} and \code{hessian}.
#' @keywords internal
S7::method(term_hessian, MarginalBreakTerm) <- function(term, eta, y, logdens,
                                                        grad, hess, psi, seed,
                                                        cols, level,
                                                        weights = NULL, ...) {
  bp <- .marg_built(term)
  v <- .marg_check_psi(term, psi)
  nm <- term_params(term)
  n <- bp$n
  if (!is.list(seed) || !length(seed)) {
    stop("'seed' must be a list with one matrix per distribution parameter.",
         call. = FALSE)
  }
  seed <- lapply(seed, as.matrix)
  mm <- ncol(seed[[1L]])
  npar <- length(seed)
  if (any(vapply(seed, nrow, 1L) != n) ||
      any(vapply(seed, ncol, 1L) != mm)) {
    stop(sprintf("every element of 'seed' must be %d by %d.", n, mm),
         call. = FALSE)
  }
  cols <- as.integer(cols)
  level <- as.integer(level)
  if (length(cols) != length(nm) || any(cols < 1L) || any(cols > mm)) {
    stop(sprintf("'cols' must give the %d columns the term's parameters",
                 length(nm)), call. = FALSE)
  }
  if (length(level) != 1L || level < 1L || level > npar) {
    stop("'level' must index one of the distribution parameters.",
         call. = FALSE)
  }
  w <- if (is.null(weights)) rep(1, n) else rep_len(as.numeric(weights), n)
  for (rs in bp$groups) {
    if (length(unique(w[rs])) > 1L) {
      stop(paste("a marginal likelihood does not factorize over a group's",
                 "observations, so an observation weight has a reading only",
                 "when it is constant within each group."), call. = FALSE)
    }
  }

  if (term@kind == "jump" && term@npsi == 1L && is.null(term@prior)) {
    return(.marg_jump1_hessian(term, eta, y, logdens, grad, hess, psi, seed,
                               cols, level, w))
  }
  if (term@kind == "jump") {
    .marg_jump_hessian(term, eta, y, logdens, grad, hess, psi, seed, cols,
                       level, w)
  } else {
    .marg_seg_hessian(term, eta, y, logdens, grad, hess, psi, seed, cols,
                      level, w)
  }
}

# One central stencil per named column on the analytic full gradient: the
# licence the toolkit's non-closed derivatives run on, one difference on the
# analytic order below. Used for the rows the closed mass curvature does not
# cover -- an explicit prior's parameters, and the node motion of the
# quadrature constructions -- and only along the term's own columns, whose
# perturbation the callbacks support in full.
.marg_fd_rows <- function(H, term, eta, y, logdens, grad, zeta, seed, cols,
                          level, w, which_own) {
  nm <- term_params(term)
  mm <- ncol(seed[[1L]])
  u0 <- numeric(mm)
  u0[cols] <- zeta
  g_at <- function(u) {
    z <- zeta
    z[] <- u[cols]
    .marg_grad_full(term, eta, y, logdens, grad, as.list(z), seed, cols,
                    level, w)
  }
  for (p in which_own) {
    d <- cols[match(p, nm)]
    h <- numericals7::fd_step(u0[d], 1L)
    up <- u0
    um <- u0
    up[d] <- up[d] + h
    um[d] <- um[d] - h
    col <- (g_at(up) - g_at(um)) / (2 * h)
    H[, d] <- col
    H[d, ] <- col
  }
  (H + t(H)) / 2
}

# The exact Hessian of the step kind with several break-points, or one
# break-point under an explicit prior: the cell sum differentiated twice.
# The conditional's second derivatives collapse into per-observation
# Hessians weighted by the pattern posterior; the first-derivative products
# run over the cells, whose per-unknown gradients are accumulated one
# observation at a time; the gaussian masses' own curvature is closed and
# separable across the coordinates, and an explicit prior's rows come from
# one stencil on the analytic gradient.
.marg_jump_hessian <- function(term, eta, y, logdens, grad, hess, psi, seed,
                               cols, level, w) {
  bp <- term@blueprint
  v <- .marg_check_psi(term, psi)
  nm <- term_params(term)
  links <- term_links(term)
  K <- term@npsi
  P <- 2^K
  n <- bp$n
  npar <- length(seed)
  mm <- ncol(seed[[1L]])
  idx <- seq_len(n)
  pb <- .marg_bits(K, v[paste0("delta", seq_len(K))])
  LD <- matrix(0, n, P)
  G <- vector("list", P)
  H2 <- vector("list", P)
  for (p in seq_len(P)) {
    LD[, p] <- as.numeric(logdens(eta + pb$shifts[p], idx))
    G[[p]] <- as.matrix(grad(eta + pb$shifts[p], idx))
    H2[[p]] <- array(as.numeric(hess(eta + pb$shifts[p], idx)),
                     c(n, npar, npar))
  }
  gamma <- .marg_jump_posterior(term, eta, y, logdens, as.list(v))
  chain <- vapply(nm, function(j)
    linkfunctions7::dlinkinv(links[[j]],
                             linkfunctions7::linkfun(links[[j]], v[[j]])),
    numeric(1))

  # per-pattern per-observation scores in the unknowns, the pattern's
  # active break-points adding the level score to their delta columns
  gm <- vector("list", P)
  for (p in seq_len(P)) {
    g0 <- matrix(0, n, mm)
    for (q in seq_len(npar)) g0 <- g0 + G[[p]][, q] * seed[[q]]
    for (k in seq_len(K)) {
      if (pb$bits[p, k]) {
        dc <- cols[match(paste0("delta", k), nm)]
        g0[, dc] <- g0[, dc] + G[[p]][, level]
      }
    }
    gm[[p]] <- g0
  }

  loglik <- numeric(n)
  gradient <- numeric(mm)
  hessian <- matrix(0, mm, mm)

  for (g in seq_along(bp$groups)) {
    rs <- bp$groups[[g]]
    ng <- length(rs)
    J <- ng + 1L
    wi <- w[rs][1L]
    pr <- .marg_jump_prior(term, bp$x[rs], v)
    A <- .marg_outer_sum(pr$lm)
    ncell <- length(A)
    # the per-cell gradient in the unknowns, one observation at a time: a
    # cell's pattern for one observation is a box, so its rows share the
    # observation's per-pattern score
    DA <- matrix(0, ncell, mm)
    tot0 <- .marg_lse(A)
    for (t in seq_len(ng)) {
      row <- rs[t]
      PA <- .marg_pattern(t, J, K)
      lf <- LD[row, ][PA + 1]
      if (K > 1L) dim(lf) <- dim(A)
      A2 <- A + lf
      tot2 <- .marg_lse(A2)
      loglik[row] <- tot2 - tot0
      tot0 <- tot2
      A <- A2
      byp <- split(seq_len(ncell), as.integer(PA))
      for (pp in names(byp)) {
        p1 <- as.integer(pp) + 1L
        DA[byp[[pp]], ] <- DA[byp[[pp]], , drop = FALSE] +
          matrix(gm[[p1]][row, ], length(byp[[pp]]), mm, byrow = TRUE)
      }
    }
    wjA <- exp(A - .marg_lse(A))
    wj <- as.numeric(wjA)
    # the masses' first derivatives, on the unconstrained scale, and the
    # closed curvature of the gaussian ones
    for (k in seq_along(pr$dlm)) {
      jk <- as.integer(.marg_outer_sum(
        c(rep(list(numeric(J)), k - 1L), list(seq_len(J)),
          rep(list(numeric(J)), length(pr$lm) - k))))
      pc <- colnames(pr$dlm[[k]])
      dc <- cols[match(pc, nm)]
      DA[, dc] <- DA[, dc] +
        pr$dlm[[k]][jk, , drop = FALSE] %*% diag(chain[match(pc, nm)],
                                                 length(pc))
    }
    gbar <- as.numeric(crossprod(DA, wj))
    gradient <- gradient + wi * gbar
    hessian <- hessian + wi * (crossprod(DA, wj * DA) - outer(gbar, gbar))
    if (is.null(term@prior)) {
      for (k in seq_len(K)) {
        iv <- .marg_intervals(bp$x[rs], v[[paste0("m", k)]],
                              v[[paste0("tau", k)]], d2 = TRUE)
        tau <- v[[paste0("tau", k)]]
        wm <- .marg_margin(wjA, k)
        i1 <- cols[match(paste0("m", k), nm)]
        i2 <- cols[match(paste0("tau", k), nm)]
        hessian[i1, i1] <- hessian[i1, i1] + wi * sum(wm * iv$dmm)
        smt <- sum(wm * iv$dmt) * tau
        hessian[i1, i2] <- hessian[i1, i2] + wi * smt
        hessian[i2, i1] <- hessian[i2, i1] + wi * smt
        hessian[i2, i2] <- hessian[i2, i2] +
          wi * (sum(wm * iv$dtt) * tau^2 + sum(wm * iv$dt) * tau)
      }
    }
  }

  # the conditional's second derivatives, collapsed over the cells: each
  # observation's per-pattern Hessian enters weighted by its pattern
  # posterior
  for (p in seq_len(P)) {
    seedp <- seed
    for (k in seq_len(K)) {
      if (pb$bits[p, k]) {
        dc <- cols[match(paste0("delta", k), nm)]
        seedp[[level]][, dc] <- seedp[[level]][, dc] + 1
      }
    }
    ws <- w * gamma[, p]
    for (aq in seq_len(npar)) {
      for (bq in seq_len(npar)) {
        hessian <- hessian +
          crossprod(seedp[[aq]], (ws * H2[[p]][, aq, bq]) * seedp[[bq]])
      }
    }
  }
  hessian <- (hessian + t(hessian)) / 2

  if (!is.null(term@prior)) {
    zeta <- vapply(nm, function(j)
      linkfunctions7::linkfun(links[[j]], v[[j]]), numeric(1))
    hessian <- .marg_fd_rows(hessian, term, eta, y, logdens, grad, zeta,
                             seed, cols, level, w,
                             c("m1", term@prior@params))
  }
  list(loglik = loglik, gradient = gradient, hessian = hessian)
}

# The Hessian of the continuous kinds: the node sum differentiated twice at
# fixed nodes -- the conditional's second derivatives collapsing into
# per-node per-observation Hessians, the first-derivative products over the
# node posterior -- with the prior's rows, where the panels below the data
# move with its parameters, from one stencil on the analytic gradient.
.marg_seg_hessian <- function(term, eta, y, logdens, grad, hess, psi, seed,
                              cols, level, w) {
  bp <- term@blueprint
  v <- .marg_check_psi(term, psi)
  nm <- term_params(term)
  links <- term_links(term)
  n <- bp$n
  npar <- length(seed)
  mm <- ncol(seed[[1L]])
  chain <- vapply(nm, function(j)
    linkfunctions7::dlinkinv(links[[j]],
                             linkfunctions7::linkfun(links[[j]], v[[j]])),
    numeric(1))
  ps <- .marg_seg_posterior(term, eta, y, logdens, as.list(v))
  loglik <- numeric(n)
  gradient <- numeric(mm)
  hessian <- matrix(0, mm, mm)

  for (st in ps$states) {
    rs <- st$rs
    ng <- length(rs)
    wi <- w[rs][1L]
    C <- length(st$w)
    ei <- rep(eta[rs], C) + as.numeric(st$sh$shift)
    ii <- rep(rs, C)
    LD <- matrix(as.numeric(logdens(ei, ii)), ng, C)
    Gf <- as.matrix(grad(ei, ii))
    Hf <- array(as.numeric(hess(ei, ii)), c(ng * C, npar, npar))
    GL <- matrix(Gf[, level], ng, C)

    # the contributions, from the same states
    A <- st$nd$lw
    tot <- .marg_lse(A)
    for (t in seq_len(ng)) {
      A <- A + LD[t, ]
      tot2 <- .marg_lse(A)
      loglik[rs[t]] <- tot2 - tot
      tot <- tot2
    }

    # the per-node partials of the shift in the term's own coefficients,
    # and the level column each node's shift augments
    ownm <- list(beta = if (term@linear) matrix(bp$x[rs], ng, C),
                 gamma1 = st$sh$hinge,
                 delta1 = if (term@kind == "jseg") st$sh$step)
    ownm <- ownm[!vapply(ownm, is.null, logical(1))]

    # the per-node gradient in the unknowns
    DA <- matrix(0, C, mm)
    for (q in seq_len(npar)) {
      Gq <- matrix(Gf[, q], ng, C)
      DA <- DA + crossprod(Gq, seed[[q]][rs, , drop = FALSE])
    }
    for (p in names(ownm)) {
      dc <- cols[match(p, nm)]
      DA[, dc] <- DA[, dc] + colSums(GL * ownm[[p]]) * chain[match(p, nm)]
    }
    cpsi <- colSums(GL * st$sh$dshift_dpsi)
    im <- cols[match("m1", nm)]
    it <- cols[match("tau1", nm)]
    DA[, im] <- DA[, im] +
      (st$nd$glw_m + cpsi * st$nd$dpsi_m) * chain[match("m1", nm)]
    DA[, it] <- DA[, it] +
      (st$nd$glw_t + cpsi * st$nd$dpsi_t) * chain[match("tau1", nm)]
    gbar <- as.numeric(crossprod(DA, st$w))
    gradient <- gradient + wi * gbar
    hessian <- hessian + wi * (crossprod(DA, st$w * DA) - outer(gbar, gbar))

    # the conditional's second derivatives at each node, the level row of
    # the per-observation Jacobian augmented by the node's own partials
    for (c0 in seq_len(C)) {
      wc <- wi * st$w[c0]
      if (wc == 0) next
      rows0 <- (c0 - 1L) * ng + seq_len(ng)
      seedc <- lapply(seed, function(s) s[rs, , drop = FALSE])
      aug <- matrix(0, ng, mm)
      for (p in names(ownm)) {
        aug[, cols[match(p, nm)]] <- ownm[[p]][, c0] * chain[match(p, nm)]
      }
      aug[, im] <- aug[, im] +
        st$sh$dshift_dpsi[, c0] * st$nd$dpsi_m[c0] * chain[match("m1", nm)]
      aug[, it] <- aug[, it] +
        st$sh$dshift_dpsi[, c0] * st$nd$dpsi_t[c0] * chain[match("tau1", nm)]
      seedc[[level]] <- seedc[[level]] + aug
      for (aq in seq_len(npar)) {
        for (bq in seq_len(npar)) {
          hessian <- hessian +
            crossprod(seedc[[aq]], (wc * Hf[rows0, aq, bq]) * seedc[[bq]])
        }
      }
    }
  }
  hessian <- (hessian + t(hessian)) / 2

  zeta <- vapply(nm, function(j)
    linkfunctions7::linkfun(links[[j]], v[[j]]), numeric(1))
  hessian <- .marg_fd_rows(hessian, term, eta, y, logdens, grad, zeta, seed,
                           cols, level, w, c("m1", "tau1"))
  list(loglik = loglik, gradient = gradient, hessian = hessian)
}

# the fully propagated Hessian of the one-break-point gaussian step model:
# the interval sum differentiated twice
.marg_jump1_hessian <- function(term, eta, y, logdens, grad, hess, psi, seed,
                                cols, level, w) {
  bp <- term@blueprint
  v <- .marg_check_psi(term, psi)
  m <- v[["m1"]]
  tau <- v[["tau1"]]
  delta <- v[["delta1"]]
  n <- bp$n
  npar <- length(seed)
  mm <- ncol(seed[[1L]])
  idx <- seq_len(n)
  LF0 <- as.numeric(logdens(eta, idx))
  LF1 <- as.numeric(logdens(eta + delta, idx))
  G0 <- as.matrix(grad(eta, idx))
  G1 <- as.matrix(grad(eta + delta, idx))
  H0 <- array(as.numeric(hess(eta, idx)), c(n, npar, npar))
  H1 <- array(as.numeric(hess(eta + delta, idx)), c(n, npar, npar))
  if (nrow(G0) != n || ncol(G0) != npar) {
    stop(sprintf("'grad' must return an %d by %d matrix.", n, npar),
         call. = FALSE)
  }

  # the per-observation emission scores in the caller's unknowns, one per
  # side of the break-point: the shifted side's delta column carries the
  # level equation's own score
  gm0 <- matrix(0, n, mm)
  gm1 <- matrix(0, n, mm)
  for (q in seq_len(npar)) {
    gm0 <- gm0 + G0[, q] * seed[[q]]
    gm1 <- gm1 + G1[, q] * seed[[q]]
  }
  ic <- cols[match(c("m1", "tau1", "delta1"), term_params(term))]
  gm1[, ic[3L]] <- gm1[, ic[3L]] + G1[, level]

  loglik <- numeric(n)
  gradient <- numeric(mm)
  hessian <- matrix(0, mm, mm)
  P1 <- numeric(n)

  for (rs in bp$groups) {
    ng <- length(rs)
    wi <- w[rs][1L]
    iv <- .marg_intervals(bp$x[rs], m, tau, d2 = TRUE)
    cvec <- c(0, cumsum(LF0[rs] - LF1[rs])) + sum(LF1[rs])
    a <- iv$lm + cvec
    logL <- .marg_lse(a)
    wj <- exp(a - logL)
    P1[rs] <- cumsum(wj)[seq_len(ng)]

    # the interval-level gradients, on the unconstrained scale: the
    # conditional part by the same one-observation-per-interval swap as the
    # values, the mass part chained onto tau's log chart
    M <- gm0[rs, , drop = FALSE] - gm1[rs, , drop = FALSE]
    CS <- apply(M, 2L, cumsum)
    if (is.null(dim(CS))) CS <- matrix(CS, nrow = 1L)
    base <- colSums(gm1[rs, , drop = FALSE])
    DA <- matrix(base, ng + 1L, mm, byrow = TRUE)
    DA[-1L, ] <- DA[-1L, ] + CS
    DA[, ic[1L]] <- DA[, ic[1L]] + iv$dm
    DA[, ic[2L]] <- DA[, ic[2L]] + iv$dt * tau
    gbar <- colSums(wj * DA)
    gradient <- gradient + wi * gbar
    hessian <- hessian + wi * (crossprod(DA, wj * DA) - outer(gbar, gbar))

    # the masses' own curvature, in the (m, zeta_tau) block: the log chart
    # contributes tau^2 times the second derivative plus tau times the first
    smm <- sum(wj * iv$dmm)
    smt <- sum(wj * iv$dmt) * tau
    stt <- sum(wj * iv$dtt) * tau^2 + sum(wj * iv$dt) * tau
    i1 <- ic[1L]
    i2 <- ic[2L]
    hessian[i1, i1] <- hessian[i1, i1] + wi * smm
    hessian[i1, i2] <- hessian[i1, i2] + wi * smt
    hessian[i2, i1] <- hessian[i2, i1] + wi * smt
    hessian[i2, i2] <- hessian[i2, i2] + wi * stt

    # the predictive contributions, from the same states
    A <- iv$lm
    tot <- .marg_lse(A)
    for (t in seq_len(ng)) {
      row <- rs[t]
      A <- A + c(rep(LF1[row], t), rep(LF0[row], ng + 1L - t))
      tot2 <- .marg_lse(A)
      loglik[row] <- tot2 - tot
      tot <- tot2
    }
  }

  # sum_j w_j d2c_j collapses over the intervals: each observation's second
  # derivative enters weighted by the posterior probability of its side
  seed1 <- seed
  seed1[[level]][, ic[3L]] <- seed1[[level]][, ic[3L]] + 1
  for (s in 0:1) {
    Hs <- if (s == 1L) H1 else H0
    sd_s <- if (s == 1L) seed1 else seed
    ws <- w * (if (s == 1L) P1 else 1 - P1)
    for (aq in seq_len(npar)) {
      for (bq in seq_len(npar)) {
        hessian <- hessian +
          crossprod(sd_s[[aq]], (ws * Hs[, aq, bq]) * sd_s[[bq]])
      }
    }
  }

  hessian <- (hessian + t(hessian)) / 2
  list(loglik = loglik, gradient = gradient, hessian = hessian)
}

S7::method(print, MarginalBreakTerm) <- function(x, ...) {
  built <- length(x@blueprint) > 0L
  cat(sprintf(paste0("<MarginalBreakTerm> '%s' (%s): %d latent break-point%s",
                     " per group, integrated out%s\n"), x@label, x@kind,
              x@npsi, if (x@npsi > 1L) "s" else "",
              if (built) sprintf(" (%d groups)", length(x@blueprint$groups))
              else " (specification)"))
  cat("  parameters: ", paste(term_params(x), collapse = ", "), "\n",
      sep = "")
  if (!is.null(x@prior)) {
    cat("  prior: ", attr(S7::S7_class(x@prior), "name"), "\n", sep = "")
  }
  invisible(x)
}
