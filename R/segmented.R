#' @include term_classes.R generics.R nonlinear.R
NULL

#' @title S7 Class for Break-Point Terms
#' @name SegTerm
#'
#' @description
#' A subclass of \code{\link{additive_term}} for a covariate whose effect
#' changes at estimated break-points: a change of slope
#' (\code{\link{seg}}), a change of level (\code{\link{jump}}), or both at
#' the same points (\code{\link{jseg}}). The design block is the working
#' one of the iteration that estimates the break-points, and is
#' recomputed by \code{\link{term_refresh}} as they move.
#'
#' @inheritParams additive_term
#' @param kind Which of the three constructions.
#' @param var The covariate expression.
#' @param npsi The number of break-points.
#' @param linear Whether the block carries the linear effect.
#' @param subformulas A named list of one-sided formulas, one per developed
#'   parameter.
#' @param spec The resolved construction settings.
#'
#' @return An object of class \code{SegTerm}.
#'
#' @seealso \code{\link{seg}}, \code{\link{jump}}, \code{\link{jseg}}
#' @examples
#' S7::S7_inherits(seg(x), SegTerm)
#' @export
SegTerm <- S7::new_class(
  name = "SegTerm",
  parent = additive_term,
  properties = list(
    kind = S7::class_character,
    var = S7::class_any,
    npsi = S7::class_integer,
    linear = S7::class_logical,
    subformulas = S7::class_list,
    spec = S7::class_list
  )
)

#' Segmented Term: a Broken Line with Estimated Break-Points
#'
#' @description
#' A covariate whose slope changes at \eqn{K} break-points estimated with
#' everything else, the relationship staying continuous
#' (\cite{muggeo2003}). Written in the equation of any parameter of any
#' distribution, the term contributes
#'
#' \deqn{\beta x_i + \sum_{k=1}^{K} \gamma_k\,(x_i - \psi_k)_{+},
#'   \qquad (u)_{+} = \max(0, u),}
#'
#' to that parameter's linear predictor, so that in a model
#' \eqn{g(\theta_i) = \eta_i} with the rest of the equation supplying
#' \eqn{z_i'\alpha},
#'
#' \deqn{\eta_i = z_i'\alpha + \beta x_i
#'   + \sum_{k=1}^{K} \gamma_k\,(x_i - \psi_k)_{+}.}
#'
#' \eqn{\beta} is the slope before the first break-point and
#' \eqn{\gamma_k} the change of slope at \eqn{\psi_k}, so the slope over
#' the \eqn{k}-th segment is \eqn{\beta + \gamma_1 + \cdots + \gamma_k}.
#' The gaussian additive model is the case \eqn{g = \mathrm{id}} on the
#' mean of a \code{gaussian1_distrib}; nothing here is particular to it.
#'
#' @details
#' \subsection{How the break-points are estimated}{
#' The contribution is differentiable in \eqn{\psi_k} away from the
#' break-point, with
#' \eqn{\partial/\partial\psi_k = -\gamma_k\,\mathbb{1}(x>\psi_k)}, so the
#' design block is the ordinary Jacobian, the break-point is an ordinary
#' coefficient, and the increment a Gauss-Newton step solves for is the
#' update of \cite{muggeo2003}: his working variables \eqn{U} and \eqn{V}
#' are the block's columns and his \eqn{\psi \leftarrow \psi +
#' \gamma/\delta} is that step. Refreshing the term at the current
#' coefficients (\code{\link{term_refresh}}) and fitting is one iteration
#' of it, which is the same contract \code{\link{nl}} uses.
#'
#' The objective has local optima in the break-points and the iteration
#' converges from within a basin around where it starts, so
#' \code{\link{seg_start}} chooses the starting positions on a grid rather
#' than at a conventional point. A break-point is confined to the interval
#' between the 5th and the 95th percentile of the covariate: outside it
#' the indicator is constant, the truncated line and that constant are
#' linearly dependent, and the block is singular rather than merely
#' ill-conditioned. A run ending against the limit has not located a
#' break-point, and \code{\link{seg_psi}} then reports the limit itself.
#'
#' The iteration can settle into a cycle of period two in the break-point,
#' in which case the rule of \code{\link{seg_converged}} is never met while
#' the objective has long since stopped moving; a caller that can evaluate
#' the objective should stop on its relative change, as \code{segmented}
#' does.
#' }
#'
#' @section The coefficients:
#' \code{beta} (present when \code{linear = TRUE}), \code{gamma1} \ldots
#' \code{gammaK} and \code{psi1} \ldots \code{psiK}, in that order. The
#' break-points are coefficients of the working block here, which they are
#' not in the discontinuous constructions, so \code{\link{seg_psi}} is the
#' function that reports them either way.
#'
#' @section Developing a coefficient with covariates:
#' Every coefficient of the term is a parameter like any other and may be
#' developed on covariates, through a two-sided formula in \code{...} whose
#' left side names it:
#'
#' \preformatted{seg(x, psi ~ id)                   a break-point per subject
#' seg(x, psi ~ random(~1 | id))      the random-changepoint model
#' seg(x, gamma1 ~ z)                 a first slope change varying with z
#' seg(x, by = ~0 + g)                every coefficient, per level of g}
#'
#' The right side goes through \code{\link{interpret_formula}}, so it takes
#' \emph{any} term of this package, penalized ones included, and the
#' hyperparameters a sub-term declares are reported by
#' \code{\link{term_penalties}} and estimated by the fitting layer. The
#' left side names either one coefficient (\code{gamma2}) or a whole kind
#' (\code{gamma}, meaning every \eqn{\gamma_k}, each with coefficients of
#' its own over the shared design). \code{by = ~f} is the shorthand giving
#' every coefficient the same development, and does not combine with
#' per-coefficient formulas.
#'
#' A developed coefficient carries one value per observation:
#' \eqn{\gamma_{k,i} = w_i'a_k} and \eqn{\psi_{k,i} =
#' w_i'p_k}, with \eqn{w_i} the row of the sub-design. Both
#' enter the block linearly, the truncated line becoming
#' \eqn{(x_i-\psi_{k,i})_{+}w_i'} and the Jacobian column
#' \eqn{-\gamma_{k,i}\,\mathbb{1}(x_i>\psi_{k,i})\,w_i'}, so the
#' continuous construction accepts a development of any coefficient and any
#' combination of them.
#'
#' A penalty on the changes themselves -- the lasso that selects how many
#' break-points are really there -- is the development on an intercept
#' alone, \code{seg(x, npsi = 4, gamma ~ 0 + lasso(~1))}: the \eqn{K}
#' changes are then one penalized block under one hyperparameter, the
#' entries of a subformula shared by every coefficient of a kind being
#' pooled into one. The \code{0 +} is not decoration. A subformula follows
#' the intercept convention of a formula, so \code{gamma ~ lasso(~1)}
#' carries the subformula's own unpenalized intercept beside the penalized
#' column, which for an intercept-only development is the same column
#' twice.
#'
#' @section The linear effect:
#' With \code{linear = TRUE}, the default, the term carries \eqn{\beta x}
#' itself, so \code{y ~ seg(x)} is the whole relationship. Writing the
#' covariate again in the same equation is then exactly collinear with that
#' column, and \code{\link{interpret_formula}} removes it with a warning.
#' \code{linear = FALSE} is the explicit way to keep the linear effect
#' outside the term, and \code{y ~ x + seg(x, linear = FALSE)} is the same
#' model written the other way round.
#'
#' @param x The covariate, an expression evaluated in the data.
#' @param ... Two-sided formulas developing the term's own coefficients;
#'   see the section above. Cannot be combined with \code{by}.
#' @param npsi The number of break-points. Defaults to 1.
#' @param psi Optional starting positions; defaults to evenly spaced
#'   quantiles of the covariate. Where a break-point carries a development
#'   they seed it, each starting vector solving \eqn{Wp_k \approx
#'   \psi_k^{0}}.
#' @param by A one-sided formula giving every coefficient of the term the
#'   same development, e.g. \code{by = ~0 + g} for an independent set per
#'   level of a factor. A bare variable is rejected; write the formula.
#' @param linear Whether the block carries the linear effect \eqn{\beta x}.
#'   Defaults to \code{TRUE}.
#' @param label A single non-empty string prefixed to the coefficient
#'   names.
#'
#' @return An object of class \code{\link{SegTerm}} (a specification; see
#'   \code{\link{term_build}}).
#'
#' @references
#' Muggeo, V. M. R. (2003). Estimating regression models with unknown
#' break-points. \emph{Statistics in Medicine}, 22(19), 3055--3071.
#'
#' Muggeo, V. M. R., Atkins, D. C., Gallop, R. J. and Dimidjian, S. (2014).
#' Segmented mixed models with random changepoints: a maximum likelihood
#' approach with application to treatment for depression study.
#' \emph{Statistical Modelling}, 14(4), 293--313.
#'
#' @examples
#' set.seed(1)
#' dd <- data.frame(x = sort(runif(200, 0, 10)))
#' dd$y <- 1 + 0.5 * dd$x + 2 * pmax(dd$x - 6, 0) + rnorm(200, sd = 0.3)
#'
#' built <- term_build(seg(x), dd)
#' term_coef_names(built)
#' seg_psi(built)
#'
#' @seealso \code{\link{jump}}, \code{\link{jseg}}, \code{\link{seg_psi}},
#'   \code{\link{seg_start}}, \code{\link{seg_step}}, \code{\link{nl}}
#' @export
seg <- function(x, ..., npsi = 1, psi = NULL, by = NULL, linear = TRUE,
                label = "seg") {
  # a continuous term has no scaling factor, its working block being
  # bounded already; the value is carried so that one blueprint serves the
  # three constructions
  .seg_spec("seg", substitute(x), npsi, psi, linear, 0.05, label,
            list(...), .seg_by_formula(substitute(by), parent.frame()))
}

#' Stepmented Term: a Level that Changes at Estimated Break-Points
#'
#' @description
#' A covariate whose effect is a constant that steps to a new value at
#' \eqn{K} break-points estimated with everything else
#' (\cite{fasola2018}). Written in the equation of any parameter of any
#' distribution, the term contributes
#'
#' \deqn{\sum_{k=1}^{K} \delta_k\,\mathbb{1}(x_i \geq \psi_k)}
#'
#' to that parameter's linear predictor, so that in a model
#' \eqn{g(\theta_i) = \eta_i} with the rest of the equation supplying
#' \eqn{z_i'\alpha},
#'
#' \deqn{\eta_i = z_i'\alpha
#'   + \sum_{k=1}^{K} \delta_k\,\mathbb{1}(x_i \geq \psi_k).}
#'
#' There is no linear effect: the relationship between \eqn{x} and the
#' predictor is a step function, the intercept of the equation being its
#' level before the first break-point and \eqn{\delta_k} the change of
#' level at \eqn{\psi_k}. A model that is linear in \eqn{x} \emph{and}
#' steps is written by putting the covariate in the equation as well,
#' \code{y ~ x + jump(x)}, or by \code{\link{jseg}} where the slope is to
#' change at the same points.
#'
#' @details
#' \subsection{How the break-points are estimated}{
#' A step is not differentiable in \eqn{\psi} at all, and yet needs no
#' search over candidate positions. The identity
#' \deqn{\mathbb{1}(x>\psi) = \frac{1}{2}
#'   + \frac{x-\psi}{2\lvert x-\psi \rvert}}
#' holds exactly for \eqn{x \neq \psi} and is linear in \eqn{\psi} once the
#' weight \eqn{1/(2\lvert x-\psi^{0}\rvert)} is held at the previous
#' iterate. Writing \eqn{W = 1/(2\lvert x-\psi^{0}\rvert)} and \eqn{Z = xW
#' + 1/2}, a step of size \eqn{\delta} at \eqn{\psi} is \eqn{\delta Z + gW}
#' with \eqn{g = -\delta\psi}, so a linear fit on \eqn{(Z, W)} returns the
#' break-point as \eqn{\psi = -g/\delta} (\cite{fasola2018}). The
#' break-point is therefore not a coefficient here but a quantity read off
#' two of them, which is why refreshing the term recovers it before
#' rebuilding the weights, and why \code{\link{seg_psi}} exists.
#' }
#'
#' \subsection{The scaling schedule}{
#' The weight is unbounded as \eqn{x} approaches \eqn{\psi}, and since
#' \eqn{Z - \psi W} is the step itself, \eqn{Z} is \eqn{\psi W} plus a
#' quantity of order one: an unbounded \eqn{W} makes the two columns
#' numerically collinear and drowns the signal the fit reads. The remedy of
#' \cite{fasola2018} is to move the observations rather than to cap the
#' weight. With a scaling factor \eqn{c} the two intervals \eqn{[x_{(1)},
#' \psi]} and \eqn{(\psi, x_{(n)}]} are mapped onto
#' \deqn{[x_{(1)},\, \psi - c(\psi - x_{(1)})]
#'   \quad\text{and}\quad
#'   (\psi + c(x_{(n)} - \psi),\, x_{(n)}],}
#' which leaves a gap of relative width \eqn{c} around \eqn{\psi} and
#' bounds \eqn{W} without altering the model: the working covariates are
#' computed on the rescaled covariate, while the reported contribution
#' stays on the original one.
#'
#' The factor is not a constant. It governs how far the break-point may
#' travel in one step, so a large value lets the estimate leave a spurious
#' optimum and a small one is faithful to the step function. \code{c0} is
#' its starting value and \code{\link{term_refresh}} halves it whenever the
#' break-point reverses direction, the signal that the iteration has begun
#' to circle an optimum rather than travel towards one. The run has
#' converged when every break-point moves less than a hundredth of the
#' distance between consecutive distinct observations, which
#' \code{\link{seg_converged}} reports.
#'
#' A break-point is confined to the interval between the 5th and the 95th
#' percentile of the covariate, outside which the indicator is constant and
#' the block singular, and the starting positions are chosen on a grid by
#' \code{\link{seg_start}}: measured over eight samples and four starting
#' positions, the fraction of runs recovering the break-point is 0 to 0.5
#' from a single conventional start and 1 from the grid.
#' }
#'
#' @section The coefficients:
#' \code{delta1} \ldots \code{deltaK}, the changes of level, followed by
#' \code{g1} \ldots \code{gK}, the auxiliary coefficients of the identity
#' above, from which \eqn{\psi_k = -g_k/\delta_k} is read. The
#' break-points are \emph{not} coefficients; ask \code{\link{seg_psi}}.
#'
#' @section Developing a coefficient with covariates:
#' As for \code{\link{seg}}, a two-sided formula in \code{...} whose left
#' side names \code{delta}, \code{deltaK} or \code{psi}, \code{psiK}
#' develops that coefficient on covariates, and \code{by = ~f} gives every
#' coefficient the same development. The right side goes through
#' \code{\link{interpret_formula}} and takes any term of this package.
#'
#' What the discontinuous construction can carry is narrower than what the
#' continuous one can, and the reason is the read-off rather than the
#' model. The auxiliary coefficient is \eqn{g_k = -\delta_k\psi_k}, a
#' product of the two quantities, so it stays linear in the unknowns only
#' when one of the two factors is a single number:
#'
#' \itemize{
#'   \item \strong{the break-point developed}, the change of level a single
#'     number: exact for any design, \eqn{g_k = -\delta_k p_k} being
#'     linear in \eqn{p_k}. \code{jump(x, psi ~ id)} is a break-point
#'     per subject with a shared step size.
#'   \item \strong{a development on group indicators}, where the design has
#'     one column per group and each observation belongs to one: the
#'     product collapses group by group and the read-off is exact within
#'     each. \code{jump(x, by = ~0 + g)} is an independent step and
#'     break-point per level.
#'   \item \strong{anything else} is rejected rather than approximated,
#'     the product of two developments needing the outer product of their
#'     designs and no unconstrained fit returning it as one.
#' }
#'
#' A sub-term carrying a penalty is rejected on a developed break-point,
#' since the estimated coefficients are \eqn{-\delta_k p_k}, the
#' development scaled by the step size, and a penalty would act on that
#' rather than on the development. A penalty on the changes of level
#' themselves is \code{jump(x, npsi = 4, delta ~ 0 + lasso(~1))}, the
#' \code{0 +} removing the subformula's own unpenalized intercept, which
#' would otherwise be the same column twice.
#'
#' @inheritParams seg
#' @param c0 The starting value of the scaling factor that separates the
#'   observations from the break-point, as a fraction of the distance to
#'   the ends of the range. Defaults to \code{0.05}, the value
#'   \cite{fasola2018} recommend; see Details.
#'
#' @return An object of class \code{\link{SegTerm}} (a specification; see
#'   \code{\link{term_build}}).
#'
#' @references
#' Fasola, S., Muggeo, V. M. R. and Kuchenhoff, H. (2018). A heuristic,
#' iterative algorithm for change-point detection in abrupt change
#' models. \emph{Computational Statistics}, 33, 997--1015.
#'
#' @examples
#' set.seed(1)
#' dd <- data.frame(x = sort(runif(200, 0, 10)))
#' dd$y <- 1 + 2 * (dd$x > 6) + rnorm(200, sd = 0.3)
#'
#' built <- term_build(jump(x), dd)
#' term_coef_names(built)
#' seg_psi(built)
#'
#' @seealso \code{\link{seg}}, \code{\link{jseg}}, \code{\link{seg_psi}},
#'   \code{\link{seg_start}}, \code{\link{seg_step}}
#' @export
jump <- function(x, ..., npsi = 1, psi = NULL, by = NULL, c0 = 0.05,
                 label = "jump") {
  .seg_spec("jump", substitute(x), npsi, psi, FALSE, c0, label,
            list(...), .seg_by_formula(substitute(by), parent.frame()))
}

#' Segmented-with-Jump Term: Slope and Level Both Changing
#'
#' @description
#' A covariate whose slope \emph{and} level both change at \eqn{K}
#' break-points estimated with everything else, the union of
#' \code{\link{seg}} and \code{\link{jump}} at the same points. Written in
#' the equation of any parameter of any distribution, the term contributes
#'
#' \deqn{\beta x_i + \sum_{k=1}^{K}
#'   \bigl[\delta_k\,\mathbb{1}(x_i \geq \psi_k)
#'         + \gamma_k\,(x_i - \psi_k)_{+}\bigr]}
#'
#' to that parameter's linear predictor, so that in a model
#' \eqn{g(\theta_i) = \eta_i} with the rest of the equation supplying
#' \eqn{z_i'\alpha},
#'
#' \deqn{\eta_i = z_i'\alpha + \beta x_i + \sum_{k=1}^{K}
#'   \bigl[\delta_k\,\mathbb{1}(x_i \geq \psi_k)
#'         + \gamma_k\,(x_i - \psi_k)_{+}\bigr].}
#'
#' \eqn{\beta} is the slope before the first break-point, \eqn{\gamma_k}
#' the change of slope at \eqn{\psi_k} and \eqn{\delta_k} the change of
#' level there.
#'
#' @details
#' Both devices of the two constructions are used at once: the truncated
#' line is differentiable in the break-point and the step is linearized by
#' the identity of \cite{fasola2018}, on the rescaled covariate its scaling
#' schedule provides (see \code{\link{jump}}). What is particular to the
#' joint term is the read-off, since the truncated line depends on the
#' break-point as well and reading only the step's pair would discard what
#' the change of slope says about it. Linearizing both parts about the
#' previous position, and writing \eqn{h} for the increment, \eqn{U} for
#' the truncated line and \eqn{I = Z - \psi^{0}W} for the indicator,
#'
#' \deqn{\gamma\,(x-\psi)_{+} \approx \gamma U - \gamma h I,
#'   \qquad \delta\,\mathbb{1}(x>\psi)
#'     = \delta Z - \delta\psi W,}
#'
#' so the fitted coefficients of \eqn{Z} and \eqn{W} are \eqn{a = \delta -
#' \gamma h} and \eqn{b = \gamma h\psi^{0} - \delta\psi^{0} - \delta h},
#' and the increment solves
#'
#' \deqn{\gamma h^{2} + a h + (b + a\psi^{0}) = 0,}
#'
#' the root of smaller modulus. At \eqn{\gamma = 0} the quadratic
#' degenerates to \eqn{h = -(b + a\psi^{0})/a}, that is \eqn{\psi = -b/a},
#' so the pure step is the case this contains rather than an exception to
#' it.
#'
#' Everything \code{\link{jump}} documents about the scaling schedule, the
#' confinement of a break-point and the choice of starting positions
#' applies here unchanged.
#'
#' @section The coefficients:
#' \code{beta} (present when \code{linear = TRUE}), \code{gamma1} \ldots
#' \code{gammaK}, the changes of slope, then \code{delta1} \ldots
#' \code{deltaK}, the changes of level, then \code{g1} \ldots \code{gK},
#' the auxiliary coefficients from which the break-points are read. Ask
#' \code{\link{seg_psi}} for the break-points.
#'
#' @section Developing a coefficient with covariates:
#' \code{beta} and \code{gammaK} take a development on any design, both
#' entering the block linearly. \code{deltaK} and \code{psiK} enter the
#' quadratic above, which is a scalar identity: it splits observation by
#' observation only where the design has one column per group and each
#' observation belongs to one, so \code{jseg(x, by = ~0 + g)} is an
#' independent set of everything per level and a development of
#' \code{delta} or \code{psi} on any other design is rejected. The
#' componentwise reading that would remain diverges whenever a change of
#' level passes near zero mid-iteration, which on a joint model it
#' routinely does; it is rejected rather than shipped diverging.
#'
#' @inheritParams seg
#' @inheritParams jump
#'
#' @return An object of class \code{\link{SegTerm}} (a specification; see
#'   \code{\link{term_build}}).
#'
#' @references
#' Muggeo, V. M. R. (2003). Estimating regression models with unknown
#' break-points. \emph{Statistics in Medicine}, 22(19), 3055--3071.
#'
#' Fasola, S., Muggeo, V. M. R. and Kuchenhoff, H. (2018). A heuristic,
#' iterative algorithm for change-point detection in abrupt change
#' models. \emph{Computational Statistics}, 33, 997--1015.
#'
#' @examples
#' set.seed(1)
#' dd <- data.frame(x = sort(runif(200, 0, 10)))
#' dd$y <- 1 + 0.3 * dd$x + 2 * (dd$x > 6) +
#'   1.5 * pmax(dd$x - 6, 0) + rnorm(200, sd = 0.3)
#'
#' built <- term_build(jseg(x), dd)
#' term_coef_names(built)
#'
#' @seealso \code{\link{seg}}, \code{\link{jump}}, \code{\link{seg_psi}},
#'   \code{\link{seg_start}}, \code{\link{seg_step}}
#' @export
jseg <- function(x, ..., npsi = 1, psi = NULL, by = NULL, linear = TRUE,
                 c0 = 0.05, label = "jseg") {
  .seg_spec("jseg", substitute(x), npsi, psi, linear, c0, label,
            list(...), .seg_by_formula(substitute(by), parent.frame()))
}

# `by` is a one-sided formula and not a variable: the shorthand it stands
# for develops every coefficient, and a bare variable would have to be
# read as `~0 + g`, which is a guess about the coding. A symbol is
# rejected with the formula it means, since the value would otherwise be
# looked up in the calling frame and fail there, three frames from the
# mistake.
.seg_by_formula <- function(by, env) {
  if (is.null(by)) return(NULL)
  if (is.call(by) && identical(by[[1L]], as.name("~"))) {
    if (length(by) != 2L) {
      stop("a formula 'by' must be one-sided, e.g. by = ~0 + g.",
           call. = FALSE)
    }
    return(stats::as.formula(by, env = env))
  }
  stop(sprintf(paste("'by' must be a one-sided formula giving every",
                     "coefficient of the term the same development;",
                     "write by = ~0 + %s for an independent set per level."),
               deparse(by)), call. = FALSE)
}

# The term's own coefficients, in the order the block carries them. They
# are what a subformula may name, and `beta` is among them because the
# term owns the linear effect: the covariate cannot appear again in the
# same equation, so a linear effect varying with a covariate has nowhere
# else to be written.
.seg_params <- function(kind, npsi, linear) {
  k <- seq_len(npsi)
  c(if (linear) "beta" else character(0),
    if (kind %in% c("seg", "jseg")) paste0("gamma", k) else character(0),
    if (kind %in% c("jump", "jseg")) paste0("delta", k) else character(0),
    paste0("psi", k))
}

# The stems a subformula may name in place of a numbered coefficient,
# meaning every coefficient of that kind.
.seg_stems <- function(kind, npsi, linear) {
  c(if (linear) "beta" else character(0),
    if (kind %in% c("seg", "jseg")) "gamma" else character(0),
    if (kind %in% c("jump", "jseg")) "delta" else character(0),
    "psi")
}

# The subformulas of a break-point term, resolved to one entry per
# coefficient: two-sided formulas in `...` whose left side names a
# coefficient or a stem, or `by = ~f` giving every coefficient the same
# development. The two forms do not mix, `by` being the shorthand for all
# of what the other spells out.
.seg_gather <- function(dots, by, kind, npsi, linear) {
  params <- .seg_params(kind, npsi, linear)
  stems <- .seg_stems(kind, npsi, linear)
  out <- list()

  if (!is.null(by)) {
    if (length(dots)) {
      stop(paste("'by = ~f' gives every coefficient the same development;",
                 "it does not combine with per-coefficient formulas."),
           call. = FALSE)
    }
    if (!inherits(by, "formula") || length(by) != 2L) {
      stop(paste("'by' must be a one-sided formula giving every coefficient",
                 "the same development, e.g. by = ~0 + g."), call. = FALSE)
    }
    for (p in params) out[[p]] <- by
    return(out)
  }

  if (!length(dots)) return(out)
  nms <- names(dots)
  for (i in seq_along(dots)) {
    d <- dots[[i]]
    # a named argument in `...` is a mistyped or removed one, and is
    # reported as such rather than swallowed by the formula check
    if (!is.null(nms) && nzchar(nms[i])) {
      stop(sprintf("'%s' is not an argument of a break-point term.", nms[i]),
           call. = FALSE)
    }
    if (!inherits(d, "formula") || length(d) != 3L || !is.name(d[[2L]])) {
      stop(paste("arguments in '...' must be two-sided formulas developing",
                 "a coefficient of the term, e.g. psi ~ g."), call. = FALSE)
    }
    lhs <- as.character(d[[2L]])
    rhs <- stats::as.formula(call("~", d[[3L]]), env = environment(d))
    targets <- if (lhs %in% params) lhs else if (lhs %in% stems) {
      grep(sprintf("^%s[0-9]+$", lhs), params, value = TRUE)
    } else {
      stop(sprintf(paste("'%s' is not a coefficient of a %s term; the",
                         "coefficients are %s, and a stem (%s) names every",
                         "one of a kind."),
                   lhs, kind, paste(params, collapse = ", "),
                   paste(setdiff(stems, params), collapse = ", ")),
           call. = FALSE)
    }
    for (p in targets) {
      if (!is.null(out[[p]])) {
        stop(sprintf("'%s' is developed twice.", p), call. = FALSE)
      }
      out[[p]] <- rhs
    }
  }
  out[intersect(params, names(out))]
}

.seg_spec <- function(kind, var, npsi, psi, linear, c0, label, dots, by) {
  if (!is.numeric(c0) || length(c0) != 1L || is.na(c0) ||
      c0 <= 0 || c0 >= 1) {
    stop("'c0' must be a single number strictly between 0 and 1.",
         call. = FALSE)
  }
  if (!is.numeric(npsi) || length(npsi) != 1L || is.na(npsi) || npsi < 1 ||
      npsi != round(npsi)) {
    stop("'npsi' must be a single integer of at least 1.", call. = FALSE)
  }
  if (!is.null(psi) && (!is.numeric(psi) || length(psi) != npsi ||
                        anyNA(psi))) {
    stop("'psi' must give one starting position per break-point.",
         call. = FALSE)
  }
  if (!is.logical(linear) || length(linear) != 1L || is.na(linear)) {
    stop("'linear' must be TRUE or FALSE.", call. = FALSE)
  }
  if (!is.character(label) || length(label) != 1L || is.na(label) ||
      !nzchar(label)) {
    stop("'label' must be a single non-empty character string.",
         call. = FALSE)
  }
  npsi <- as.integer(npsi)
  subs <- .seg_gather(dots, by, kind, npsi, linear)
  SegTerm(label = label, kind = kind, var = var, npsi = npsi,
          linear = linear, subformulas = subs,
          spec = list(psi = psi, c0 = c0),
          X = NULL, coef_names = character(0),
          blueprint = list(), penalty = NULL)
}

# The stem of each own coefficient's name. It is the parameter's own name
# except for a break-point of a discontinuous construction, whose
# coefficients are the auxiliary g of the identity of Fasola et al. and not
# the position: the position is read off, and seg_psi() is what reports it.
.seg_coef_stem <- function(kind, p) {
  if (kind != "seg" && grepl("^psi[0-9]+$", p)) sub("^psi", "g", p) else p
}

# A design with one column per group and every observation in exactly one
# of them. The product of two developments over such a design collapses
# group by group, which is what lets a discontinuous construction carry
# both a change of level and a break-point that vary.
.seg_is_partition <- function(Z) {
  Z <- as.matrix(Z)
  all(Z == 0 | Z == 1) && all(rowSums(Z) == 1)
}

# What each construction can carry, and why. The continuous one accepts a
# development of anything: every coefficient enters its block linearly. The
# discontinuous ones read the break-point off g = -delta * psi, a product
# of two of the unknowns, so it stays linear only where one factor is a
# single number or where the design collapses the product; and the joint
# construction reads it off a quadratic that also carries the change of
# slope, so that one is constrained as well.
.seg_check_devs <- function(kind, npsi, dev) {
  if (kind == "seg") return(invisible(NULL))
  for (k in seq_len(npsi)) {
    pp <- paste0("psi", k)
    involved <- c(if (kind == "jseg") paste0("gamma", k), paste0("delta", k), pp)
    # a development of one column is a single number under another name and
    # constrains nothing; what has to collapse is a coefficient that really
    # varies
    wide <- involved[vapply(involved, function(p) {
      !is.null(dev[[p]]) && ncol(dev[[p]]$Z) > 1L
    }, logical(1))]
    if (!length(wide)) next
    # a pure step reads g = -delta * psi, which stays linear in the
    # development whenever the change of level is a single number
    if (kind == "jump" && identical(wide, pp)) next
    keys <- unique(vapply(wide, function(p) dev[[p]]$key, character(1)))
    if (length(keys) == 1L && pp %in% wide && isTRUE(dev[[pp]]$partition)) next
    stop(sprintf(paste("'%s' is developed in a %s term, whose break-point is",
                       "read off a product of the unknowns. That splits",
                       "observation by observation only when every one of %s",
                       "that varies carries the SAME development and its",
                       "design has one column per group with each observation",
                       "in one of them, as by = ~0 + g does.%s"),
                 wide[1L], kind,
                 paste(sprintf("'%s'", involved), collapse = ", "),
                 if (kind == "jump")
                   " Developing the break-point alone is exact on any design."
                 else " Use seg, or develop nothing."),
         call. = FALSE)
  }
  invisible(NULL)
}

# The sub-design of every developed coefficient, its built sub-terms and
# the penalties they declare. A coefficient with no subformula has no
# design: it is a single number, which the block reads as a constant column
# without forming one.
.seg_designs <- function(term, data) {
  params <- .seg_params(term@kind, term@npsi, term@linear)
  dev <- stats::setNames(vector("list", length(params)), params)
  for (p in params) {
    sf <- term@subformulas[[p]]
    if (is.null(sf)) next
    sub <- .nl_submodel(p, sf, data)
    dev[[p]] <- list(Z = sub$Z, subs = sub$terms, pens = sub$penalties,
                     coef_names = sub$coef_names,
                     key = paste(deparse(sf), collapse = ""),
                     partition = .seg_is_partition(sub$Z))
  }
  .seg_check_devs(term@kind, term@npsi, dev)
  if (term@kind != "seg") {
    for (k in seq_len(term@npsi)) {
      d <- dev[[paste0("psi", k)]]
      if (!is.null(d) && length(d$pens)) {
        stop(paste("a sub-term carrying a penalty is rejected on a developed",
                   "break-point of a discontinuous construction: the",
                   "estimated coefficients are -delta * gamma, the",
                   "development scaled by the change of level, so the",
                   "penalty would act on that rather than on the",
                   "development."), call. = FALSE)
      }
    }
  }
  dev
}

# The value of one developed coefficient at every observation, and the
# positions of one break-point. A scalar coefficient is recycled rather
# than multiplied into a column of ones.
.seg_pval <- function(bp, coef, p, n) {
  z <- bp$Z[[p]]
  b <- coef[bp$index[[p]]]
  if (is.null(z)) rep(b, n) else as.numeric(z %*% b)
}

# The coefficient vector of a break-point's development, or the position
# itself where it carries none. A continuous term holds it directly; the
# discontinuous ones read it off, from g = -delta * psi for a pure step and
# from the quadratic of the joint construction, both componentwise, which
# is exact under the conditions .seg_check_devs() enforces.
.seg_pk <- function(bp, coef, k, prev = NULL) {
  cf <- coef[bp$index[[paste0("psi", k)]]]
  if (bp$kind == "seg") return(cf)
  dk <- coef[bp$index[[paste0("delta", k)]]]
  dk[abs(dk) < 1e-12] <- 1e-12
  if (bp$kind == "jseg" && !is.null(prev)) {
    gk <- coef[bp$index[[paste0("gamma", k)]]]
    cc <- cf + dk * prev
    disc <- dk^2 - 4 * gk * cc
    lin_h <- -cc / dk
    r <- sqrt(pmax(disc, 0))
    h1 <- (-dk + r) / (2 * gk)
    h2 <- (-dk - r) / (2 * gk)
    h <- ifelse(abs(h1) <= abs(h2), h1, h2)
    ok <- disc >= 0 & abs(gk) > 1e-10 & is.finite(h)
    return(prev + ifelse(ok, h, lin_h))
  }
  -cf / dk
}

# The positions of every break-point, one column each, and the coefficient
# vectors they were read from. Confining the position is what keeps the
# block non-singular: outside the data the indicator is constant and the
# columns are linearly dependent, not merely ill-conditioned.
.seg_positions <- function(bp, coef, n, Z = bp$Z, prev = bp$pk) {
  K <- bp$npsi
  pk <- vector("list", K)
  psi <- matrix(0, n, K)
  for (k in seq_len(K)) {
    pk[[k]] <- .seg_pk(bp, coef, k, prev[[k]])
    z <- Z[[paste0("psi", k)]]
    v <- if (is.null(z)) rep(pk[[k]], n) else as.numeric(z %*% pk[[k]])
    psi[, k] <- pmin(pmax(v, bp$lim[1L]), bp$lim[2L])
  }
  list(psi = psi, pk = pk)
}

# The rescaled covariate of Fasola, Muggeo and Kuchenhoff: the two
# intervals on either side of the break-point are shrunk towards the
# ends of the range, leaving a gap of width c around psi.
.seg_rescale <- function(xv, psi, cs, lo, hi) {
  keep <- 1 - cs
  ifelse(xv <= psi, lo + (xv - lo) * keep,
         (psi + cs * (hi - psi)) + (xv - psi) * keep)
}

# The working block and the contribution. Every coefficient of the term is
# a multiplier times its own design, and which multiplier is the whole of
# the difference between the three constructions:
#
#   beta     x
#   gamma_k  (x - psi_k)_+
#   delta_k  x~ W + 1/2                (the identity of Fasola et al.)
#   psi_k    -gamma_k(x) 1(x > psi_k)  continuous, the ordinary Jacobian
#            W                         discontinuous, the auxiliary column
#
# A coefficient with no development has a design of one constant column,
# which is never formed: its block is the multiplier itself.
.seg_block <- function(bp, xv, coef, cscale, Z = bp$Z, prev = bp$pk) {
  n <- length(xv)
  K <- bp$npsi
  pos <- .seg_positions(bp, coef, n, Z, prev)
  psi <- pos$psi
  value <- numeric(n)
  cols <- list()
  put <- function(p, mult) {
    z <- Z[[p]]
    cols[[p]] <<- if (is.null(z)) matrix(mult, ncol = 1L) else mult * z
  }
  val <- function(p) .seg_pval(list(Z = Z, index = bp$index), coef, p, n)

  if (bp$linear) {
    put("beta", xv)
    value <- value + val("beta") * xv
  }
  gam <- vector("list", K)
  if (bp$kind %in% c("seg", "jseg")) {
    for (k in seq_len(K)) {
      gam[[k]] <- val(paste0("gamma", k))
      tr <- pmax(xv - psi[, k], 0)
      put(paste0("gamma", k), tr)
      value <- value + gam[[k]] * tr
    }
  }
  if (bp$kind == "seg") {
    for (k in seq_len(K)) {
      put(paste0("psi", k), -gam[[k]] * (xv > psi[, k]))
    }
  } else {
    W <- vector("list", K)
    for (k in seq_len(K)) {
      del <- val(paste0("delta", k))
      xs <- .seg_rescale(xv, psi[, k], cscale[k], bp$lo, bp$hi)
      W[[k]] <- 1 / (2 * abs(xs - psi[, k]))
      put(paste0("delta", k), xs * W[[k]] + 0.5)
      value <- value + del * (xv > psi[, k])
    }
    for (k in seq_len(K)) put(paste0("psi", k), W[[k]])
  }

  list(X = .nl_bind(cols[bp$params]), value = value, psi = psi, pk = pos$pk)
}

# The same block through the compiled kernel, which is taken wherever no
# coefficient carries a development. That is the ordinary case and the one
# the kernel was written for; it measures 1.2 to 3.2 times the R form over
# n from 1e3 to 1e6 and one to five break-points. .seg_block stays the twin
# the tests compare it against, at a tolerance rather than by identity, a
# compiler being free to contract a multiply-add.
.seg_block_cpp <- function(bp, xv, coef, cscale, prev = bp$pk) {
  n <- length(xv)
  K <- bp$npsi
  pos <- .seg_positions(bp, coef, n, bp$Z, prev)
  off <- if (bp$linear) 1L else 0L
  gam <- if (bp$kind %in% c("seg", "jseg")) coef[off + seq_len(K)] else numeric(0)
  d <- if (bp$kind == "jseg") K else 0L
  del <- if (bp$kind == "seg") numeric(0) else coef[off + d + seq_len(K)]
  # the positions are read from the coefficient vectors rather than off the
  # matrix of values, which has no rows to read where the subset is empty
  psi1 <- pmin(pmax(vapply(pos$pk, function(v) v[1L], numeric(1)),
                    bp$lim[1L]), bp$lim[2L])
  out <- seg_block_cpp(switch(bp$kind, seg = 0L, jump = 1L, jseg = 2L),
                       xv, psi1, gam, del, cscale,
                       if (bp$linear) coef[1L] else 0, bp$linear,
                       bp$lo, bp$hi)
  list(X = out$X, value = out$value, psi = pos$psi, pk = pos$pk)
}

.seg_assemble <- function(bp, xv, coef, cscale = bp$cscale, Z = bp$Z,
                          prev = bp$pk) {
  if (bp$developed) .seg_block(bp, xv, coef, cscale, Z, prev)
  else .seg_block_cpp(bp, xv, coef, cscale, prev)
}

.seg_x <- function(expr, data) {
  v <- as.numeric(eval(expr, data, baseenv()))
  if (length(v) != nrow(data)) {
    stop("the covariate of a break-point term must give one value per row.",
         call. = FALSE)
  }
  if (anyNA(v)) {
    stop("the covariate of a break-point term must not contain missing values.",
         call. = FALSE)
  }
  v
}

# The starting coefficients of one development: the vector solving
# W p ~ psi0, the constant start projected onto the design. The solve is
# exact first -- a regularized one displaces the start by its own ridge,
# which a development of ~1 must not inherit -- and the regularized
# fallback answers for a rank-deficient design. A branch decided by a
# tryCatch is the shape the toolkit distrusts, and it is admissible here
# for the reason distrib_start() records: a starting value is allowed to be
# approximate, where a contract quantity is not.
.seg_proj_start <- function(Z, v) {
  q <- ncol(Z)
  zz <- Matrix::crossprod(Z)
  zy <- Matrix::crossprod(Z, rep(v, nrow(Z)))
  out <- tryCatch(as.numeric(Matrix::solve(zz, zy)),
                  error = function(e) NULL, warning = function(w) NULL)
  if (is.null(out) || !all(is.finite(out))) {
    sc <- max(Matrix::diag(zz), 1)
    out <- as.numeric(Matrix::solve(zz + 1e-8 * sc * Matrix::Diagonal(q), zy))
  }
  out
}

S7::method(term_build, SegTerm) <- function(term, data, ...) {
  xv <- .seg_x(term@var, data)
  n <- length(xv)
  params <- .seg_params(term@kind, term@npsi, term@linear)
  dev <- .seg_designs(term, data)

  rr <- range(xv)
  if (diff(rr) <= 0) {
    stop("the covariate of a break-point term must vary.", call. = FALSE)
  }

  # the interval a break-point is held in: far enough inside the data
  # that both sides carry observations
  lim <- as.numeric(stats::quantile(xv, c(0.05, 0.95), names = FALSE))

  # The convergence tolerance of Fasola et al., a hundredth of the
  # distance between consecutive distinct observations: below it the
  # objective, a step function, cannot change. They write the SMALLEST
  # such distance, which is the same number on the evenly spaced
  # covariates of their examples and is of order n^-2 on a random one,
  # so the median is taken instead and the two rules agree wherever
  # theirs is usable.
  ux <- sort(unique(xv))
  delta <- 0.01 * if (length(ux) > 1L) stats::median(diff(ux)) else diff(rr)

  start_psi <- if (!is.null(term@spec$psi)) term@spec$psi else {
    as.numeric(stats::quantile(xv, seq_len(term@npsi) / (term@npsi + 1),
                               names = FALSE))
  }

  # the layout: one stretch of coefficients per own parameter, in the
  # order the block carries them
  Z <- stats::setNames(vector("list", length(params)), params)
  index <- stats::setNames(vector("list", length(params)), params)
  cn <- character(0)
  pos <- 0L
  for (p in params) {
    stem <- .seg_coef_stem(term@kind, p)
    d <- dev[[p]]
    q <- if (is.null(d)) 1L else ncol(d$Z)
    Z[[p]] <- if (is.null(d)) NULL else d$Z
    index[[p]] <- pos + seq_len(q)
    pos <- pos + q
    cn <- c(cn, if (is.null(d)) paste(term@label, stem, sep = ".")
            else paste(term@label, stem, d$coef_names, sep = "."))
  }

  # the starting coefficients: unit changes, the break-points at their
  # starting positions, and g = -delta * psi where a step reads them off.
  # A development starts at the vector projecting the constant start onto
  # its own design.
  coef0 <- numeric(pos)
  start_at <- function(p, v) {
    d <- dev[[p]]
    coef0[index[[p]]] <<- if (is.null(d)) v else .seg_proj_start(d$Z, v)
  }
  for (k in seq_len(term@npsi)) {
    if (term@kind %in% c("seg", "jseg")) start_at(paste0("gamma", k), 1)
    if (term@kind %in% c("jump", "jseg")) start_at(paste0("delta", k), 1)
    start_at(paste0("psi", k),
             if (term@kind == "seg") start_psi[k] else -start_psi[k])
  }

  bp <- list(kind = term@kind, npsi = term@npsi, linear = term@linear,
             params = params, Z = Z, index = index,
             subs = lapply(dev, function(d) if (is.null(d)) NULL else d$subs),
             developed = any(!vapply(dev, is.null, logical(1))),
             lo = rr[1L], hi = rr[2L], lim = lim, delta = delta,
             cscale = rep(term@spec$c0, term@npsi),
             sgn = rep(0, term@npsi), step = rep(NA_real_, term@npsi),
             nref = 0L, var = term@var, coef = coef0, xv = xv,
             pk = vector("list", term@npsi))

  asm <- .seg_assemble(bp, xv, coef0)
  X <- asm$X
  colnames(X) <- cn

  # The penalties a development declares, keyed by the coefficient it
  # develops. Where one subformula covers every coefficient of a kind the
  # entries are POOLED into one, so that `gamma ~ 0 + lasso(~1)` is a
  # penalized block over the K changes under one hyperparameter, which is
  # what selecting how many break-points there are asks for; per-parameter
  # subformulas name one coefficient each and there is nothing to pool.
  bp$penalties <- .seg_penalties(term, dev, index)

  bp$value <- asm$value
  bp$psi <- asm$psi
  bp$pk <- asm$pk
  term@X <- X
  term@coef_names <- cn
  term@blueprint <- bp
  term@penalty <- NULL
  term
}

# The penalty entries of the developments. A sub-term shared by several
# coefficients of one kind declares one penalty per coefficient, and K
# hyperparameters over the same quantity is not what a shared subformula
# says: the entries are pooled, the penalty rebuilt over the union through
# the sub-term's own factory. A sub-term that has no factory -- anything
# but the penalized constructors -- keeps one entry per coefficient, its
# penalty being sized to its own block.
.seg_penalties <- function(term, dev, index) {
  out <- list()
  devd <- names(dev)[!vapply(dev, is.null, logical(1))]
  kinds <- sub("[0-9]+$", "", devd)
  pooled <- character(0)
  for (g in unique(kinds)) {
    ps <- devd[kinds == g]
    keys <- vapply(ps, function(p) dev[[p]]$key, character(1))
    ent <- dev[[ps[1L]]]$pens
    if (length(ps) < 2L || length(unique(keys)) > 1L || !length(ent)) next
    fac <- lapply(ent, function(e) .seg_sub_factory(dev[[ps[1L]]]$subs, e$name))
    if (any(vapply(fac, is.null, logical(1)))) next
    for (i in seq_along(ent)) {
      idx <- unlist(lapply(ps, function(p) index[[p]][ent[[i]]$index]),
                    use.names = FALSE)
      out[[length(out) + 1L]] <- list(
        name = paste0(g, "::", ent[[i]]$name), index = idx,
        penalty = fac[[i]](length(idx)),
        # the pooled entry is ONE penalty over the union, so it holds what
        # the shared sub-term holds -- they are the same declaration read
        # once per coefficient of the kind
        fixed = ent[[i]]$fixed, n_values = ent[[i]]$n_values,
        min_ratio = ent[[i]]$min_ratio)
    }
    pooled <- c(pooled, ps)
  }
  for (p in setdiff(devd, pooled)) {
    for (e in dev[[p]]$pens) {
      out[[length(out) + 1L]] <- list(name = paste0(p, "::", e$name),
                                      index = index[[p]][e$index],
                                      penalty = e$penalty, fixed = e$fixed,
                                      n_values = e$n_values,
                                      min_ratio = e$min_ratio,
                                      search = e$search)
    }
  }
  out
}

# The factory of the sub-term an entry came from, where it has one. The
# property is asked for rather than the class tested, which is what lets a
# term written later be pooled without an edit here.
.seg_sub_factory <- function(subs, name) {
  lb <- sub("::.*$", "", name)
  tm <- subs[[lb]]
  if (is.null(tm) || !("factory" %in% S7::prop_names(tm))) return(NULL)
  tm@factory
}

#' @title Penalties of a Break-Point Term
#' @name term_penalties.SegTerm
#' @description
#' One entry per penalty the developments of the term's own coefficients
#' declare, naming the coefficients it covers. The entries of a subformula
#' shared by every coefficient of one kind are pooled into one, so that
#' \code{gamma ~ 0 + lasso(~1)} is a single penalized block over the changes
#' under one hyperparameter. The list is empty for a term whose
#' coefficients carry no development, and for a specification, whose
#' coefficients there is nothing yet to index: a penalty is attached at
#' build, as it is for every penalized term here.
#' @param term A built \code{\link{SegTerm}}.
#' @param ... Unused.
#' @return A list of entries, as \code{\link{term_penalties}} documents.
#' @keywords internal
S7::method(term_penalties, SegTerm) <- function(term, ...) {
  pens <- term@blueprint$penalties
  if (is.null(pens)) list() else pens
}

S7::method(term_block_contract, SegTerm) <- function(term, coef = NULL, A,
                                                     ...) {
  bp <- term@blueprint
  if (!length(bp)) stop("the term is not built.", call. = FALSE)
  cf <- if (is.null(coef)) bp$coef else as.numeric(coef)
  A <- as.matrix(A)
  ncoef <- length(term_coef_names(term))
  n <- length(bp$xv)
  if (nrow(A) != n || ncol(A) != ncoef) {
    stop(sprintf(paste0("'A' must be %d by %d, the observations by the term's",
                        " coefficients."), n, ncoef), call. = FALSE)
  }
  out <- numeric(ncoef)
  # ⚠️ THE CONTINUOUS CONSTRUCTION ONLY. In a jump or a jseg the position is
  # READ OFF a product of the unknowns, psi = -g/delta, so the derivative of a
  # column in the coefficients runs through that read-off rather than through
  # a development's design, and the weight W = 1/(2|x~ - psi|) those two carry
  # has a derivative in the break-point that is unbounded -- measured at h,
  # h/4 and h/16 the quotient reads 3.6e4, 1.4e5 and 5.8e5. Their block is a
  # working LINEARIZATION with a frozen weight rather than a Jacobian, which
  # is the same fact that makes term_converged() answer differently for them.
  # Zeros here are what the base class gives, and are honest: a partial answer
  # that looked complete would be worse. Measured, this is the whole of the
  # gap on a jseg -- the contraction over its change columns reads 0 against a
  # brute-force 1.30e+02.
  if (!identical(bp$kind, "seg")) return(out)
  # every column of the block is a multiplier times that coefficient's own
  # design, so the contraction over a parameter's columns collapses to one
  # weight per observation before any derivative is taken
  sfun <- function(p) {
    idx <- bp$index[[p]]
    z <- bp$Z[[p]]
    Ap <- A[, idx, drop = FALSE]
    if (is.null(z)) as.numeric(Ap[, 1L]) else
      as.numeric(rowSums(Ap * as.matrix(z)))
  }
  place <- function(p, w) {
    idx <- bp$index[[p]]
    z <- bp$Z[[p]]
    out[idx] <<- out[idx] +
      if (is.null(z)) sum(w) else as.numeric(crossprod(as.matrix(z), w))
  }

  xv <- bp$xv
  pos <- .seg_positions(bp, cf, n)
  for (k in seq_len(bp$npsi)) {
    pk <- paste0("psi", k)
    gk <- paste0("gamma", k)
    psi <- pos$psi[, k]
    ind <- as.numeric(xv > psi)
    # where the position is against a limit it does not move with its
    # coefficients, and neither does anything downstream of it
    zp <- bp$Z[[pk]]
    v <- if (is.null(zp)) rep(pos$pk[[k]], n) else
      as.numeric(as.matrix(zp) %*% pos$pk[[k]])
    free <- as.numeric(v > bp$lim[1L] & v < bp$lim[2L])

    # the truncated line (x - psi)_+ has derivative -1(x > psi) in the
    # break-point: bounded, and existing everywhere but at the point itself
    place(pk, -sfun(gk) * ind * free)
    # the break-point column is -gamma(x) 1(x > psi). Its derivative in the
    # CHANGE is the same indicator; in the break-point it is zero almost
    # everywhere, the indicator being a step, and that is the value taken.
    place(gk, -sfun(pk) * ind)
  }
  out
}


S7::method(term_refresh, SegTerm) <- function(term, coef, ...) {
  .assert_built(term)
  bp <- term@blueprint
  coef <- as.numeric(coef)
  if (length(coef) != ncol(term@X)) {
    stop(sprintf("'coef' must have length %d.", ncol(term@X)), call. = FALSE)
  }

  # The scaling schedule of Fasola et al.: the factor is halved whenever
  # the break-point reverses direction, which is the signal that the
  # iteration has begun to circle an optimum rather than travel towards
  # one. A large factor moves the estimate far and lets it leave a
  # spurious optimum; a small one is faithful to the step function.
  # A developed break-point moves once per observation, so the schedule is
  # driven by the mean movement for the direction and by the largest for
  # the step, and the conditioning bound is taken at the observation
  # nearest an end of the range.
  psi_new <- .seg_positions(bp, coef, length(bp$xv))$psi
  dm <- psi_new - bp$psi
  d <- colMeans(dm)
  stepv <- apply(abs(dm), 2L, max)
  dfar <- apply(psi_new, 2L, function(p) max(pmax(p - bp$lo, bp$hi - p)))
  dnear <- apply(psi_new, 2L, function(p) min(pmin(p - bp$lo, bp$hi - p)))
  amax <- apply(abs(psi_new), 2L, max)

  s <- sign(d)
  flip <- s != 0 & bp$sgn != 0 & s != bp$sgn
  bp$cscale[flip] <- bp$cscale[flip] / 2
  # The factor is not allowed to collapse. Writing D and d for the
  # distances from the break-point to the further and the nearer end of
  # the range, the weight spans D/(c d) and Z = psi W + 1(x > psi) shears
  # the pair by a further |psi|/D, so the condition number of the block
  # is of order max(D, |psi|)/(c d). Holding that below eps^-1/2, which
  # leaves half the digits of a double to a QR of the design, is
  #     c >= sqrt(eps) * max(D, |psi|) / d.
  # The bound is on the DESIGN: a caller forming the normal equations
  # squares it and needs a factor a thousand times larger, which is why
  # the working model is fitted by a QR of X, as `segmented` does.
  cmin <- sqrt(.Machine$double.eps) * pmax(dfar, amax) / dnear
  bp$cscale <- pmax(bp$cscale, cmin)
  bp$sgn[s != 0] <- s[s != 0]
  # The first refresh evaluates the block at the starting coefficients,
  # so the break-point has not moved and the difference is zero by
  # construction rather than because the iteration has finished.
  bp$nref <- bp$nref + 1L
  bp$step <- if (bp$nref > 1L) stepv else rep(NA_real_, length(stepv))

  asm <- .seg_assemble(bp, bp$xv, coef)
  X <- asm$X
  colnames(X) <- term@coef_names
  bp$coef <- coef
  bp$value <- asm$value
  bp$psi <- asm$psi
  bp$pk <- asm$pk
  term@X <- X
  term@blueprint <- bp
  term
}

# The developments' designs on other rows, each sub-term reapplied through
# its own blueprint
.seg_z_at <- function(bp, newdata) {
  if (!bp$developed) return(bp$Z)
  lapply(stats::setNames(bp$params, bp$params), function(p) {
    if (is.null(bp$subs[[p]])) NULL
    else .nl_bind(lapply(bp$subs[[p]], term_predict, newdata = newdata))
  })
}

S7::method(term_value, SegTerm) <- function(term, coef = NULL, newdata = NULL,
                                            ...) {
  .assert_built(term)
  bp <- term@blueprint
  if (!is.null(newdata)) {
    # the break-points are the ones the term carries, as term_predict()
    # reapplies them: refreshing here would read them off the new rows
    xv <- .seg_x(bp$var, newdata)
    cf <- if (is.null(coef)) bp$coef else as.numeric(coef)
    return(.seg_assemble(bp, xv, cf, bp$cscale, .seg_z_at(bp, newdata),
                         bp$pk)$value)
  }
  if (is.null(coef)) return(bp$value)
  term_refresh(term, coef)@blueprint$value
}

S7::method(term_predict, SegTerm) <- function(term, newdata, ...) {
  .assert_built(term)
  bp <- term@blueprint
  xv <- .seg_x(bp$var, newdata)
  X <- .seg_assemble(bp, xv, bp$coef, bp$cscale, .seg_z_at(bp, newdata),
                     bp$pk)$X
  colnames(X) <- term@coef_names
  X
}
#' The Break-Points of a Break-Point Term
#'
#' @description
#' The estimated positions of the break-points. For a continuous term
#' these are coefficients; for a discontinuous one they are read off two
#' coefficients as \eqn{-g_k/\delta_k}, so this is the function that
#' reports them either way.
#'
#' @param term A built \code{\link{SegTerm}}.
#' @param coef Optional coefficients; defaults to the ones the term was
#'   last refreshed at.
#'
#' @return A numeric vector with one entry per break-point, or, where a
#'   break-point carries a development, a matrix with one row per
#'   observation and one column per break-point.
#'
#' @examples
#' set.seed(1)
#' dd <- data.frame(x = sort(runif(100, 0, 10)))
#' seg_psi(term_build(seg(x, npsi = 2), dd))
#'
#' @seealso \code{\link{seg}}, \code{\link{jump}}, \code{\link{jseg}}
#' @export
seg_psi <- function(term, coef = NULL) {
  if (!S7::S7_inherits(term, SegTerm)) {
    stop("'term' must be a break-point term.", call. = FALSE)
  }
  .assert_built(term)
  tm <- if (is.null(coef)) term else term_refresh(term, coef)
  .seg_psi_report(tm@blueprint)
}

# The positions as a reader wants them: one number per break-point where
# none of them is developed, since the matrix then repeats one row.
.seg_psi_report <- function(bp) {
  psi <- bp$psi
  devd <- any(vapply(paste0("psi", seq_len(bp$npsi)),
                     function(p) !is.null(bp$Z[[p]]), logical(1)))
  if (devd) psi else as.numeric(psi[1L, ])
}

#' The Progress of a Break-Point Iteration
#'
#' @description
#' \code{seg_step} returns how far each break-point moved at the last
#' call to \code{\link{term_refresh}}, and \code{seg_converged} compares
#' the largest of those with the tolerance of \cite{fasola2018}, a
#' hundredth of the distance between consecutive distinct observations of
#' the covariate. A term that has been built but not yet refreshed has
#' taken no step, so \code{seg_step} returns \code{NA} and
#' \code{seg_converged} returns \code{FALSE}.
#'
#' @details
#' With \eqn{x_{(1)} < \cdots < x_{(m)}} the distinct covariate values, the
#' run stops at
#'
#' \deqn{\max_k \lvert \psi_k^{(t)} - \psi_k^{(t-1)} \rvert < \Delta,
#'   \qquad \Delta = 0.01 \cdot
#'     \operatorname{median}_{i}\, (x_{(i+1)} - x_{(i)}).}
#'
#' \cite{fasola2018} take the smallest of those gaps rather than their
#' median, which agrees with this on the evenly spaced covariates of their
#' examples and is of order \eqn{m^{-2}} on a random one, hence unreachable.
#'
#' The rule is one of resolution: below that distance the objective of a
#' discontinuous term, a step function of the break-point, cannot change.
#' It therefore tightens as the sample grows while the precision the
#' fixed point is reached at does not, so on a large sample the last
#' iterations move the break-point by a little more than the rule allows
#' and the run continues past the point where the estimate has settled.
#' A caller that can evaluate the objective should stop on its relative
#' change instead, which is what \code{segmented} does and what costs a
#' continuous term nothing: its iteration can settle into a cycle of
#' period two in the break-point, in which case this rule is never met
#' while the objective has long since stopped moving.
#'
#' @param term A built \code{\link{SegTerm}}.
#'
#' @return \code{seg_step} returns a numeric vector with one entry per
#'   break-point; \code{seg_converged} returns a single logical.
#'
#' @examples
#' set.seed(1)
#' dd <- data.frame(x = sort(runif(100, 0, 10)))
#' dd$y <- 2 * (dd$x > 6) + rnorm(100, sd = 0.3)
#' b <- term_build(jump(x, psi = 4), dd)
#' cf <- b@blueprint$coef
#' for (it in 1:30) {
#'   b <- term_refresh(b, cf)
#'   X <- term_matrix(b)
#'   cf <- as.numeric(qr.solve(crossprod(X), crossprod(X, dd$y)))
#'   if (seg_converged(b)) break
#' }
#' c(psi = seg_psi(b, cf), step = seg_step(b))
#'
#' @seealso \code{\link{seg}}, \code{\link{seg_psi}}, \code{\link{seg_start}}
#'   \code{\link{seg_start}}
#' @export
seg_step <- function(term) {
  if (!S7::S7_inherits(term, SegTerm)) {
    stop("'term' must be a break-point term.", call. = FALSE)
  }
  .assert_built(term)
  term@blueprint$step
}

#' @rdname seg_step
#' @export
seg_converged <- function(term) {
  st <- seg_step(term)
  !anyNA(st) && max(st) < term@blueprint$delta
}

#' @title Whether a Break-Point Term's Break-Points Have Settled
#' @name term_converged.SegTerm
#' @description
#' \code{\link{seg_converged}}, so that a fitting layer reads the rule the
#' construction is stopped on without knowing it is holding a break-point
#' term. A term that has not been refreshed has not moved and reports
#' \code{FALSE}, the first refresh being the one that measures nothing.
#' @param term A built \code{\link{SegTerm}}.
#' @param ... Unused.
#' @return A single logical.
#' @keywords internal
S7::method(term_converged, SegTerm) <- function(term, ...) {
  isTRUE(seg_converged(term))
}

#' Starting Positions for a Break-Point Term
#'
#' @description
#' Chooses the starting positions of a \code{\link{seg}}, \code{\link{jump}}
#' or \code{\link{jseg}} term by scoring an equally spaced grid on the
#' least-squares profile of the term's own columns, and returns the
#' specification with \code{psi} set to the best combination found.
#'
#' @details
#' \cite{fasola2018} recommend fixing the starting value by evaluating
#' the objective on a small grid spanned over the range of the covariate
#' rather than at a single conventional point, and the recommendation
#' matters more than it sounds: the objective has local optima in the
#' break-point, and the iteration converges from within a basin around
#' the position it starts at. Measured on a joint jump and change of
#' slope in 500 observations, over eight samples, the fraction of runs
#' recovering the break-point is 0 to 0.5 depending on where a single
#' start is placed and 1 from the grid.
#'
#' Writing \eqn{X(\psi)} for the design the term produces at a candidate
#' position, the position chosen is
#'
#' \deqn{\hat\psi = \arg\min_{\psi \in \mathcal{G}}
#'   \bigl\lVert y - X(\psi)\,
#'     \widehat{\beta}(\psi) \bigr\rVert^{2},
#'   \qquad \widehat{\beta}(\psi)
#'     = \arg\min_{\beta} \lVert y - X(\psi)\beta \rVert^{2},}
#'
#' over an equally spaced grid \eqn{\mathcal{G}} of \code{k} points in the
#' range of the covariate. The inner minimization is a linear fit, so the
#' whole rule costs \code{k} of them.
#'
#' The grid is scored on the residual sum of squares of an intercept,
#' the term's columns at each candidate position and, where the term
#' carries one, the linear effect. That is the exact profile for a
#' gaussian response and an adequate starting rule for any other, the
#' quantity being used to place a starting value and not to fit. The
#' positions found seed a development as well, each starting vector
#' projecting the position onto the sub-design. With several break-points
#' every increasing combination of grid points is scored, so \code{k}
#' should be kept small.
#'
#' @param spec An unbuilt \code{\link{SegTerm}}.
#' @param data A data frame in which the covariate is evaluated.
#' @param y The response, one value per row of \code{data}.
#' @param k The number of grid points. Defaults to 10.
#'
#' @return The specification, with \code{psi} set.
#'
#' @references
#' Fasola, S., Muggeo, V. M. R. and Kuchenhoff, H. (2018). A heuristic,
#' iterative algorithm for change-point detection in abrupt change
#' models. \emph{Computational Statistics}, 33, 997--1015.
#'
#' @examples
#' set.seed(1)
#' dd <- data.frame(x = sort(runif(200, 0, 10)))
#' dd$y <- 0.3 * dd$x + 2 * (dd$x > 6.5) + rnorm(200, sd = 0.3)
#' seg_start(jump(x), dd, dd$y)@spec$psi
#'
#' @seealso \code{\link{seg}}, \code{\link{jump}}, \code{\link{jseg}}
#' @export
seg_start <- function(spec, data, y, k = 10) {
  if (!S7::S7_inherits(spec, SegTerm)) {
    stop("'spec' must be a break-point term.", call. = FALSE)
  }
  if (length(spec@X) > 0L) {
    stop("'spec' must be an unbuilt term; see term_build().", call. = FALSE)
  }
  if (!is.numeric(k) || length(k) != 1L || is.na(k) || k < spec@npsi + 1) {
    stop(sprintf("'k' must be a single number of at least %d.",
                 spec@npsi + 1L), call. = FALSE)
  }
  xv <- .seg_x(spec@var, data)
  y <- as.numeric(y)
  if (length(y) != length(xv)) {
    stop("'y' must have one value per row of 'data'.", call. = FALSE)
  }

  lim <- as.numeric(stats::quantile(xv, c(0.05, 0.95), names = FALSE))
  grid <- seq(lim[1L], lim[2L], length.out = as.integer(k))
  combos <- utils::combn(length(grid), spec@npsi, simplify = FALSE)

  cols <- function(x, psi) {
    Z <- if (spec@linear) cbind(1, x) else matrix(1, length(x), 1L)
    if (spec@kind %in% c("seg", "jseg")) {
      for (p in psi) Z <- cbind(Z, pmax(x - p, 0))
    }
    if (spec@kind %in% c("jump", "jseg")) {
      for (p in psi) Z <- cbind(Z, as.numeric(x > p))
    }
    Z
  }

  v <- vapply(combos, function(ii) {
    Z <- cols(xv, grid[ii])
    out <- tryCatch(sum(qr.resid(qr(Z), y)^2), error = function(e) Inf)
    if (is.finite(out)) out else Inf
  }, numeric(1))
  best <- grid[combos[[which.min(v)]]]

  spec@spec$psi <- best
  spec
}

S7::method(print, SegTerm) <- function(x, ...) {
  if (term_is_built(x)) {
    dv <- names(x@subformulas)
    cat(sprintf("<SegTerm> '%s': %s, %d break-point%s%s
", x@label, x@kind,
                x@npsi, if (x@npsi == 1L) "" else "s",
                if (length(dv))
                  sprintf("; %s developed", paste(dv, collapse = ", ")) else ""))
    psi <- .seg_psi_report(x@blueprint)
    if (is.matrix(psi)) {
      # a developed break-point has one position per observation; what a
      # reader wants is where each one ranges
      rng <- apply(psi, 2L, function(p)
        sprintf("[%s, %s]", format(min(p), digits = 4),
                format(max(p), digits = 4)))
      cat("  at (developed): ", paste(rng, collapse = ", "), "
", sep = "")
    } else {
      cat("  at: ", paste(format(psi, digits = 4), collapse = ", "),
          "
", sep = "")
    }
  } else {
    cat(sprintf("<SegTerm> '%s': %s (specification)
", x@label, x@kind))
  }
  invisible(x)
}

#' @title Where a Break-Point Term's Coefficients Begin
#' @name term_coef_start.SegTerm
#' @description
#' The start \code{\link{term_build}} computed: unit changes, and the
#' break-points at the positions \code{psi} names or at the interior
#' quantiles of the covariate. Zero is degenerate here rather than neutral
#' -- a discontinuous term reads its break-point off \eqn{-g_k/\delta_k},
#' which at zero is the same clamped position for every one of them, and a
#' continuous term's Jacobian column vanishes -- so a fitting layer that
#' starts every coefficient at zero has to be told otherwise.
#' @param term A built \code{\link{SegTerm}}.
#' @param ... Unused.
#' @return A numeric vector, one value per column of the block.
#' @keywords internal
S7::method(term_coef_start, SegTerm) <- function(term, ...) {
  .assert_built(term)
  term@blueprint$coef
}

#' @title What a Fitted Break-Point Term Is About
#' @name term_readable.SegTerm
#'
#' @description
#' The quantities of the model the term defines, rather than the
#' coefficients of the working block it is fitted through: the linear
#' effect \eqn{\beta} where the term carries one, the changes of slope
#' \eqn{\gamma_k} and of level \eqn{\delta_k}, and the break-points
#' \eqn{\psi_k}, with the Jacobian from the coefficients so that a caller
#' holding their variance matrix can carry it across.
#'
#' @details
#' A continuous term holds every one of those as a coefficient and the
#' Jacobian is a selection. A discontinuous one does not hold the
#' break-point at all: it is read off the auxiliary pair as \eqn{\psi_k =
#' -g_k/\delta_k}, which the joint construction reaches too, its quadratic
#' degenerating to that reading once the increment has vanished. The
#' Jacobian rows are then
#'
#' \deqn{\frac{\partial \psi_k}{\partial g_k} = -\frac{1}{\delta_k},
#'   \qquad
#'   \frac{\partial \psi_k}{\partial \delta_k} = \frac{g_k}{\delta_k^{2}},}
#'
#' which is the delta method \code{segmented} reports a break-point's
#' standard error by. The joint construction reads its position from a
#' quadratic that also carries the change of slope, and the two readings
#' differ by the increment the iteration has left: the VALUE reported is
#' the one the term holds, so that this and \code{\link{seg_psi}} cannot
#' give two numbers for one quantity, and the Jacobian is the fixed
#' point's, which the quadratic degenerates to once the increment has
#' vanished.
#'
#' Every quantity is on the identity scale: a change is unbounded and a
#' break-point is a position on the covariate's own scale, held inside the
#' interval between the 5th and the 95th percentile. Where a coefficient
#' carries a development there is no single number to report -- a
#' break-point then has one value per observation -- and the method
#' returns nothing, leaving a caller to report the coefficients
#' themselves.
#'
#' @param term A built \code{\link{SegTerm}}.
#' @param zeta The term's coefficients, in the order of its block.
#' @param ... Unused.
#'
#' @return A list with \code{name}, \code{value}, \code{jacobian} and
#'   \code{scale}, as \code{\link{term_readable}} documents, or
#'   \code{NULL} where a coefficient carries a development.
#'
#' @keywords internal
S7::method(term_readable, SegTerm) <- function(term, zeta, ...) {
  .assert_built(term)
  bp <- term@blueprint
  if (isTRUE(bp$developed)) return(NULL)
  cf <- as.numeric(zeta)
  K <- bp$npsi
  np <- length(cf)
  nm <- character(0)
  val <- numeric(0)
  rows <- list()
  at <- function(p) bp$index[[p]][1L]
  push <- function(label, v, j) {
    nm <<- c(nm, label)
    val <<- c(val, v)
    rows[[length(rows) + 1L]] <<- j
  }
  sel <- function(i) { e <- numeric(np); e[i] <- 1; e }

  if (bp$linear) push("beta", cf[at("beta")], sel(at("beta")))
  if (bp$kind %in% c("seg", "jseg")) {
    for (k in seq_len(K)) {
      i <- at(paste0("gamma", k))
      push(paste0("gamma", k), cf[i], sel(i))
    }
  }
  if (bp$kind %in% c("jump", "jseg")) {
    for (k in seq_len(K)) {
      i <- at(paste0("delta", k))
      push(paste0("delta", k), cf[i], sel(i))
    }
  }
  # the positions the TERM holds, read by the same function it reads them
  # with, so that this and seg_psi() cannot report two numbers for one
  # quantity
  pos <- .seg_positions(bp, cf, 1L, bp$Z, bp$pk)$psi[1L, ]
  for (k in seq_len(K)) {
    ip <- at(paste0("psi", k))
    if (bp$kind == "seg") {
      push(paste0("psi", k), pos[k], sel(ip))
    } else {
      id <- at(paste0("delta", k))
      d <- cf[id]
      if (abs(d) < 1e-12) d <- sign(d + (d == 0)) * 1e-12
      j <- numeric(np)
      j[ip] <- -1 / d
      j[id] <- cf[ip] / d^2
      push(paste0("psi", k), pos[k], j)
    }
  }
  J <- do.call(rbind, rows)
  dimnames(J) <- list(nm, term@coef_names)
  list(name = nm, value = val, jacobian = J,
       scale = stats::setNames(
         rep(list(linkfunctions7::identity_link()), length(nm)), nm))
}
