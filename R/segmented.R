#' @include term_classes.R generics.R nonlinear.R
NULL

#' @title S7 Class for Break-Point Terms
#' @name SegTerm
#'
#' @description
#' The subclass of [additive_term()] for a covariate whose effect changes at
#' estimated break-points: a change of slope ([seg()]), a change of level
#' ([jump()]), or both at the same points ([jseg()]). Its design block is the
#' **working** block of the iteration that estimates the break-points, and
#' [term_refresh()] recomputes it as they move.
#'
#' @details
#' # The six properties of its own
#'
#' `kind` is `"seg"`, `"jump"` or `"jseg"`, and it decides almost everything
#' else: which coefficients the term has, whether its block is a Jacobian, and
#' which developments it will accept.
#'
#' `var` is the covariate expression, `npsi` the number of break-points, and
#' `linear` whether the block carries the linear effect \eqn{\beta x}, which
#' `seg` and `jseg` do by default and `jump` never does.
#'
#' `subformulas` holds one right-hand side per developed coefficient, and
#' `spec` the resolved construction settings: the starting positions, the
#' scaling factor `c0`, the restart budget `n_boot` and any smoother.
#'
#' # Two kinds of block, and what distinguishes them
#'
#' For [seg()] the block is the exact Jacobian: the contribution is
#' differentiable in \eqn{\psi_k} away from the break-point, so the
#' break-point is an ordinary coefficient and a linear fit on the block is a
#' Gauss-Newton step. [term_jacobian_block()] answers `TRUE`.
#'
#' For [jump()] and [jseg()] it is a working linearization with the weight
#' \eqn{W = 1/(2\lvert\tilde{x} - \psi\rvert)} frozen at the previous
#' break-point, and the break-point is **read off** two coefficients rather
#' than incremented. [term_jacobian_block()] answers `FALSE`, and a fitting
#' layer routes such a term to exact working fits alternating with committed
#' read-offs. Smoothing the step ([seg()]'s `smoothed` argument) turns the
#' break-point back into an ordinary parameter, and the answer back to `TRUE`.
#'
#' [seg_psi()] is the one function that reports the position under either
#' construction.
#'
#' @inheritParams additive_term
#' @param kind One of `"seg"`, `"jump"` or `"jseg"`.
#' @param var The covariate expression, unevaluated.
#' @param npsi The number of break-points, an integer of at least 1.
#' @param linear Whether the block carries the linear effect.
#' @param subformulas A named list of one-sided formulas, one per developed
#'   coefficient. Empty where none is.
#' @param spec A named list of the resolved construction settings.
#'
#' @return An S7 object of class `SegTerm`, inheriting from [additive_term()]
#'   and [model_term()], with the six properties above beside the ten they
#'   supply.
#'
#' @seealso [seg()], [jump()] and [jseg()] for the three constructions;
#'   [seg_psi()] for the positions; [term_refresh()] for the block as they
#'   move; [seg_start()] for where to begin.
#'
#' @examples
#' set.seed(1)
#' d <- data.frame(x = sort(runif(120, 0, 10)))
#' d$y <- 1 + 0.5 * d$x + 2 * pmax(d$x - 6, 0) + rnorm(120, sd = 0.4)
#'
#' # All three constructions are this class, and differ in `kind`.
#' vapply(list(seg(x), jump(x), jseg(x)), function(t) t@kind, character(1))
#'
#' # Which decides the coefficients the term carries.
#' lapply(list(seg = seg(x), jump = jump(x), jseg = jseg(x)),
#'        function(t) term_coef_names(term_build(t, d)))
#'
#' # And whether its block is a Jacobian or a working linearization.
#' vapply(list(seg = seg(x), jump = jump(x), jseg = jseg(x)),
#'        term_jacobian_block, logical(1))
#'
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
#' mean of a `gaussian1_distrib`; nothing here is particular to it.
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
#' coefficients ([term_refresh()]) and fitting is one iteration
#' of it, which is the same contract [nl()] uses.
#'
#' The objective has local optima in the break-points and the iteration
#' converges from within a basin around where it starts, so
#' [seg_start()] chooses the starting positions on a grid rather
#' than at a conventional point. A break-point is confined to the interval
#' between the 5th and the 95th percentile of the covariate: outside it
#' the indicator is constant, the truncated line and that constant are
#' linearly dependent, and the block is exactly singular, not merely
#' ill-conditioned. A run ending against the limit has not located a
#' break-point, and [seg_psi()] then reports the limit itself.
#'
#' The iteration can settle into a cycle of period two in the break-point,
#' in which case the rule of [seg_converged()] is never met while
#' the objective has long since stopped moving; a caller that can evaluate
#' the objective should stop on its relative change, as `segmented`
#' does.
#' }
#'
#' @section The coefficients:
#' `beta` (present when `linear = TRUE`), `gamma1` \ldots
#' `gammaK` and `psi1` \ldots `psiK`, in that order. The
#' break-points are coefficients of the working block here, which they are
#' not in the discontinuous constructions, so [seg_psi()] is the
#' function that reports them either way.
#'
#' @section Developing a coefficient with covariates:
#' Every coefficient of the term is a parameter like any other and may be
#' developed on covariates, through a two-sided formula in `...` whose
#' left side names it:
#'
#' \preformatted{seg(x, psi ~ id)                   a break-point per subject
#' seg(x, psi ~ random(~1 | id))      the random-changepoint model
#' seg(x, gamma1 ~ z)                 a first slope change varying with z
#' seg(x, by = ~0 + g)                every coefficient, per level of g}
#'
#' The right side goes through [interpret_formula()], so it takes
#' *any* term of this package, penalized ones included, and the
#' hyperparameters a sub-term declares are reported by
#' [term_penalties()] and estimated by the fitting layer. The
#' left side names either one coefficient (`gamma2`) or a whole kind
#' (`gamma`, meaning every \eqn{\gamma_k}, each with coefficients of
#' its own over the shared design). `by = ~f` is the shorthand giving
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
#' A penalty on the changes themselves, the lasso that selects how many
#' break-points are really there, is the development on an intercept
#' alone, `seg(x, npsi = 4, gamma ~ 0 + lasso(~1))`: the \eqn{K}
#' changes are then one penalized block under one hyperparameter, the
#' entries of a subformula shared by every coefficient of a kind being
#' pooled into one. The `0 +` is not decoration. A subformula follows
#' the intercept convention of a formula, so `gamma ~ lasso(~1)`
#' carries the subformula's own unpenalized intercept beside the penalized
#' column, which for an intercept-only development is the same column
#' twice.
#'
#' @section The linear effect:
#' With `linear = TRUE`, the default, the term carries \eqn{\beta x}
#' itself, so `y ~ seg(x)` is the whole relationship. Writing the
#' covariate again in the same equation is then exactly collinear with that
#' column, and [interpret_formula()] removes it with a warning.
#' `linear = FALSE` is the explicit way to keep the linear effect
#' outside the term, and `y ~ x + seg(x, linear = FALSE)` is the same
#' model written the other way round.
#'
#' @param x The covariate, an expression evaluated in the data.
#' @param ... Two-sided formulas developing the term's own coefficients;
#'   see the section above. Cannot be combined with `by`.
#' @param npsi The number of break-points. Defaults to 1.
#' @param psi Optional starting positions; defaults to evenly spaced
#'   quantiles of the covariate. Where a break-point carries a development
#'   they seed it, each starting vector solving \eqn{Wp_k \approx
#'   \psi_k^{0}}.
#' @param by A one-sided formula giving every coefficient of the term the
#'   same development, e.g. `by = ~0 + g` for an independent set per
#'   level of a factor. A bare variable is rejected; write the formula.
#' @param linear Whether the block carries the linear effect \eqn{\beta x}.
#'   Defaults to `TRUE`.
#' @param smoothed `NULL` (the default: the construction exactly as
#'   documented above) or a \pkg{penalties7}
#'   [penalties7::abs_smoother()], e.g.
#'   `penalties7::smooth_probit()`. The smoother replaces the step and
#'   the hinge by their smooth versions, \eqn{(1 + s'(u))/2} and
#'   \eqn{(u + s(u))/2}, so every break-point becomes an ordinary
#'   parameter of a \eqn{C^\infty} model: there is no working
#'   parametrization, no auxiliary coefficient and no scaling schedule
#'   (`c0` is ignored, with a message), the block is the true Jacobian
#'   and the term is fitted by Gauss-Newton like [nl()]. A
#'   development of a break-point, `psi ~ random(~1 | id)` and penalized
#'   ones included, is then legal for every kind, the read-off
#'   that constrained the discontinuous constructions having gone. The
#'   smoother's width is resolved at build from the covariate's spacing
#'   (the median gap between distinct values, within groups where a
#'   break-point development supplies a partition) unless the object
#'   carries one, and is reported: it is the width of the transition, the
#'   bent-cable reading, and the smoothing bias it buys is confined to a
#'   window of that width (probit, quintic) or decays as \eqn{c/(4|u|)}
#'   (hyperbolic). The objective is still multimodal in the positions --
#'   smoothing rounds the local optima and does not remove them, so the
#'   profile start and the `n_boot` restarts stay necessary; a
#'   smoothed fit from a bad start has been measured converging to an
#'   absurd local optimum while reporting success.
#' @param marginal Whether the break-point is a latent variable per group,
#'   integrated out of the likelihood exactly instead of being an estimated
#'   position. `FALSE`, the default, is the construction documented
#'   above. `TRUE` requires the subformula `psi ~ random(~1 | g)`
#'   and returns a structural term of the likelihood shape
#'   ([MarginalBreakTerm()]); see the section of
#'   [jump()], whose step model is where the marginal buys the
#'   most: for a `seg` term the native random-changepoint fit
#'   (`psi ~ random(~1 | g)` without `marginal`) is measured
#'   equivalent and remains the recommended route, the marginal being the
#'   exact-likelihood alternative. One break-point here: the conditional is
#'   smooth in the position, so the integral runs on a Gauss-Kronrod panel
#'   per interval between a group's ordered observations, and several
#'   latents would need a product quadrature whose component count no
#'   fitting layer can carry.
#' @param n_boot How many bootstrap restarts the fitting layer runs after
#'   the iteration first converges (Wood 2001, the device `segmented`
#'   runs by default): each restart re-estimates on a bootstrap resample
#'   from the current break-points and then on the data again from where
#'   the resample ended, keeping the better fit. The objective has local
#'   optima in the break-points, and with several of them a grid start does
#'   not reach the right basin from every placement. Defaults to 10,
#'   `segmented`'s own default; 0 disables. The term itself only
#'   declares the value; running the restarts is the fitting layer's,
#'   as with a penalty's hyperparameters.
#' @param label A single non-empty string prefixed to the coefficient
#'   names.
#'
#' @return An object of class [SegTerm()] (a specification; see
#'   [term_build()]).
#'
#' @references
#' Muggeo, V. M. R. (2003). Estimating regression models with unknown
#' break-points. *Statistics in Medicine*, 22(19), 3055--3071.
#'
#' Muggeo, V. M. R., Atkins, D. C., Gallop, R. J. and Dimidjian, S. (2014).
#' Segmented mixed models with random changepoints: a maximum likelihood
#' approach with application to treatment for depression study.
#' *Statistical Modelling*, 14(4), 293--313.
#'
#' Wood, S. N. (2001). Minimizing model fitting objectives that contain
#' spurious local minima by bootstrap restarting. *Biometrics*, 57(1),
#' 240--244.
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
#' @seealso [jump()], [jseg()], [seg_psi()],
#'   [seg_start()], [seg_step()], [nl()]
#' @export
seg <- function(x, ..., npsi = 1, psi = NULL, by = NULL, linear = TRUE,
                smoothed = NULL, marginal = FALSE, n_boot = 10,
                label = "seg") {
  .seg_check_marginal(marginal)
  if (marginal) {
    return(.marg_spec("seg", substitute(x), npsi, psi, substitute(by),
                      list(...), smoothed, FALSE, !missing(n_boot), label,
                      linear = linear))
  }
  # a continuous term has no scaling factor, its working block being
  # bounded already; the value is carried so that one blueprint serves the
  # three constructions
  .seg_spec("seg", substitute(x), npsi, psi, linear, 0.05, n_boot, label,
            list(...), .seg_by_formula(substitute(by), parent.frame()),
            smoothed)
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
#' level at \eqn{\psi_k}. A model that is linear in \eqn{x} *and*
#' steps is written by putting the covariate in the equation as well,
#' `y ~ x + jump(x)`, or by [jseg()] where the slope is to
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
#' rebuilding the weights, and why [seg_psi()] exists.
#' }
#'
#' \subsection{The scaling schedule}{
#' The weight is unbounded as \eqn{x} approaches \eqn{\psi}, and since
#' \eqn{Z - \psi W} is the step itself, \eqn{Z} is \eqn{\psi W} plus a
#' quantity of order one: an unbounded \eqn{W} makes the two columns
#' numerically collinear and drowns the signal the fit reads. The remedy of
#' \cite{fasola2018} is to move the observations instead of capping the
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
#' optimum and a small one is faithful to the step function. `c0` is
#' its starting value and [term_refresh()] halves it whenever the
#' break-point reverses direction, the signal that the iteration has begun
#' to circle an optimum instead of traveling toward one. The run has
#' converged when every break-point moves less than a hundredth of the
#' distance between consecutive distinct observations, which
#' [seg_converged()] reports.
#'
#' A break-point is confined to the interval between the 5th and the 95th
#' percentile of the covariate, outside which the indicator is constant and
#' the block singular, and the starting positions are chosen on a grid by
#' [seg_start()]: measured over eight samples and four starting
#' positions, the fraction of runs recovering the break-point is 0 to 0.5
#' from a single conventional start and 1 from the grid.
#' }
#'
#' @section The coefficients:
#' `delta1` \ldots `deltaK`, the changes of level, followed by
#' `g1` \ldots `gK`, the auxiliary coefficients of the identity
#' above, from which \eqn{\psi_k = -g_k/\delta_k} is read. The
#' break-points are *not* coefficients; ask [seg_psi()].
#'
#' @section Developing a coefficient with covariates:
#' As for [seg()], a two-sided formula in `...` whose left
#' side names `delta`, `deltaK` or `psi`, `psiK`
#' develops that coefficient on covariates, and `by = ~f` gives every
#' coefficient the same development. The right side goes through
#' [interpret_formula()] and takes any term of this package.
#'
#' What the discontinuous construction can carry is narrower than what the
#' continuous one can, and the reason is the read-off, never the
#' model. The auxiliary coefficient is \eqn{g_k = -\delta_k\psi_k}, a
#' product of the two quantities, so it stays linear in the unknowns only
#' when one of the two factors is a single number:
#'
#' - **the break-point developed**, the change of level a single
#'   number: exact for any design, \eqn{g_k = -\delta_k p_k} being
#'   linear in \eqn{p_k}. `jump(x, psi ~ id)` is a break-point
#'   per subject with a shared step size.
#' - **a development on group indicators**, where the design has
#'   one column per group and each observation belongs to one: the
#'   product collapses group by group and the read-off is exact within
#'   each. `jump(x, by = ~0 + g)` is an independent step and
#'   break-point per level.
#' - **anything else** is rejected, never approximated,
#'   the product of two developments needing the outer product of their
#'   designs and no unconstrained fit returning it as one.
#'
#' A sub-term carrying a penalty is rejected on a developed break-point,
#' since the estimated coefficients are \eqn{-\delta_k p_k}, the
#' development scaled by the step size, and a penalty would act on that
#' and never on the development. A penalty on the changes of level
#' themselves is `jump(x, npsi = 4, delta ~ 0 + lasso(~1))`, the
#' `0 +` removing the subformula's own unpenalized intercept, which
#' would otherwise be the same column twice.
#'
#' @section The marginal construction:
#' `jump(x, psi ~ random(~1 | g), marginal = TRUE)` treats each
#' break-point as a latent variable per group, \eqn{\psi_{ik} = m_k +
#' u_{ik}} with independent \eqn{u_{ik} \sim N(0, \tau_k^2)}, and
#' integrates them out of the likelihood exactly. The conditional is
#' constant on the product partition of the intervals between a group's
#' ordered covariate values, giving \eqn{(n_i+1)^K} cells, but the sum is
#' never taken over the cells: the side process \eqn{S_t = \{k : \psi_k
#' \le x_{(t)}\}} is monotone on the subset lattice with independent
#' coordinates, so it is a hidden Markov chain on the \eqn{2^K} side
#' patterns whose transition factors over the coordinates, and the exact
#' forward recursion costs \eqn{n K 2^K}: a coordinate that flips at a
#' step contributes its interval's prior mass as the transition weight,
#' one that never flips its upper-tail mass at the end. Up to eight
#' break-points are covered; what prices more is not the recursion but the
#' \eqn{2^K} components a fitting layer evaluates the family at.
#'
#' The prior is part of the likelihood: \eqn{(m_k, \tau_k, \delta_k)} are
#' ordinary parameters estimated by maximum likelihood, with no penalty, no
#' marginal criterion and no smoothing constant, and the term is structural
#' ([MarginalBreakTerm()], the likelihood shape of the contract)
#' itself. With one break-point the prior may be any
#' continuous \pkg{distributions7} distribution with its location fixed at
#' zero, through `random(distrib = )`: the interval masses are
#' differences of its cdf and their derivatives come from the cdf surface
#' built for the censored likelihoods, with that surface's own caveats
#' (closed for the gaussian and, in location and scale, for the t).
#' `smoothed` and `c0` do not apply and are ignored with a
#' message; the posterior of each group's break-points is read by
#' [term_latent()]. On the step model the marginal is the route
#' that resolves what a smoothed mode cannot: the conditional is a step
#' in the position, so a Laplace approximation has no curvature to read --
#' and it is the robust route on non-gaussian families.
#'
#' @inheritParams seg
#' @param marginal Whether the break-points are latent variables per group,
#'   integrated out of the likelihood exactly; see the section below.
#'   Requires the subformula `psi ~ random(~1 | g)`, whose grouping
#'   carries the latents. Defaults to `FALSE`, the construction
#'   documented above.
#' @param c0 The starting value of the scaling factor that separates the
#'   observations from the break-point, as a fraction of the distance to
#'   the ends of the range. Defaults to `0.05`, the value
#'   \cite{fasola2018} recommend; see Details.
#'
#' @return An object of class [SegTerm()] (a specification; see
#'   [term_build()]), or of class [MarginalBreakTerm()]
#'   with `marginal = TRUE`.
#'
#' @references
#' Fasola, S., Muggeo, V. M. R. and Kuchenhoff, H. (2018). A heuristic,
#' iterative algorithm for change-point detection in abrupt change
#' models. *Computational Statistics*, 33, 997--1015.
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
#' @seealso [seg()], [jseg()], [seg_psi()],
#'   [seg_start()], [seg_step()]
#' @export
jump <- function(x, ..., npsi = 1, psi = NULL, by = NULL, c0 = 0.05,
                 smoothed = NULL, marginal = FALSE, n_boot = 10,
                 label = "jump") {
  .seg_check_marginal(marginal)
  if (marginal) {
    return(.marg_spec("jump", substitute(x), npsi, psi, substitute(by),
                      list(...), smoothed, !missing(c0), !missing(n_boot),
                      label, linear = FALSE))
  }
  if (!is.null(smoothed) && !missing(c0)) {
    message(paste("'c0' is the scaling schedule of the working",
                  "parametrization, and a smoothed term has none:",
                  "it is ignored."))
  }
  .seg_spec("jump", substitute(x), npsi, psi, FALSE, c0, n_boot, label,
            list(...), .seg_by_formula(substitute(by), parent.frame()),
            smoothed)
}

#' Segmented-with-Jump Term: Slope and Level Both Changing
#'
#' @description
#' A covariate whose slope *and* level both change at \eqn{K}
#' break-points estimated with everything else, the union of
#' [seg()] and [jump()] at the same points. Written in
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
#' schedule provides (see [jump()]). What is particular to the
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
#' so the pure step is a case this contains, never an exception to
#' it.
#'
#' Everything [jump()] documents about the scaling schedule, the
#' confinement of a break-point and the choice of starting positions
#' applies here unchanged.
#'
#' @section The coefficients:
#' `beta` (present when `linear = TRUE`), `gamma1` \ldots
#' `gammaK`, the changes of slope, then `delta1` \ldots
#' `deltaK`, the changes of level, then `g1` \ldots `gK`,
#' the auxiliary coefficients from which the break-points are read. Ask
#' [seg_psi()] for the break-points.
#'
#' @section Developing a coefficient with covariates:
#' `beta` and `gammaK` take a development on any design, both
#' entering the block linearly. `deltaK` and `psiK` enter the
#' quadratic above, which is a scalar identity: it splits observation by
#' observation only where the design has one column per group and each
#' observation belongs to one, so `jseg(x, by = ~0 + g)` is an
#' independent set of everything per level and a development of
#' `delta` or `psi` on any other design is rejected. The
#' componentwise reading that would remain diverges whenever a change of
#' level passes near zero mid-iteration, which on a joint model it
#' routinely does, so it is rejected instead of shipped diverging.
#'
#' @section The marginal construction:
#' `jseg(x, psi ~ random(~1 | g), marginal = TRUE)` integrates a
#' latent break-point per group out of the likelihood exactly, on the
#' Gauss-Kronrod panels of [seg()]'s marginal with the change of
#' level entering each interval's conditional as a constant. On gaussian
#' data the mode-based routes are measured equivalent; on non-gaussian
#' families the marginal is the robust one, the measured comparison on a
#' Poisson panel recovering the per-group positions where the smoothed
#' modes lose them.
#'
#' @inheritParams seg
#' @inheritParams jump
#'
#' @return An object of class [SegTerm()] (a specification; see
#'   [term_build()]), or of class [MarginalBreakTerm()]
#'   with `marginal = TRUE`.
#'
#' @references
#' Muggeo, V. M. R. (2003). Estimating regression models with unknown
#' break-points. *Statistics in Medicine*, 22(19), 3055--3071.
#'
#' Fasola, S., Muggeo, V. M. R. and Kuchenhoff, H. (2018). A heuristic,
#' iterative algorithm for change-point detection in abrupt change
#' models. *Computational Statistics*, 33, 997--1015.
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
#' @seealso [seg()], [jump()], [seg_psi()],
#'   [seg_start()], [seg_step()]
#' @export
jseg <- function(x, ..., npsi = 1, psi = NULL, by = NULL, linear = TRUE,
                 c0 = 0.05, smoothed = NULL, marginal = FALSE, n_boot = 10,
                 label = "jseg") {
  .seg_check_marginal(marginal)
  if (marginal) {
    return(.marg_spec("jseg", substitute(x), npsi, psi, substitute(by),
                      list(...), smoothed, !missing(c0), !missing(n_boot),
                      label, linear = linear))
  }
  if (!is.null(smoothed) && !missing(c0)) {
    message(paste("'c0' is the scaling schedule of the working",
                  "parametrization, and a smoothed term has none:",
                  "it is ignored."))
  }
  .seg_spec("jseg", substitute(x), npsi, psi, linear, c0, n_boot, label,
            list(...), .seg_by_formula(substitute(by), parent.frame()),
            smoothed)
}

.seg_check_marginal <- function(marginal) {
  if (!is.logical(marginal) || length(marginal) != 1L || is.na(marginal)) {
    stop("'marginal' must be TRUE or FALSE.", call. = FALSE)
  }
  invisible(NULL)
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

.seg_spec <- function(kind, var, npsi, psi, linear, c0, n_boot, label, dots,
                      by, smoothed = NULL) {
  if (!is.null(smoothed) &&
      !S7::S7_inherits(smoothed, penalties7::abs_smoother)) {
    stop(paste("'smoothed' must be NULL or a penalties7 abs_smoother, e.g.",
               "penalties7::smooth_probit()."), call. = FALSE)
  }
  if (!is.numeric(c0) || length(c0) != 1L || is.na(c0) ||
      c0 <= 0 || c0 >= 1) {
    stop("'c0' must be a single number strictly between 0 and 1.",
         call. = FALSE)
  }
  if (!is.numeric(n_boot) || length(n_boot) != 1L || is.na(n_boot) ||
      n_boot < 0 || n_boot != round(n_boot)) {
    stop("'n_boot' must be a single non-negative integer.", call. = FALSE)
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
          spec = list(psi = psi, c0 = c0, n_boot = as.integer(n_boot),
                      smoothed = smoothed),
          X = NULL, coef_names = character(0),
          blueprint = list(), penalty = NULL)
}

# The stem of each own coefficient's name. It is the parameter's own name
# except for a break-point of a discontinuous construction, whose
# coefficients are the auxiliary g of the identity of Fasola et al. and not
# the position: the position is read off, and seg_psi() is what reports it.
# A SMOOTHED term holds every break-point directly, whatever its kind, so
# the stem stays psi there.
.seg_coef_stem <- function(kind, p, smoothed = FALSE) {
  if (!smoothed && kind != "seg" && grepl("^psi[0-9]+$", p)) {
    sub("^psi", "g", p)
  } else {
    p
  }
}

# A design with one column per group and every observation in exactly one
# of them. The product of two developments over such a design collapses
# group by group, which is what lets a discontinuous construction carry
# both a change of level and a break-point that vary.
.seg_is_partition <- function(Z) {
  Z <- as.matrix(Z)
  all(Z == 0 | Z == 1) && all(rowSums(Z) == 1)
}

# The columns of a development's design that partition the observations,
# after dropping constant columns: `psi ~ random(~1 | id)` carries the
# subformula's unpenalized intercept beside the indicators, and the
# grouping a width is resolved within is the indicators'. NULL where no
# such set exists (a numeric development, say).
.seg_partition_cols <- function(Z) {
  Z <- as.matrix(Z)
  keep <- which(apply(Z, 2L, function(z) {
    all(z %in% c(0, 1)) && any(z != z[1L])
  }))
  if (!length(keep)) return(NULL)
  if (all(rowSums(Z[, keep, drop = FALSE]) == 1)) keep else NULL
}

# What each construction can carry, and why. The continuous one accepts a
# development of anything: every coefficient enters its block linearly. The
# discontinuous ones read the break-point off g = -delta * psi, a product
# of two of the unknowns, so it stays linear only where one factor is a
# single number or where the design collapses the product; and the joint
# construction reads it off a quadratic that also carries the change of
# slope, so that one is constrained as well.
.seg_check_devs <- function(kind, npsi, dev, smoothed = FALSE) {
  # a smoothed term holds every break-point directly, so there is no
  # read-off to constrain and a development of anything is legal, as it is
  # for the continuous construction
  if (kind == "seg" || smoothed) return(invisible(NULL))
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
  smoothed <- !is.null(term@spec$smoothed)
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
  .seg_check_devs(term@kind, term@npsi, dev, smoothed)
  if (term@kind != "seg" && !smoothed) {
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
  # a smoothed term holds the break-point directly whatever its kind
  if (bp$kind == "seg" || !is.null(bp$smooth)) return(cf)
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

# The smoothed derivatives at one set of points: s^(order) with the width
# the build resolved, per observation where it was resolved per group. On
# other rows a per-group width is recomputed from the break-point
# development's design, which is what carries the grouping.
.seg_sw <- function(sm, u, w, order) {
  sm@s[[order + 1L]](u, w)
}

.seg_smooth_w <- function(sm, n, Z) {
  if (is.null(sm$w_group)) return(rep(sm$width, length.out = n))
  z <- as.matrix(Z[[sm$group_param]])[, sm$group_cols, drop = FALSE]
  as.numeric(z %*% sm$w_group)
}

# The block and the contribution of a SMOOTHED term: the break-points are
# ordinary parameters and the block is the true Jacobian of
#
#   f = beta x + sum_k [ gamma_k H(u_k) + delta_k S(u_k) ],  u_k = x - psi_k,
#
# with H(u) = (u + s(u))/2 the smooth hinge and S(u) = (1 + s'(u))/2 the
# smooth step. The columns are
#
#   beta     x
#   gamma_k  H(u_k)
#   delta_k  S(u_k)
#   psi_k    -( gamma_k S(u_k) + delta_k s''(u_k)/2 )
#
# since H' = S and S' = s''/2. Everything is the closed form of the
# smoother's own derivatives, so the same machinery serves the three kinds
# with a coefficient absent simply contributing nothing.
.seg_block_smooth <- function(bp, xv, coef, Z = bp$Z, prev = bp$pk) {
  n <- length(xv)
  K <- bp$npsi
  sm <- bp$smooth
  w <- .seg_smooth_w(sm, n, Z)
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
  has_g <- bp$kind %in% c("seg", "jseg")
  has_d <- bp$kind %in% c("jump", "jseg")
  for (k in seq_len(K)) {
    u <- xv - psi[, k]
    s0 <- .seg_sw(sm$sm, u, w, 0L)
    s1 <- .seg_sw(sm$sm, u, w, 1L)
    s2 <- .seg_sw(sm$sm, u, w, 2L)
    H <- (u + s0) / 2
    S <- (1 + s1) / 2
    dpsi <- numeric(n)
    if (has_g) {
      gam <- val(paste0("gamma", k))
      put(paste0("gamma", k), H)
      value <- value + gam * H
      dpsi <- dpsi + gam * S
    }
    if (has_d) {
      del <- val(paste0("delta", k))
      put(paste0("delta", k), S)
      value <- value + del * S
      dpsi <- dpsi + del * s2 / 2
    }
    put(paste0("psi", k), -dpsi)
  }
  list(X = .nl_bind(cols[bp$params]), value = value, psi = psi, pk = pos$pk)
}

.seg_assemble <- function(bp, xv, coef, cscale = bp$cscale, Z = bp$Z,
                          prev = bp$pk) {
  if (!is.null(bp$smooth)) .seg_block_smooth(bp, xv, coef, Z, prev)
  else if (bp$developed) .seg_block(bp, xv, coef, cscale, Z, prev)
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

#' @title Build a Break-Point Term
#' @name term_build.SegTerm
#'
#' @description
#' Evaluates the covariate, resolves any development of the term's own
#' coefficients, places the break-points at their starting positions and builds
#' the working block there. It also settles the confinement limits, the
#' scaling schedule and, where the term is smoothed, the transition width.
#'
#' @details
#' # What is settled here
#'
#' The covariate must vary. The break-points are confined to the interval
#' between the 5th and the 95th percentile of it: outside that the indicator is
#' constant, the truncated line and that constant are linearly dependent, and
#' the block is exactly singular. A run ending against the limit has found no
#' break-point, and [seg_psi()] then reports the limit itself.
#'
#' The starting positions are `psi` where it was given and evenly spaced
#' quantiles of the covariate otherwise. Zero is degenerate here, so
#' [term_coef_start()] returns these rather than a vector of zeros.
#'
#' Each subformula's right-hand side goes through [interpret_formula()] and its
#' terms are built, so their blueprints are recorded and reapplied at
#' prediction. What a development may carry depends on `kind`: the continuous
#' construction takes any of them, and the discontinuous ones only a scalar
#' change of level or a design that partitions the rows, their read-off being a
#' product of two unknowns.
#'
#' # The smoothed form
#'
#' With `smoothed` the step and the hinge are replaced by a
#' [penalties7::abs_smoother()]'s versions, and the transition width is
#' resolved here from the covariate's spacing, within groups where a
#' development of the break-point supplies a partition. The width is checked
#' against the derived floor \eqn{\sqrt{\epsilon}D} and reported by [print()],
#' a number that changes what the term means having to be legible.
#'
#' @param term A [SegTerm()], built or not.
#' @param data A data frame carrying the covariate and whatever the
#'   subformulas name.
#' @param ... Unused.
#'
#' @return The term with `X`, `coef_names`, `blueprint` and, where a
#'   subformula brought one, `penalty` filled.
#'
#' @seealso [seg()], [jump()], [jseg()]; [term_refresh.SegTerm()] for the block
#'   as the break-points move; [seg_start()] for a better starting position
#'   than the quantiles.
#'
#' @examples
#' set.seed(1)
#' d <- data.frame(x = sort(runif(120, 0, 10)), id = factor(rep(1:4, each = 30)))
#' d$y <- 1 + 0.5 * d$x + 2 * pmax(d$x - 6, 0) + rnorm(120, sd = 0.4)
#'
#' b <- term_build(seg(x, npsi = 1), d)
#' term_coef_names(b)
#' seg_psi(b)                     # the starting quantile
#' setNames(term_coef_start(b), term_coef_names(b))
#'
#' # A developed break-point expands into its own design's coefficients.
#' bd <- term_build(seg(x, psi ~ id), d)
#' term_coef_names(bd)
#' lapply(term_components(bd), function(z) z$index)
#'
#' @keywords internal
S7::method(term_build, SegTerm) <- function(term, data, ...) {
  xv <- .seg_x(term@var, data)
  n <- length(xv)
  params <- .seg_params(term@kind, term@npsi, term@linear)
  smoothed <- term@spec$smoothed
  dev <- .seg_designs(term, data)

  rr <- range(xv)
  if (diff(rr) <= 0) {
    stop("the covariate of a break-point term must vary.", call. = FALSE)
  }

  # The smoother's width, resolved from the covariate's spacing where the
  # object carries none: the median gap between distinct values, within
  # groups where a break-point development supplies a partition -- the
  # validity window of a Laplace approximation being per-subject -- and
  # globally otherwise. The floor is derived, not chosen: see
  # penalties7::smoother_width_floor().
  smooth <- NULL
  if (!is.null(smoothed)) {
    gaps_of <- function(v) {
      u <- sort(unique(v))
      if (length(u) > 1L) diff(u) else numeric(0)
    }
    # the covariate varies (checked above), so the global gap set is
    # non-empty
    global_sp <- stats::median(gaps_of(xv))
    gp <- NULL
    gcols <- NULL
    for (k in seq_len(term@npsi)) {
      d <- dev[[paste0("psi", k)]]
      if (is.null(d)) next
      pc <- .seg_partition_cols(d$Z)
      if (!is.null(pc)) {
        gp <- paste0("psi", k)
        gcols <- pc
        break
      }
    }
    per <- NULL
    if (!is.null(gp)) {
      Zg <- as.matrix(dev[[gp]]$Z)[, gcols, drop = FALSE]
      per <- vapply(seq_len(ncol(Zg)), function(j) {
        g <- gaps_of(xv[Zg[, j] == 1])
        if (length(g)) stats::median(g) else global_sp
      }, numeric(1))
    }
    floorw <- penalties7::smoother_width_floor(smoothed, diff(rr))
    if (isTRUE(smoothed@per_group)) {
      if (is.null(per)) {
        stop(paste("a per-group width needs a break-point development whose",
                   "design has one column per group with each observation in",
                   "one of them, as psi ~ 0 + g or psi ~ random(~1 | g)",
                   "gives."), call. = FALSE)
      }
      wg <- if (is.null(smoothed@width)) {
        penalties7::smoother_width(smoothed, per)
      } else {
        rep(smoothed@width, length(per))
      }
      if (any(wg < floorw)) {
        stop(sprintf(paste("the smoother's width (%.3g) is below the floor",
                           "%.3g = the width at which the Jacobian column,",
                           "of order 1/h against the covariate's range,",
                           "takes the design's condition past eps^-1/2."),
                     min(wg), floorw), call. = FALSE)
      }
      smooth <- list(sm = smoothed, w_group = wg, group_param = gp,
                     group_cols = gcols, width = stats::median(wg))
    } else {
      sp <- if (!is.null(per)) stats::median(per) else global_sp
      wd <- penalties7::smoother_width(smoothed, sp)
      if (wd < floorw) {
        stop(sprintf(paste("the smoother's width (%.3g) is below the floor",
                           "%.3g = the width at which the Jacobian column,",
                           "of order 1/h against the covariate's range,",
                           "takes the design's condition past eps^-1/2."),
                     wd, floorw), call. = FALSE)
      }
      smooth <- list(sm = smoothed, w_group = NULL, group_param = NULL,
                     group_cols = NULL, width = wd)
    }
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
    stem <- .seg_coef_stem(term@kind, p, !is.null(smoothed))
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
             if (term@kind == "seg" || !is.null(smoothed)) start_psi[k]
             else -start_psi[k])
  }

  bp <- list(kind = term@kind, npsi = term@npsi, linear = term@linear,
             params = params, Z = Z, index = index,
             subs = lapply(dev, function(d) if (is.null(d)) NULL else d$subs),
             developed = any(!vapply(dev, is.null, logical(1))),
             lo = rr[1L], hi = rr[2L], lim = lim, delta = delta,
             cscale = rep(term@spec$c0, term@npsi),
             sgn = rep(0, term@npsi), step = rep(NA_real_, term@npsi),
             nref = 0L, var = term@var, coef = coef0, xv = xv,
             pk = vector("list", term@npsi), smooth = smooth)

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
#' `gamma ~ 0 + lasso(~1)` is a single penalized block over the changes
#' under one hyperparameter. The list is empty for a term whose
#' coefficients carry no development, and for a specification, whose
#' coefficients there is nothing yet to index: a penalty is attached at
#' build, as it is for every penalized term here.
#' @param term A built [SegTerm()].
#' @param ... Unused.
#' @return A list of entries, as [term_penalties()] documents.
#' @keywords internal
S7::method(term_penalties, SegTerm) <- function(term, ...) {
  pens <- term@blueprint$penalties
  if (is.null(pens)) list() else pens
}

# The continuous construction's block is the Jacobian of its contribution;
# the discontinuous ones carry the frozen weight of the identity of Fasola
# et al., so their block is a working linearization and the fixed-point
# iteration it belongs to is not a descent method on the model's objective.
# A SMOOTHED term is a Jacobian whatever its kind: the break-points are
# ordinary parameters and there is no frozen weight. The question is asked
# of the specification, which a fitting layer routes on before the term is
# built.
#' @title Whether a Break-Point Term's Block Is a Jacobian
#' @name term_jacobian_block.SegTerm
#'
#' @description
#' `TRUE` for [seg()] and for any smoothed construction, `FALSE` for the sharp
#' [jump()] and [jseg()]. It is the predicate a fitting layer routes on: a
#' Jacobian licenses a Gauss-Newton step and a line search on the model's own
#' objective, and a working linearization does not.
#'
#' @details
#' [seg()] is differentiable in its break-points away from them, so its block
#' is the exact derivative. [jump()] and [jseg()] freeze the weight
#' \eqn{W = 1/(2\lvert\tilde{x}-\psi\rvert)} at the previous position and read
#' the break-point off two coefficients, so their block is a working
#' linearization: the fixed-point iteration on it is not a descent method on
#' the model's objective, its early steps going uphill on purpose under a large
#' scaling factor, and forcing a sufficient decrease stalls it.
#'
#' `smoothed` changes the answer for all three. Replacing the step and the
#' hinge by an [penalties7::abs_smoother()]'s smooth versions makes every
#' break-point an ordinary parameter with a true derivative, so a smoothed
#' `jump` answers `TRUE` and is routed like an [nl()] term.
#'
#' @param term A [SegTerm()], built or not: the answer is a property of the
#'   construction.
#' @param ... Unused.
#'
#' @return A single logical.
#'
#' @seealso [term_jacobian_block()] for the generic and what the answer
#'   decides, [term_converged()] for the verdict a working linearization gets
#'   instead of a score.
#'
#' @examples
#' # The continuous construction differentiates; the two discontinuous
#' # ones report a working linearization.
#' vapply(list(seg = seg(x), jump = jump(x), jseg = jseg(x)),
#'        term_jacobian_block, logical(1))
#'
#' # Smoothing the step makes every break-point an ordinary parameter.
#' vapply(list(seg = seg(x, smoothed = penalties7::smooth_probit()),
#'             jump = jump(x, smoothed = penalties7::smooth_probit()),
#'             jseg = jseg(x, smoothed = penalties7::smooth_probit())),
#'        term_jacobian_block, logical(1))
#'
#' @keywords internal
S7::method(term_jacobian_block, SegTerm) <- function(term, ...) {
  identical(term@kind, "seg") || !is.null(term@spec$smoothed)
}

# How a smoothed block moves with the coefficients, in the two shapes the
# generics ask: contracted against a caller's weight (per coefficient) and
# taken along one direction (per entry of the block). Everything is the
# closed form of the smoother's derivatives one order up from the block --
# with u = x - psi, S = (1 + s')/2, P = s''/2 and T = s'''/2,
#
#   d(gamma col)/d psi = -S      d(psi col)/d gamma = -S
#   d(delta col)/d psi = -P      d(psi col)/d delta = -P
#   d(psi col)/d psi   = gamma P + delta T
#
# with the derivative in a break-point's own coefficients gated by `free`
# where the position sits against a confinement limit, exactly as the
# continuous construction gates its indicator.
# `order` is the highest derivative of the smoother wanted: three for the
# first-order generics, four for term_block_deriv2(), which needs one more.
# The fourth is not computed by default because these run at every step the
# fit takes and nothing on that path reads it.
.seg_smooth_parts <- function(bp, cf, k, order = 3L) {
  n <- length(bp$xv)
  sm <- bp$smooth
  wv <- .seg_smooth_w(sm, n, bp$Z)
  pos <- .seg_positions(bp, cf, n)
  u <- bp$xv - pos$psi[, k]
  pk <- paste0("psi", k)
  zp <- bp$Z[[pk]]
  v <- if (is.null(zp)) rep(pos$pk[[k]], n) else
    as.numeric(as.matrix(zp) %*% pos$pk[[k]])
  out <- list(
    S = (1 + .seg_sw(sm$sm, u, wv, 1L)) / 2,
    P = .seg_sw(sm$sm, u, wv, 2L) / 2,
    T = .seg_sw(sm$sm, u, wv, 3L) / 2,
    free = as.numeric(v > bp$lim[1L] & v < bp$lim[2L])
  )
  if (order >= 4L) out$Q <- .seg_sw(sm$sm, u, wv, 4L) / 2
  out
}

.seg_smooth_contract <- function(bp, cf, A, out) {
  n <- length(bp$xv)
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
  val <- function(p) .seg_pval(list(Z = bp$Z, index = bp$index), cf, p, n)
  has_g <- bp$kind %in% c("seg", "jseg")
  has_d <- bp$kind %in% c("jump", "jseg")
  for (k in seq_len(bp$npsi)) {
    pk <- paste0("psi", k)
    pt <- .seg_smooth_parts(bp, cf, k)
    dself <- numeric(n)
    if (has_g) {
      gk <- paste0("gamma", k)
      place(pk, -pt$S * sfun(gk) * pt$free)
      place(gk, -pt$S * sfun(pk))
      dself <- dself + val(gk) * pt$P
    }
    if (has_d) {
      dk <- paste0("delta", k)
      place(pk, -pt$P * sfun(dk) * pt$free)
      place(dk, -pt$P * sfun(pk))
      dself <- dself + val(dk) * pt$T
    }
    place(pk, dself * sfun(pk) * pt$free)
  }
  out
}

.seg_smooth_deriv <- function(bp, cf, v, out) {
  n <- length(bp$xv)
  along <- function(p) {
    idx <- bp$index[[p]]
    z <- bp$Z[[p]]
    if (is.null(z)) rep(v[idx], n) else as.numeric(as.matrix(z) %*% v[idx])
  }
  put <- function(p, w) {
    idx <- bp$index[[p]]
    z <- bp$Z[[p]]
    out[, idx] <<- out[, idx] + if (is.null(z)) w else as.matrix(z) * w
  }
  val <- function(p) .seg_pval(list(Z = bp$Z, index = bp$index), cf, p, n)
  has_g <- bp$kind %in% c("seg", "jseg")
  has_d <- bp$kind %in% c("jump", "jseg")
  for (k in seq_len(bp$npsi)) {
    pk <- paste0("psi", k)
    pt <- .seg_smooth_parts(bp, cf, k)
    dself <- numeric(n)
    if (has_g) {
      gk <- paste0("gamma", k)
      put(gk, -pt$S * pt$free * along(pk))
      put(pk, -pt$S * along(gk))
      dself <- dself + val(gk) * pt$P
    }
    if (has_d) {
      dk <- paste0("delta", k)
      put(dk, -pt$P * pt$free * along(pk))
      put(pk, -pt$P * along(dk))
      dself <- dself + val(dk) * pt$T
    }
    put(pk, dself * pt$free * along(pk))
  }
  out
}

# The block's SECOND derivative in the coefficients, contracted in two
# directions. It is the same closed forms one order further up, and the only
# new quantity is Q = s''''/2: the block reads the smoother to order two, the
# first derivative to order three, this to order four.
#
# With u = x - psi and the columns H(u), S(u) and -(gamma S + delta P), and
# writing pv and pu for the two directions carried onto the break-point --
# each with the confinement gate, one factor per differentiation in a
# break-point coefficient --
#
#   d2(gamma col) = P pv pu                 (H'' = S' = P)
#   d2(delta col) = T pv pu                 (S'' = P' = T)
#   d2(psi col)   = P (gamma_v pu + gamma_u pv) + T (delta_v pu + delta_u pv)
#                   - (gamma T + delta Q) pv pu
#
# The psi column's three pieces come from differentiating its own first
# derivative -gamma_v S - delta_v P + (gamma P + delta T) pv once more: the
# first two move only through u, and the third moves through gamma, delta and
# u at once. Exchanging the two directions maps each piece to itself, which
# is the symmetry the generic promises.
.seg_smooth_deriv2 <- function(bp, cf, v, u, out) {
  n <- length(bp$xv)
  along <- function(vec, p) {
    idx <- bp$index[[p]]
    z <- bp$Z[[p]]
    if (is.null(z)) rep(vec[idx], n) else as.numeric(z %*% vec[idx])
  }
  put <- function(p, w) {
    idx <- bp$index[[p]]
    z <- bp$Z[[p]]
    out[, idx] <<- out[, idx] + if (is.null(z)) w else as.matrix(z * w)
  }
  val <- function(p) .seg_pval(list(Z = bp$Z, index = bp$index), cf, p, n)
  has_g <- bp$kind %in% c("seg", "jseg")
  has_d <- bp$kind %in% c("jump", "jseg")
  for (k in seq_len(bp$npsi)) {
    pk <- paste0("psi", k)
    pt <- .seg_smooth_parts(bp, cf, k, order = 4L)
    pv <- pt$free * along(v, pk)
    pu <- pt$free * along(u, pk)
    dself <- numeric(n)
    if (has_g) {
      gk <- paste0("gamma", k)
      put(gk, pt$P * pv * pu)
      put(pk, pt$P * (along(v, gk) * pu + along(u, gk) * pv))
      dself <- dself + val(gk) * pt$T
    }
    if (has_d) {
      dk <- paste0("delta", k)
      put(dk, pt$T * pv * pu)
      put(pk, pt$T * (along(v, dk) * pu + along(u, dk) * pv))
      dself <- dself + val(dk) * pt$Q
    }
    put(pk, -dself * pv * pu)
  }
  out
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
  # brute-force 1.30e+02. A SMOOTHED term of any kind is the exception: its
  # block is a true Jacobian and its second derivatives are the smoother's
  # own, one order up.
  if (!is.null(bp$smooth)) return(.seg_smooth_contract(bp, cf, A, out))
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


S7::method(term_block_deriv, SegTerm) <- function(term, coef = NULL, v, ...) {
  bp <- term@blueprint
  if (!length(bp)) stop("the term is not built.", call. = FALSE)
  cf <- if (is.null(coef)) bp$coef else as.numeric(coef)
  ncoef <- length(term_coef_names(term))
  n <- length(bp$xv)
  v <- as.numeric(v)
  if (length(v) != ncoef) {
    stop(sprintf("'v' must have length %d, the term's coefficients.", ncoef),
         call. = FALSE)
  }
  out <- matrix(0, n, ncoef)
  # the continuous construction only, for the reason term_block_contract()
  # records: a discontinuous one reads its position off a product of the
  # unknowns and carries a weight whose derivative in the break-point is
  # unbounded. A smoothed term of any kind has the closed forms instead.
  if (!is.null(bp$smooth)) return(.seg_smooth_deriv(bp, cf, v, out))
  if (!identical(bp$kind, "seg")) return(out)

  along <- function(p) {
    idx <- bp$index[[p]]
    z <- bp$Z[[p]]
    if (is.null(z)) rep(v[idx], n) else as.numeric(as.matrix(z) %*% v[idx])
  }
  put <- function(p, w) {
    idx <- bp$index[[p]]
    z <- bp$Z[[p]]
    out[, idx] <<- out[, idx] + if (is.null(z)) w else as.matrix(z) * w
  }

  xv <- bp$xv
  pos <- .seg_positions(bp, cf, n)
  for (k in seq_len(bp$npsi)) {
    pk <- paste0("psi", k)
    gk <- paste0("gamma", k)
    psi <- pos$psi[, k]
    ind <- as.numeric(xv > psi)
    zp <- bp$Z[[pk]]
    vv <- if (is.null(zp)) rep(pos$pk[[k]], n) else
      as.numeric(as.matrix(zp) %*% pos$pk[[k]])
    free <- as.numeric(vv > bp$lim[1L] & vv < bp$lim[2L])
    # the change column moves with the position, and the position column with
    # the change: the same two pieces term_block_contract() carries, read the
    # other way round
    put(gk, -ind * free * along(pk))
    put(pk, -ind * along(gk))
  }
  out
}

#' @title How a Break-Point Term's Columns Divide
#' @name term_components.SegTerm
#'
#' @description
#' One entry per coefficient of the term, `beta`, `gamma1` ... , `delta1` ... ,
#' `psi1` ... , giving the columns that coefficient owns and the sub-terms
#' developing it. A coefficient with no development owns one column; a
#' developed one owns as many as its own design has.
#'
#' @details
#' It is what lets a consumer report a fitted break-point term coefficient by
#' coefficient, and what [term_penalties()] uses to place a sub-term's penalty
#' on exactly the coordinates it covers. Those coordinates are **named** as a
#' subset of the term's own parameters and never selected with a map: a
#' separable penalty under a selection map is the generalized-lasso problem,
#' which has no proximal operator.
#'
#' @param term A built [SegTerm()]. An unbuilt one gives an empty list.
#' @param ... Unused.
#'
#' @return A named list, one entry per coefficient and named by it, each with
#'   `name`, `index`, `subs` and `sub_index` as [term_components()] describes.
#'
#' @seealso [term_components()] for the contract, [term_penalties()] for the
#'   penalties the sub-terms bring.
#'
#' @examples
#' set.seed(1)
#' d <- data.frame(x = sort(runif(120, 0, 10)), id = factor(rep(1:4, each = 30)))
#' d$y <- 1 + 0.5 * d$x + 2 * pmax(d$x - 6, 0) + rnorm(120, sd = 0.4)
#'
#' # A break-point per subject: psi1 owns four columns, the others one each.
#' b <- term_build(seg(x, psi ~ id), d)
#' term_coef_names(b)
#' lapply(term_components(b), function(z) z$index)
#'
#' @keywords internal
S7::method(term_components, SegTerm) <- function(term, ...) {
  bp <- term@blueprint
  if (!length(bp)) return(list())
  # the same division the block is assembled from: `index` per parameter --
  # beta, the changes, the break-points -- and `subs` where one of them is
  # developed. A break-point developed by `psi ~ random(~1 | g)` carries TWO
  # sub-terms, an unpenalized intercept and the random block, and the entry
  # reports both.
  out <- lapply(bp$params, function(p) {
    s <- bp$subs[[p]]
    ss <- if (is.null(s)) list() else s
    list(name = p, index = bp$index[[p]], subs = ss,
         sub_index = component_sub_index(bp$index[[p]], ss))
  })
  stats::setNames(out, bp$params)
}

S7::method(term_block_deriv2, SegTerm) <- function(term, coef = NULL, v, u,
                                                   ...) {
  bp <- term@blueprint
  if (!length(bp)) stop("the term is not built.", call. = FALSE)
  cf <- if (is.null(coef)) bp$coef else as.numeric(coef)
  ncoef <- length(term_coef_names(term))
  n <- length(bp$xv)
  v <- as.numeric(v)
  u <- as.numeric(u)
  for (arg in list(list(v, "v"), list(u, "u"))) {
    if (length(arg[[1L]]) != ncoef) {
      stop(sprintf("'%s' must have length %d, the term's coefficients.",
                   arg[[2L]], ncoef), call. = FALSE)
    }
  }
  out <- matrix(0, n, ncoef)
  # A SMOOTHED term of any kind has the closed forms, its block being a true
  # Jacobian; the sharp constructions answer zeros, and the two reasons are
  # different. For the continuous one the second derivative really is zero
  # away from the break-points: the truncated line's derivative in the
  # position is the indicator, whose own derivative is a point mass, and the
  # position column is linear in the change. For a jump or a jseg the block is
  # a working LINEARIZATION with a frozen weight rather than a Jacobian, so
  # the question the generic asks is not the one their columns answer -- which
  # is why term_block_deriv() and term_block_contract() already answer zeros
  # there.
  if (!is.null(bp$smooth)) return(.seg_smooth_deriv2(bp, cf, v, u, out))
  out
}


# The relabeling itself: NULL where there is nothing to do. The per-k
# coefficient stretches are permuted along with every per-break-point
# field of the blueprint, so a lineage keeps its own previous position,
# its scaling factor and its direction under the new label.
.seg_reorder <- function(bp, coef, psi_new) {
  K <- bp$npsi
  if (K < 2L) return(NULL)
  perk <- c(if (bp$kind %in% c("seg", "jseg")) paste0("gamma", seq_len(K)),
            if (bp$kind %in% c("jump", "jseg")) paste0("delta", seq_len(K)),
            paste0("psi", seq_len(K)))
  if (any(vapply(perk, function(p) !is.null(bp$Z[[p]]), logical(1)))) {
    return(NULL)
  }
  pos <- psi_new[1L, ]
  o <- order(pos)
  if (!is.unsorted(pos, strictly = FALSE)) return(NULL)
  relabel <- function(pref) {
    for (j in seq_len(K)) {
      coef[bp$index[[paste0(pref, j)]]] <<-
        old[bp$index[[paste0(pref, o[j])]]]
    }
  }
  old <- coef
  if (bp$kind %in% c("seg", "jseg")) relabel("gamma")
  if (bp$kind %in% c("jump", "jseg")) relabel("delta")
  relabel("psi")
  bp$psi <- bp$psi[, o, drop = FALSE]
  bp$pk <- bp$pk[o]
  bp$cscale <- bp$cscale[o]
  bp$sgn <- bp$sgn[o]
  list(coef = coef, bp = bp, psi_new = psi_new[, o, drop = FALSE])
}

#' @title Recompute a Break-Point Term at New Coefficients
#' @name term_refresh.SegTerm
#'
#' @description
#' Moves the break-points to where the coefficients put them, rebuilds the
#' working block there, and advances the scaling schedule. It is one step of
#' the iteration that estimates the break-points, and the block it returns is
#' what the next linear fit is taken on.
#'
#' @details
#' # How the new positions are found
#'
#' For `seg` the break-point is a coefficient, so the new position is read
#' straight out of `coef`. For `jump` and `jseg` it is **read off** two of
#' them: the identity \eqn{1(x > \psi) = 1/2 + (x-\psi)/(2|x-\psi|)} makes the
#' step linear in \eqn{\psi} once the weight is frozen, and
#' \eqn{\psi = -g_k/\delta_k} follows. `jseg` reads a quadratic instead, its
#' truncated line depending on the position as well.
#'
#' Every new position is clamped into the confinement interval the build
#' settled, and crossed break-points are relabeled so that each keeps
#' travelling with its own scaling factor and direction.
#'
#' # The scaling schedule
#'
#' The rescaling of Fasola, Muggeo and Kuchenhoff (2018) opens a gap of
#' relative width \eqn{c} around each break-point, and the factor is **halved
#' whenever that break-point reverses direction**, which is the signal that the
#' iteration has begun to circle an optimum instead of traveling toward one.
#' The factor is both the conditioning device and the step control, so it must
#' advance once per committed step: refreshing inside a line search would
#' anneal at every trial point, and refreshing from the specification would
#' freeze it at `c0`. [seg_reheat()] resets it.
#'
#' @param term A built [SegTerm()].
#' @param coef The coefficients to move to, of length `ncol(term@X)`. Any
#'   other length throws with the required length named.
#' @param ... Unused.
#'
#' @return The term with its block, its break-points, its contribution and its
#'   scaling schedule all recomputed at `coef`.
#'
#' @seealso [term_refresh()] for the generic, [seg_psi()] for the positions
#'   after a refresh, [seg_step()] and [term_converged()] for the verdict,
#'   [seg_reheat()] to start the schedule again.
#'
#' @examples
#' set.seed(1)
#' d <- data.frame(x = sort(runif(120, 0, 10)))
#' d$y <- 1 + 0.5 * d$x + 2 * pmax(d$x - 6, 0) + rnorm(120, sd = 0.4)
#' b <- term_build(seg(x, npsi = 1), d)
#'
#' # Moving the break-point coefficient moves the block.
#' cf <- b@blueprint$coef
#' r <- term_refresh(b, c(cf[1], cf[2], 6))
#' c(before = seg_psi(b), after = seg_psi(r))
#' max(abs(term_matrix(r) - term_matrix(b))) > 0
#'
#' # A wrong length is refused.
#' try(term_refresh(b, c(1, 2)))
#'
#' @keywords internal
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
  # Two break-points that cross have exchanged roles, and left to
  # themselves they go on chasing the same feature: the fit that follows
  # reads two nearly identical columns and the iteration settles on a
  # collision. The lineages are relabeled so the positions come out
  # ordered -- each (change, level, position) triple travels with its own
  # scaling factor and direction -- which is what `segmented` does at
  # every iteration. The contribution is a sum over k, so relabeling
  # moves no value. A term whose per-break-point coefficients carry a
  # development is left alone: its positions are one per observation and
  # a single order need not exist.
  psi_new <- .seg_positions(bp, coef, length(bp$xv))$psi
  ro <- .seg_reorder(bp, coef, psi_new)
  if (!is.null(ro)) {
    coef <- ro$coef
    bp <- ro$bp
    psi_new <- ro$psi_new
  }
  dm <- psi_new - bp$psi
  d <- colMeans(dm)
  stepv <- apply(abs(dm), 2L, max)
  dfar <- apply(psi_new, 2L, function(p) max(pmax(p - bp$lo, bp$hi - p)))
  dnear <- apply(psi_new, 2L, function(p) min(pmin(p - bp$lo, bp$hi - p)))
  amax <- apply(abs(psi_new), 2L, max)

  # a smoothed term has no working parametrization and therefore no
  # schedule: the block is the Jacobian and the step is a Gauss-Newton one,
  # so only the step record (for seg_step) and the relabeling above apply
  if (!is.null(bp$smooth)) {
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
    return(term)
  }

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

#' @title The Contribution of a Break-Point Term
#' @name term_value.SegTerm
#'
#' @description
#' The values the term contributes to the linear predictor at its current
#' break-points: the broken line for [seg()], the step for [jump()], both for
#' [jseg()]. It is what a Gauss-Newton step needs beside the block, and for
#' the continuous construction it is **not** the block times the coefficients.
#'
#' @details
#' For `seg` the block is the Jacobian, so \eqn{X\beta} is the linearization
#' and the two differ by whatever that drops: a step at the break-point, in a
#' construction that is continuous. For `jump` the columns satisfy
#' \eqn{X\beta = } the contribution exactly, the identity behind the read-off
#' making them so, and the two agree to the last bit.
#'
#' With `newdata` the contribution is returned on those rows, at the
#' break-points the term **carries**. They are not re-derived from the new
#' rows, exactly as [term_predict()] reapplies rather than rebuilds.
#'
#' @param term A built [SegTerm()].
#' @param coef Optional coefficients to evaluate at; `NULL`, the default,
#'   uses the ones [term_refresh()] last committed.
#' @param newdata An optional data frame; the contribution comes back on its
#'   rows instead of the fitting ones.
#' @param ... Unused.
#'
#' @return A numeric vector of one value per observation.
#'
#' @seealso [term_value()] for the generic, [term_refresh.SegTerm()] for the
#'   block at the same coefficients, [seg_psi()] for the positions.
#'
#' @examples
#' set.seed(1)
#' d <- data.frame(x = sort(runif(120, 0, 10)))
#' d$y <- 1 + 0.5 * d$x + 2 * pmax(d$x - 6, 0) + rnorm(120, sd = 0.4)
#'
#' # seg: the block is the Jacobian, so X beta is the linearization and
#' # differs from the contribution.
#' b <- term_refresh(term_build(seg(x, npsi = 1), d), c(0.5, 2, 6))
#' max(abs(as.numeric(term_matrix(b) %*% b@blueprint$coef) - term_value(b)))
#'
#' # jump: the columns reproduce the contribution exactly.
#' bj <- term_build(jump(x, npsi = 1), d)
#' max(abs(as.numeric(term_matrix(bj) %*% bj@blueprint$coef) - term_value(bj)))
#'
#' # And the seg contribution really is the broken line at psi.
#' cf <- b@blueprint$coef
#' max(abs(term_value(b) - (cf[1] * d$x + cf[2] * pmax(d$x - cf[3], 0))))
#'
#' @keywords internal
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

#' @title A Break-Point Term's Block at New Rows
#' @name term_predict.SegTerm
#'
#' @description
#' The working block evaluated at `newdata`, at the break-points and the
#' coefficients the term currently carries. The positions are not re-derived
#' from the new rows, and any development of a coefficient is reapplied through
#' its sub-terms' own blueprints.
#'
#' @details
#' Re-deriving the break-points would give a different model at every set of
#' rows, which is the general reason [term_predict()] reapplies instead of
#' rebuilding. For a break-point term it is sharper than usual: the positions
#' are what the fit estimated, and a subset of the data may not even contain
#' one.
#'
#' [term_value()] with `newdata` is what gives the contribution there. For
#' `seg` the two differ, the block being a Jacobian.
#'
#' @param term A built [SegTerm()]. An unbuilt one throws
#'   `"the term has not been built; call term_build(term, data) first."`.
#' @param newdata A data frame carrying the covariate and whatever the
#'   subformulas name.
#' @param ... Unused.
#'
#' @return A block of `nrow(newdata)` rows and [term_npar()] columns, with the
#'   term's coefficient names as column names.
#'
#' @seealso [term_predict()] for the generic, [term_value()] for the
#'   contribution, [seg_psi()] for the positions it is evaluated at.
#'
#' @examples
#' set.seed(1)
#' d <- data.frame(x = sort(runif(120, 0, 10)))
#' d$y <- 1 + 0.5 * d$x + 2 * pmax(d$x - 6, 0) + rnorm(120, sd = 0.4)
#' b <- term_refresh(term_build(seg(x, npsi = 1), d), c(0.5, 2, 6))
#'
#' # On the fitting data it is the block itself.
#' max(abs(term_predict(b, d) - term_matrix(b)))
#'
#' # At other rows the break-point is the fitted one, not one re-derived
#' # from those rows: here every row is below it.
#' nd <- data.frame(x = c(1, 2, 3))
#' term_predict(b, nd)
#' seg_psi(b)
#'
#' @keywords internal
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
#' @param term A built [SegTerm()].
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
#' @seealso [seg()], [jump()], [jseg()]
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
#' `seg_step` returns how far each break-point moved at the last
#' call to [term_refresh()], and `seg_converged` compares
#' the largest of those with the tolerance of \cite{fasola2018}, a
#' hundredth of the distance between consecutive distinct observations of
#' the covariate. A term that has been built but not yet refreshed has
#' taken no step, so `seg_step` returns `NA` and
#' `seg_converged` returns `FALSE`.
#'
#' @details
#' With \eqn{x_{(1)} < \cdots < x_{(m)}} the distinct covariate values, the
#' run stops at
#'
#' \deqn{\max_k \lvert \psi_k^{(t)} - \psi_k^{(t-1)} \rvert < \Delta,
#'   \qquad \Delta = 0.01 \cdot
#'     \operatorname{median}_{i}\, (x_{(i+1)} - x_{(i)}).}
#'
#' \cite{fasola2018} take the smallest of those gaps instead of their
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
#' change instead, as `segmented` does, at a cost of a
#' continuous term nothing: its iteration can settle into a cycle of
#' period two in the break-point, in which case this rule is never met
#' while the objective has long since stopped moving.
#'
#' @param term A built [SegTerm()].
#'
#' @return `seg_step` returns a numeric vector with one entry per
#'   break-point; `seg_converged` returns a single logical.
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
#' @seealso [seg()], [seg_psi()], [seg_start()]
#'   [seg_start()]
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

#' Reset a Break-Point Term's Scaling Schedule
#'
#' @description
#' The term with its scaling factors back at `c0`, the directions and
#' the step record cleared, and the break-points where they are. A restart
#' needs it: the schedule is a state of the iteration that only ever
#' tightens, so an iteration resumed from a converged fit inherits factors
#' at their floor and cannot travel. Measured, bootstrap restarts without
#' the reset returned the incumbent unchanged ten times out of ten.
#'
#' @param term A built break-point term (see [term_build()]).
#'
#' @return The term, ready to iterate afresh from its current positions.
#'
#' @examples
#' set.seed(1)
#' dd <- data.frame(x = sort(runif(100, 0, 10)))
#' b <- term_build(jump(x, psi = 4), dd)
#' b <- term_refresh(b, b@blueprint$coef)
#' identical(seg_reheat(b)@blueprint$cscale, b@blueprint$cscale * 0 + 0.05)
#'
#' @seealso [seg_step()], [term_refresh()]
#' @export
seg_reheat <- function(term) {
  if (!S7::S7_inherits(term, SegTerm)) {
    stop("'term' must be a break-point term.", call. = FALSE)
  }
  .assert_built(term)
  bp <- term@blueprint
  bp$cscale <- rep(term@spec$c0, bp$npsi)
  bp$sgn <- rep(0, bp$npsi)
  bp$step <- rep(NA_real_, bp$npsi)
  bp$nref <- 0L
  term@blueprint <- bp
  term
}

#' Move a Break-Point Term to Given Positions
#'
#' @description
#' The term with its break-points placed at `psi`, the changes kept,
#' the scaling schedule fresh and the block rebuilt, ready to iterate
#' from there. A restart proposal needs it: a bootstrap resample perturbs
#' the objective by \eqn{1/\sqrt{n}} and stops escaping a deep basin as the
#' sample grows, so the restarting loop also proposes fresh positions
#' outright, which this is the operation behind.
#'
#' @details
#' The positions are sorted and confined to the term's own interval. For a
#' discontinuous construction the auxiliary coefficients are set to
#' \eqn{g_k = -\delta_k\psi_k}, a change of level at zero replaced by one,
#' so the read-off returns exactly the positions given. A term whose
#' per-break-point coefficients carry a development is rejected: its
#' positions are one per observation and a single vector does not place
#' them.
#'
#' @param term A built break-point term (see [term_build()]).
#' @param psi A numeric vector, one position per break-point.
#'
#' @return The term at the new positions.
#'
#' @examples
#' set.seed(1)
#' dd <- data.frame(x = sort(runif(100, 0, 10)))
#' b <- term_build(jump(x, psi = 4), dd)
#' seg_psi(seg_relocate(b, 6))
#'
#' @seealso [seg_reheat()], [seg_start()]
#' @export
seg_relocate <- function(term, psi) {
  if (!S7::S7_inherits(term, SegTerm)) {
    stop("'term' must be a break-point term.", call. = FALSE)
  }
  .assert_built(term)
  bp <- term@blueprint
  K <- bp$npsi
  if (!is.numeric(psi) || length(psi) != K || anyNA(psi)) {
    stop(sprintf("'psi' must give %d finite positions.", K), call. = FALSE)
  }
  perk <- c(if (bp$kind %in% c("seg", "jseg")) paste0("gamma", seq_len(K)),
            if (bp$kind %in% c("jump", "jseg")) paste0("delta", seq_len(K)),
            paste0("psi", seq_len(K)))
  if (any(vapply(perk, function(p) !is.null(bp$Z[[p]]), logical(1)))) {
    stop(paste("a term whose per-break-point coefficients carry a",
               "development cannot be relocated by a single vector."),
         call. = FALSE)
  }
  psi <- sort(pmin(pmax(as.numeric(psi), bp$lim[1L]), bp$lim[2L]))
  cf <- bp$coef
  for (k in seq_len(K)) {
    pk <- paste0("psi", k)
    if (bp$kind == "seg" || !is.null(bp$smooth)) {
      cf[bp$index[[pk]]] <- psi[k]
    } else {
      dk <- cf[bp$index[[paste0("delta", k)]]]
      if (abs(dk) < 1e-8) {
        dk <- 1
        cf[bp$index[[paste0("delta", k)]]] <- dk
      }
      cf[bp$index[[pk]]] <- -dk * psi[k]
    }
  }
  n <- length(bp$xv)
  bp$pk <- lapply(seq_len(K), function(k) psi[k])
  bp$psi <- matrix(psi, n, K, byrow = TRUE)
  bp$cscale <- rep(term@spec$c0, K)
  bp$sgn <- rep(0, K)
  bp$step <- rep(NA_real_, K)
  bp$nref <- 0L
  asm <- .seg_assemble(bp, bp$xv, cf)
  X <- asm$X
  colnames(X) <- term@coef_names
  bp$coef <- cf
  bp$value <- asm$value
  bp$psi <- asm$psi
  bp$pk <- asm$pk
  term@X <- X
  term@blueprint <- bp
  term
}

#' @title Whether a Break-Point Term's Break-Points Have Settled
#' @name term_converged.SegTerm
#' @description
#' [seg_converged()], so that a fitting layer reads the rule the
#' construction is stopped on without knowing it is holding a break-point
#' term. A term that has not been refreshed has not moved and reports
#' `FALSE`, the first refresh being the one that measures nothing.
#' @param term A built [SegTerm()].
#' @param ... Unused.
#' @return A single logical.
#' @keywords internal
S7::method(term_converged, SegTerm) <- function(term, ...) {
  isTRUE(seg_converged(term))
}

#' Starting Positions for a Break-Point Term
#'
#' @description
#' Chooses the starting positions of a [seg()], [jump()]
#' or [jseg()] term by scoring an equally spaced grid on the
#' least-squares profile of the term's own columns, and returns the
#' specification with `psi` set to the best combination found.
#'
#' @details
#' \cite{fasola2018} recommend fixing the starting value by evaluating
#' the objective on a small grid spanned over the range of the covariate
#' instead of at a single conventional point, and the recommendation
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
#' over an equally spaced grid \eqn{\mathcal{G}} of `k` points in the
#' range of the covariate. The inner minimization is a linear fit, so the
#' whole rule costs `k` of them.
#'
#' The grid is scored on the residual sum of squares of an intercept,
#' the term's columns at each candidate position and, where the term
#' carries one, the linear effect. That is the exact profile for a
#' gaussian response and an adequate starting rule for any other, the
#' quantity being used to place a starting value, never to fit. The
#' positions found seed a development as well, each starting vector
#' projecting the position onto the sub-design. With several break-points
#' every increasing combination of grid points is scored, so `k`
#' should be kept small.
#'
#' @param spec An unbuilt [SegTerm()].
#' @param data A data frame in which the covariate is evaluated.
#' @param y The response, one value per row of `data`.
#' @param k The number of grid points. Defaults to 10.
#'
#' @return The specification, with `psi` set.
#'
#' @references
#' Fasola, S., Muggeo, V. M. R. and Kuchenhoff, H. (2018). A heuristic,
#' iterative algorithm for change-point detection in abrupt change
#' models. *Computational Statistics*, 33, 997--1015.
#'
#' @examples
#' set.seed(1)
#' dd <- data.frame(x = sort(runif(200, 0, 10)))
#' dd$y <- 0.3 * dd$x + 2 * (dd$x > 6.5) + rnorm(200, sd = 0.3)
#' seg_start(jump(x), dd, dd$y)@spec$psi
#'
#' @seealso [seg()], [jump()], [jseg()]
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

#' Polish a Break-Point Term's Positions on the Exact Profile
#'
#' @description
#' Coordinate descent over the break-point positions on the exact profile:
#' one position at a time is swept over a grid with the others held, the
#' least-squares fit at fixed positions being an ordinary linear model, and
#' the sweeps repeat until no position moves. The term comes back relocated
#' ([seg_relocate()]) at the best positions found.
#'
#' @details
#' This is [seg_start()]'s device made conditional, and it exists
#' for the failure a start cannot fix: with several break-points the
#' iteration can capture two of them on one feature and press the third
#' against its confinement limit, and neither a bootstrap excursion --
#' whose perturbation is of order \eqn{1/\sqrt{n}}, nor a random
#' relocation escapes it, the basin of the right configuration being
#' narrow. A grid sweep of one position with the others held walks
#' straight to the missing feature, because the profile at fixed positions
#' is exact and unimodal around it.
#'
#' The profile is least squares of `y` on the term's own columns at
#' fixed positions, so it is exact for a gaussian response and a starting
#' rule for any other, the argument [seg_start()] already makes.
#' A term whose per-break-point coefficients carry a development is
#' rejected, its positions being one per observation.
#'
#' With `weights`, the profile is weighted least squares. A restarting
#' loop sweeps the profile of a bootstrap resample that way, the
#' multinomial counts of sampling the rows with replacement as weights --
#' which moves the profile's optima the way refitting the resample would,
#' at the cost of grid-many linear fits instead of a whole model fit: the
#' non-convexity of these models lives entirely in the positions, so
#' exploring the positions is the whole of what an excursion is for.
#'
#' @param term A built break-point term (see [term_build()]).
#' @param y A numeric vector, one value per observation of the build data:
#'   the response, net of whatever the caller wants held.
#' @param k How many grid points per sweep. Defaults to 50.
#' @param sweeps At most how many passes over the positions. Defaults
#'   to 10; the descent usually stops moving in two or three.
#' @param weights Optional non-negative weights, one per observation; the
#'   profile is then weighted least squares.
#'
#' @return The term at the polished positions.
#'
#' @examples
#' set.seed(1)
#' dd <- data.frame(x = sort(runif(300, 0, 10)))
#' dd$y <- 2 * (dd$x > 3) - 1.5 * (dd$x > 7) + rnorm(300, sd = 0.3)
#' b <- term_build(jump(x, npsi = 2, psi = c(1, 2)), dd)
#' seg_psi(seg_polish(b, dd$y))
#'
#' @seealso [seg_start()], [seg_relocate()]
#' @export
seg_polish <- function(term, y, k = 50, sweeps = 10, weights = NULL) {
  pr <- .seg_profile(term, y, weights)
  grid <- unique(as.numeric(stats::quantile(
    pr$xv, seq(0.02, 0.98, length.out = as.integer(k)), names = FALSE)))
  grid <- pmin(pmax(grid, pr$lim[1L]), pr$lim[2L])
  psi <- pr$psi
  best <- pr$rss(psi)
  K <- length(psi)
  for (s in seq_len(as.integer(sweeps))) {
    moved <- FALSE
    for (j in seq_len(K)) {
      for (g in grid) {
        cand <- psi
        cand[j] <- g
        v <- pr$rss(cand)
        if (v < best - 1e-10 * (best + 1)) {
          best <- v
          psi <- cand
          moved <- TRUE
        }
      }
    }
    if (!moved) break
  }
  seg_relocate(term, psi)
}

#' The Exact Profile of a Break-Point Term at Its Positions
#'
#' @description
#' The residual sum of squares of `y` on the term's own columns at its
#' current positions, the number [seg_polish()] descends on,
#' read at one point. A restarting loop screens proposals with it: two
#' configurations of positions are compared by their profiles at the cost
#' of two linear fits, where refitting the model to compare them costs a
#' whole fit each.
#'
#' @inheritParams seg_polish
#'
#' @return A single number; `Inf` where the fit fails.
#'
#' @examples
#' set.seed(1)
#' dd <- data.frame(x = sort(runif(100, 0, 10)))
#' dd$y <- 2 * (dd$x > 6) + rnorm(100, sd = 0.3)
#' b <- term_build(jump(x, psi = 4), dd)
#' seg_profile_rss(b, dd$y) > seg_profile_rss(seg_relocate(b, 6), dd$y)
#'
#' @seealso [seg_polish()], [seg_relocate()]
#' @export
seg_profile_rss <- function(term, y, weights = NULL) {
  pr <- .seg_profile(term, y, weights)
  pr$rss(pr$psi)
}

# The shared machinery of the two: validation, the fixed-position columns
# and the (weighted) least-squares evaluator.
.seg_profile <- function(term, y, weights = NULL) {
  if (!S7::S7_inherits(term, SegTerm)) {
    stop("'term' must be a break-point term.", call. = FALSE)
  }
  .assert_built(term)
  bp <- term@blueprint
  K <- bp$npsi
  perk <- c(if (bp$kind %in% c("seg", "jseg")) paste0("gamma", seq_len(K)),
            if (bp$kind %in% c("jump", "jseg")) paste0("delta", seq_len(K)),
            paste0("psi", seq_len(K)))
  if (any(vapply(perk, function(p) !is.null(bp$Z[[p]]), logical(1)))) {
    stop(paste("a term whose per-break-point coefficients carry a",
               "development has no single profile."), call. = FALSE)
  }
  xv <- bp$xv
  y <- as.numeric(y)
  if (length(y) != length(xv)) {
    stop("'y' must have one value per observation of the build data.",
         call. = FALSE)
  }
  sw <- NULL
  if (!is.null(weights)) {
    weights <- as.numeric(weights)
    if (length(weights) != length(xv) || anyNA(weights) ||
        any(weights < 0)) {
      stop("'weights' must be non-negative, one per observation.",
           call. = FALSE)
    }
    sw <- sqrt(weights)
    y <- sw * y
  }
  cols <- function(psi) {
    Z <- if (bp$linear) cbind(1, xv) else matrix(1, length(xv), 1L)
    if (bp$kind %in% c("seg", "jseg")) {
      for (p in psi) Z <- cbind(Z, pmax(xv - p, 0))
    }
    if (bp$kind %in% c("jump", "jseg")) {
      for (p in psi) Z <- cbind(Z, as.numeric(xv > p))
    }
    if (is.null(sw)) Z else Z * sw
  }
  rss <- function(psi) {
    out <- tryCatch(sum(qr.resid(qr(cols(psi)), y)^2),
                    error = function(e) Inf)
    if (is.finite(out)) out else Inf
  }
  list(xv = xv, lim = bp$lim, psi = as.numeric(bp$psi[1L, ]), rss = rss)
}

#' @title Print a Break-Point Term
#' @name print.SegTerm
#'
#' @description
#' Prints the label, the construction and how many break-points the term
#' carries, followed by where they currently sit. A smoothed term adds a line
#' naming its smoother and the transition width, and a term with a developed
#' coefficient says which one.
#'
#' @details
#' The form is
#'
#' ```
#' <SegTerm> 'seg': seg, 1 break-point
#'   at: 6.501
#' ```
#'
#' The positions come from [seg_psi()], so they are the current ones and move
#' with every [term_refresh()]. Where a break-point carries a development
#' there is no single number, a position then being one value per observation,
#' and the range is printed instead.
#'
#' The smoother line reports the width because it changes what the term means:
#' the transition is a real feature of the model, not an implementation
#' detail, and a fit at one width is not a fit at another.
#'
#' @param x A [SegTerm()], built or not.
#' @param ... Unused, and accepted so that the signature matches [print()]'s.
#'
#' @return `x`, invisibly. Called for the lines it writes.
#'
#' @seealso [seg_psi()] for the positions, [seg()], [jump()] and [jseg()].
#'
#' @examples
#' set.seed(1)
#' d <- data.frame(x = sort(runif(120, 0, 10)), id = factor(rep(1:4, each = 30)))
#' d$y <- 1 + 0.5 * d$x + 2 * pmax(d$x - 6, 0) + rnorm(120, sd = 0.4)
#'
#' # A specification, and the same term built at its starting quantile.
#' seg(x, npsi = 2)
#' term_build(seg(x, npsi = 2), d)
#'
#' # A smoothed term names its smoother and its width.
#' term_build(jump(x, smoothed = penalties7::smooth_probit()), d)
#'
#' # A developed break-point has a range instead of a number.
#' term_build(seg(x, psi ~ id), d)
#'
#' @keywords internal
S7::method(print, SegTerm) <- function(x, ...) {
  if (term_is_built(x)) {
    dv <- names(x@subformulas)
    cat(sprintf("<SegTerm> '%s': %s, %d break-point%s%s
", x@label, x@kind,
                x@npsi, if (x@npsi == 1L) "" else "s",
                if (length(dv))
                  sprintf("; %s developed", paste(dv, collapse = ", ")) else ""))
    sm <- x@blueprint$smooth
    if (!is.null(sm)) {
      # the width is the transition's, the bent-cable reading, so it is
      # reported rather than treated as a detail
      cat(sprintf("  smoothed (%s, %s = %s%s)
", sm$sm@smoother_name,
                  sm$sm@width_name, format(sm$width, digits = 3),
                  if (!is.null(sm$w_group)) sprintf(", per group [%s, %s]",
                    format(min(sm$w_group), digits = 3),
                    format(max(sm$w_group), digits = 3)) else ""))
    }
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
#' The start [term_build()] computed: unit changes, and the
#' break-points at the positions `psi` names or at the interior
#' quantiles of the covariate. Zero is degenerate here, never neutral: a
#' discontinuous term reads its break-point off \eqn{-g_k/\delta_k},
#' which at zero is the same clamped position for every one of them, and a
#' continuous term's Jacobian column vanishes, so a fitting layer that
#' starts every coefficient at zero has to be told otherwise.
#' @param term A built [SegTerm()].
#' @param target Unused: a break-point term already reads the covariate's
#'   interior quantiles at [term_build()].
#' @param ... Unused.
#' @return A numeric vector, one value per column of the block.
#' @keywords internal
S7::method(term_coef_start, SegTerm) <- function(term, target = NULL, ...) {
  .assert_built(term)
  term@blueprint$coef
}

#' @title What a Fitted Break-Point Term Is About
#' @name term_readable.SegTerm
#'
#' @description
#' The quantities of the model the term defines, in place of the
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
#' which is the delta method `segmented` reports a break-point's
#' standard error by. The joint construction reads its position from a
#' quadratic that also carries the change of slope, and the two readings
#' differ by the increment the iteration has left: the value reported is
#' the one the term holds, so that this and [seg_psi()] cannot
#' give two numbers for one quantity, and the Jacobian is the fixed
#' point's, which the quadratic degenerates to once the increment has
#' vanished.
#'
#' Every quantity is on the identity scale: a change is unbounded and a
#' break-point is a position on the covariate's own scale, held inside the
#' interval between the 5th and the 95th percentile. Where a coefficient
#' carries a development there is no single number to report, a
#' break-point then having one value per observation, and the method
#' returns nothing, leaving a caller to report the coefficients
#' themselves.
#'
#' @param term A built [SegTerm()].
#' @param zeta The term's coefficients, in the order of its block.
#' @param ... Unused.
#'
#' @return A list with `name`, `value`, `jacobian` and
#'   `scale`, as [term_readable()] documents, or
#'   `NULL` where a coefficient carries a development.
#'
#' @keywords internal
S7::method(term_readable, SegTerm) <- function(term, zeta, ...) {
  .assert_built(term)
  bp <- term@blueprint
  cf <- as.numeric(zeta)
  # QUANTITY BY QUANTITY and not all or nothing. A developed parameter is a
  # vector of coefficients over covariates and has no one value to report,
  # so it is skipped; the parameters beside it are unaffected and are
  # exactly what a reader wants, the working coefficients a discontinuous
  # construction is fitted through being no part of the model. A jump's
  # position is read off its change of level, so it is reported only where
  # both are numbers.
  dev <- function(p) !is.null(bp$subs[[p]])
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

  if (bp$linear && !dev("beta")) push("beta", cf[at("beta")], sel(at("beta")))
  if (bp$kind %in% c("seg", "jseg")) {
    for (k in seq_len(K)) {
      if (dev(paste0("gamma", k))) next
      i <- at(paste0("gamma", k))
      push(paste0("gamma", k), cf[i], sel(i))
    }
  }
  if (bp$kind %in% c("jump", "jseg")) {
    for (k in seq_len(K)) {
      if (dev(paste0("delta", k))) next
      i <- at(paste0("delta", k))
      push(paste0("delta", k), cf[i], sel(i))
    }
  }
  # the position the TERM holds, read by the same function it reads it with,
  # so that this and seg_psi() cannot report two numbers for one quantity
  pos_at <- function(k) {
    v <- .seg_pk(bp, cf, k, bp$pk[[k]])
    min(max(v, bp$lim[1L]), bp$lim[2L])
  }
  for (k in seq_len(K)) {
    if (dev(paste0("psi", k))) next
    ip <- at(paste0("psi", k))
    if (bp$kind == "seg" || !is.null(bp$smooth)) {
      push(paste0("psi", k), pos_at(k), sel(ip))
    } else {
      if (dev(paste0("delta", k))) next
      id <- at(paste0("delta", k))
      d <- cf[id]
      if (abs(d) < 1e-12) d <- sign(d + (d == 0)) * 1e-12
      j <- numeric(np)
      j[ip] <- -1 / d
      j[id] <- cf[ip] / d^2
      push(paste0("psi", k), pos_at(k), j)
    }
  }
  if (!length(nm)) return(NULL)
  J <- do.call(rbind, rows)
  dimnames(J) <- list(nm, term@coef_names)
  list(name = nm, value = val, jacobian = J,
       scale = stats::setNames(
         rep(list(linkfunctions7::identity_link()), length(nm)), nm))
}
