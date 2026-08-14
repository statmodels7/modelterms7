# Package index

## The classes and the generics

The contract a term implements: how it builds a block from data, how it
reproduces that block on new rows, what it carries in the way of
coefficients, hyperparameters and penalties, and the two shapes a
structural term takes instead.

- [`model_term()`](https://statmodels7.github.io/modelterms7/reference/model_term.md)
  : S7 Base Class for Model Terms
- [`additive_term()`](https://statmodels7.github.io/modelterms7/reference/additive_term.md)
  : S7 Class for Additive Terms
- [`structural_term()`](https://statmodels7.github.io/modelterms7/reference/structural_term.md)
  : S7 Class for Structural Terms
- [`term_build()`](https://statmodels7.github.io/modelterms7/reference/term_build.md)
  : Build a Term on Data
- [`term_is_built()`](https://statmodels7.github.io/modelterms7/reference/term_is_built.md)
  : Whether a Term Has Been Built
- [`term_matrix()`](https://statmodels7.github.io/modelterms7/reference/term_matrix.md)
  : Design Block of a Built Term
- [`term_penalty()`](https://statmodels7.github.io/modelterms7/reference/term_penalty.md)
  : Penalty of a Term
- [`term_penalties()`](https://statmodels7.github.io/modelterms7/reference/term_penalties.md)
  : Every Penalty a Term Carries
- [`term_npar()`](https://statmodels7.github.io/modelterms7/reference/term_npar.md)
  : Number of Parameters of a Built Term
- [`term_coef_start()`](https://statmodels7.github.io/modelterms7/reference/term_coef_start.md)
  : Where a Term's Own Coefficients Begin
- [`term_coef_names()`](https://statmodels7.github.io/modelterms7/reference/term_coef_names.md)
  : Coefficient Names of a Built Term
- [`term_smooth()`](https://statmodels7.github.io/modelterms7/reference/term_smooth.md)
  : Whether a Term's Penalized Objective Is Smooth
- [`term_predict()`](https://statmodels7.github.io/modelterms7/reference/term_predict.md)
  : Design Block on New Data

## The parametric term

The unpenalized block of a one-sided formula, with the usual
model.matrix conventions for factors, contrasts and the intercept.

- [`linpar()`](https://statmodels7.github.io/modelterms7/reference/linpar.md)
  : Unpenalized Parametric Term
- [`LinparTerm()`](https://statmodels7.github.io/modelterms7/reference/LinparTerm.md)
  : S7 Class for the Unpenalized Parametric Term

## The penalized terms

A block with a penalties7 object attached to its coefficients: ridge,
lasso, SCAD, MCP and the elastic net, over a formula or over a matrix
column of the data.

- [`ridge()`](https://statmodels7.github.io/modelterms7/reference/ridge.md)
  [`lasso()`](https://statmodels7.github.io/modelterms7/reference/ridge.md)
  [`enet()`](https://statmodels7.github.io/modelterms7/reference/ridge.md)
  [`scad()`](https://statmodels7.github.io/modelterms7/reference/ridge.md)
  [`mcp()`](https://statmodels7.github.io/modelterms7/reference/ridge.md)
  : Penalized Parametric Terms
- [`PenalizedTerm()`](https://statmodels7.github.io/modelterms7/reference/PenalizedTerm.md)
  : S7 Class for Penalized Parametric Terms

## The random-effect term

Grouped intercepts and slopes, the effects’ distribution attached as the
penalty, which is what a random effect is under penalized likelihood.

- [`random()`](https://statmodels7.github.io/modelterms7/reference/random.md)
  : Grouped Random-Effect Term
- [`RandomTerm()`](https://statmodels7.github.io/modelterms7/reference/RandomTerm.md)
  : S7 Class for Grouped Random-Effect Terms

## Smooth terms

The penalized smooths of one and of several covariates, built on basis7:
a Demmler-Reinsch reparametrization for one variable and a tensor
product, anisotropic by default, for more.

- [`s()`](https://statmodels7.github.io/modelterms7/reference/s.md) :
  Penalized Smooth of One Covariate
- [`te()`](https://statmodels7.github.io/modelterms7/reference/te.md) :
  Penalized Smooth of Several Covariates
- [`SmoothTerm()`](https://statmodels7.github.io/modelterms7/reference/SmoothTerm.md)
  : S7 Class for Smooth Terms

## Nonlinear terms

A contribution nonlinear in its own parameters, whose block is the
Jacobian recomputed as they move, so a linear fit on it is a
Gauss-Newton step.

- [`nl()`](https://statmodels7.github.io/modelterms7/reference/nl.md) :
  Nonlinear Parametric Term
- [`NlTerm()`](https://statmodels7.github.io/modelterms7/reference/NlTerm.md)
  : S7 Class for Nonlinear Parametric Terms
- [`term_refresh()`](https://statmodels7.github.io/modelterms7/reference/term_refresh.md)
  : Refresh a Term at New Coefficients
- [`term_value()`](https://statmodels7.github.io/modelterms7/reference/term_value.md)
  : The Contribution of a Term at Its Current Coefficients
- [`term_converged()`](https://statmodels7.github.io/modelterms7/reference/term_converged.md)
  : Has a Term's Own Iteration Settled?

## Break-point terms

The break-points at which an effect changes slope, level or both, each
written as a working linear model in the break-point position, with the
grid rule that chooses where to start.

- [`seg()`](https://statmodels7.github.io/modelterms7/reference/seg.md)
  : Segmented Term: a Broken Line with Estimated Break-Points
- [`jump()`](https://statmodels7.github.io/modelterms7/reference/jump.md)
  : Stepmented Term: a Level that Changes at Estimated Break-Points
- [`jseg()`](https://statmodels7.github.io/modelterms7/reference/jseg.md)
  : Segmented-with-Jump Term: Slope and Level Both Changing
- [`SegTerm()`](https://statmodels7.github.io/modelterms7/reference/SegTerm.md)
  : S7 Class for Break-Point Terms
- [`seg_psi()`](https://statmodels7.github.io/modelterms7/reference/seg_psi.md)
  : The Break-Points of a Break-Point Term
- [`seg_step()`](https://statmodels7.github.io/modelterms7/reference/seg_step.md)
  [`seg_converged()`](https://statmodels7.github.io/modelterms7/reference/seg_step.md)
  : The Progress of a Break-Point Iteration
- [`seg_start()`](https://statmodels7.github.io/modelterms7/reference/seg_start.md)
  : Starting Positions for a Break-Point Term

## Structural terms

Terms that rewrite the likelihood rather than adding to the predictor:
score-driven dynamics and a latent Markov chain, both propagating their
exact derivative beside the state.

- [`gas()`](https://statmodels7.github.io/modelterms7/reference/gas.md)
  : Score-Driven Dynamics
- [`GasTerm()`](https://statmodels7.github.io/modelterms7/reference/GasTerm.md)
  : S7 Class for Score-Driven Dynamics
- [`regime()`](https://statmodels7.github.io/modelterms7/reference/regime.md)
  : Markov Regime Switching
- [`RegimeTerm()`](https://statmodels7.github.io/modelterms7/reference/RegimeTerm.md)
  : S7 Class for Markov Regime Terms
- [`term_params()`](https://statmodels7.github.io/modelterms7/reference/term_params.md)
  : Parameters of a Structural Term
- [`term_links()`](https://statmodels7.github.io/modelterms7/reference/term_links.md)
  : Links of a Structural Term's Parameters
- [`term_start()`](https://statmodels7.github.io/modelterms7/reference/term_start.md)
  : Where a Term's Own Parameters Start
- [`term_level_design()`](https://statmodels7.github.io/modelterms7/reference/term_level_design.md)
  : The Design of a Term's Level Development
- [`term_filter()`](https://statmodels7.github.io/modelterms7/reference/term_filter.md)
  : Apply a Structural Term to a Linear Predictor
- [`term_loglik()`](https://statmodels7.github.io/modelterms7/reference/term_loglik.md)
  : Log-Likelihood Contributions of a Structural Term
- [`term_adjoint()`](https://statmodels7.github.io/modelterms7/reference/term_adjoint.md)
  : Differentiate a Structural Term Backwards
- [`term_curvature()`](https://statmodels7.github.io/modelterms7/reference/term_curvature.md)
  : Second Derivatives of a Structural Term's Predictor
- [`term_third()`](https://statmodels7.github.io/modelterms7/reference/term_third.md)
  : Third Derivatives of a Structural Term's Predictor, in One Direction
- [`term_hessian()`](https://statmodels7.github.io/modelterms7/reference/term_hessian.md)
  : The Observed Hessian of a Likelihood Mixed Over States
- [`term_posterior()`](https://statmodels7.github.io/modelterms7/reference/term_posterior.md)
  : The Smoothed State Probabilities of a Latent Markov Term
- [`term_level_param()`](https://statmodels7.github.io/modelterms7/reference/term_level_param.md)
  : Which of a Term's Parameters Acts as an Intercept
- [`term_readable()`](https://statmodels7.github.io/modelterms7/reference/term_readable.md)
  : The Quantities a Fitted Term Reports

## The formula layer

Reading a model formula into terms, by what each call evaluates to
rather than by its name, so a term class defined outside the package
needs no registration.

- [`interpret_formula()`](https://statmodels7.github.io/modelterms7/reference/interpret_formula.md)
  : Interpret a Model Formula Into Terms

## The response

Marking a censored response on the left side of a formula, the status of
each observation derived from the values given.

- [`cens()`](https://statmodels7.github.io/modelterms7/reference/cens.md)
  : Censored Response Constructor
- [`censored_response()`](https://statmodels7.github.io/modelterms7/reference/censored_response.md)
  : S7 Class for a Censored Response

## Fitted views

What a fitted term reports: effective degrees of freedom counted per
penalty, the printed summary, and the plot at supplied coefficients.

- [`edf()`](https://statmodels7.github.io/modelterms7/reference/edf.md)
  : Effective Degrees of Freedom of a Term

## Validation

The structural checks a term must pass against a data frame, including
the one a term that rebuilds instead of reapplying its blueprint fails.

- [`check_term()`](https://statmodels7.github.io/modelterms7/reference/check_term.md)
  : Numerical Validation of a Model Term
