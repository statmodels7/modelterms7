# Package index

## The classes and the generics

The two branches of the class tree and the generics that a term
implements: building a block from data, reproducing that block on new
rows, and reporting coefficients, hyperparameters and penalties. A
structural term implements a filter or a likelihood instead of a block.

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
  : Every Penalty That a Term Carries
- [`term_components()`](https://statmodels7.github.io/modelterms7/reference/term_components.md)
  : How a Term's Columns Divide Among Its Own Parameters
- [`term_hyper()`](https://statmodels7.github.io/modelterms7/reference/term_hyper.md)
  : The Hyperparameters That a Term Holds
- [`term_grid()`](https://statmodels7.github.io/modelterms7/reference/term_grid.md)
  : The Grid Sizes of a Term
- [`term_values()`](https://statmodels7.github.io/modelterms7/reference/term_values.md)
  : The Hyperparameter Values That a Term Writes Out
- [`term_path_min()`](https://statmodels7.github.io/modelterms7/reference/term_path_min.md)
  : How Far Down Its Path a Term Reaches
- [`term_search()`](https://statmodels7.github.io/modelterms7/reference/term_search.md)
  : How a Term Covers Its Own Kinked Hyperparameters
- [`term_ids()`](https://statmodels7.github.io/modelterms7/reference/term_ids.md)
  : Which of a Term's Hyperparameters Are Shared
- [`term_npar()`](https://statmodels7.github.io/modelterms7/reference/term_npar.md)
  : Number of Parameters of a Built Term
- [`term_coef_start()`](https://statmodels7.github.io/modelterms7/reference/term_coef_start.md)
  : Where a Term's Own Coefficients Begin
- [`term_coef_draw()`](https://statmodels7.github.io/modelterms7/reference/term_coef_draw.md)
  : Coefficients That a Term Draws for Itself
- [`term_jacobian_block()`](https://statmodels7.github.io/modelterms7/reference/term_jacobian_block.md)
  : Is a Term's Block the Jacobian of Its Contribution?
- [`term_kinks()`](https://statmodels7.github.io/modelterms7/reference/term_kinks.md)
  : Where a Term's Objective Has a Kink
- [`term_kink_constraints()`](https://statmodels7.github.io/modelterms7/reference/term_kink_constraints.md)
  : The Linear Constraints a Kink Places on a Term's Coefficients
- [`term_held()`](https://statmodels7.github.io/modelterms7/reference/term_held.md)
  : Coefficients Held by a Term
- [`term_tag()`](https://statmodels7.github.io/modelterms7/reference/term_tag.md)
  : The Covariance Label of a Term
- [`term_group()`](https://statmodels7.github.io/modelterms7/reference/term_group.md)
  : The Grouping by Which a Term's Coefficients Are Indexed
- [`term_within()`](https://statmodels7.github.io/modelterms7/reference/term_within.md)
  : The Within-Group Design of a Term With Grouped Effects
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
  : Ridge Penalty on a Block of Coefficients
- [`lasso()`](https://statmodels7.github.io/modelterms7/reference/lasso.md)
  : Lasso Penalty on a Block of Coefficients
- [`enet()`](https://statmodels7.github.io/modelterms7/reference/enet.md)
  : Elastic Net Penalty on a Block of Coefficients
- [`scad()`](https://statmodels7.github.io/modelterms7/reference/scad.md)
  : SCAD Penalty on a Block of Coefficients
- [`mcp()`](https://statmodels7.github.io/modelterms7/reference/mcp.md)
  : MCP Penalty on a Block of Coefficients
- [`penalized_terms`](https://statmodels7.github.io/modelterms7/reference/penalized_terms.md)
  : What the Penalized Terms Share
- [`PenalizedTerm()`](https://statmodels7.github.io/modelterms7/reference/PenalizedTerm.md)
  : S7 Class for Penalized Parametric Terms

## The random-effect term

Grouped intercepts, slopes and smooth curves, in which the distribution
of the effects is attached as the penalty. This is how a random effect
enters a penalized likelihood.

- [`random()`](https://statmodels7.github.io/modelterms7/reference/random.md)
  : Grouped Random-Effect Term
- [`RandomTerm()`](https://statmodels7.github.io/modelterms7/reference/RandomTerm.md)
  : S7 Class for Grouped Random-Effect Terms
- [`RandomSmoothTerm()`](https://statmodels7.github.io/modelterms7/reference/RandomSmoothTerm.md)
  : S7 Class for a Smooth Random Effect

## Smooth terms

The penalized smooths of one and of several covariates, built on the
smoothers of basis7: s() for one covariate and te() for a tensor
product, anisotropic by default.

- [`s()`](https://statmodels7.github.io/modelterms7/reference/s.md) :
  Penalized Smooth of One Covariate
- [`te()`](https://statmodels7.github.io/modelterms7/reference/te.md) :
  Penalized Smooth of Several Covariates
- [`SmoothTerm()`](https://statmodels7.github.io/modelterms7/reference/SmoothTerm.md)
  : S7 Class for Smooth Terms

## Shape-constrained smooths

A smooth of one covariate constrained to be monotone, convex or concave,
through a map from free coefficients to the B-spline coefficients, with
a penalty on the free coefficients.

- [`scs()`](https://statmodels7.github.io/modelterms7/reference/scs.md)
  : Shape-Constrained Smooth of One Covariate
- [`ScsTerm()`](https://statmodels7.github.io/modelterms7/reference/ScsTerm.md)
  : S7 Class for Shape-Constrained Smooth Terms
- [`shape_increasing()`](https://statmodels7.github.io/modelterms7/reference/shape_constructors.md)
  [`shape_decreasing()`](https://statmodels7.github.io/modelterms7/reference/shape_constructors.md)
  [`shape_convex()`](https://statmodels7.github.io/modelterms7/reference/shape_constructors.md)
  [`shape_concave()`](https://statmodels7.github.io/modelterms7/reference/shape_constructors.md)
  : Shape Constraints for a Smooth
- [`shape_constraint()`](https://statmodels7.github.io/modelterms7/reference/shape_constraint.md)
  [`ShapeMonotone()`](https://statmodels7.github.io/modelterms7/reference/shape_constraint.md)
  [`ShapeCurvature()`](https://statmodels7.github.io/modelterms7/reference/shape_constraint.md)
  : S7 Class for Shape Constraints
- [`shape_map()`](https://statmodels7.github.io/modelterms7/reference/shape_map.md)
  : The Map of a Shape Constraint
- [`shape_start()`](https://statmodels7.github.io/modelterms7/reference/shape_start.md)
  : The Starting Value of a Shape Constraint

## Nonlinear terms

A contribution nonlinear in its own parameters, whose block is the
Jacobian recomputed as they move, so a linear fit on it is a
Gauss-Newton step.

- [`nl()`](https://statmodels7.github.io/modelterms7/reference/nl.md) :
  Nonlinear Parametric Term
- [`nl_fderiv()`](https://statmodels7.github.io/modelterms7/reference/nl_fderiv.md)
  : The Derivatives of a Nonlinear Term in Its Own Parameters
- [`NlTerm()`](https://statmodels7.github.io/modelterms7/reference/NlTerm.md)
  : S7 Class for Nonlinear Parametric Terms
- [`term_refresh()`](https://statmodels7.github.io/modelterms7/reference/term_refresh.md)
  : Refresh a Term at New Coefficients
- [`term_block_contract()`](https://statmodels7.github.io/modelterms7/reference/term_block_contract.md)
  : How a Term's Block Moves With Its Coefficients
- [`term_block_deriv()`](https://statmodels7.github.io/modelterms7/reference/term_block_deriv.md)
  : How a Term's Block Moves Along One Direction
- [`term_block_deriv2()`](https://statmodels7.github.io/modelterms7/reference/term_block_deriv2.md)
  : The Second Derivative of a Design Block in the Coefficients
- [`term_value()`](https://statmodels7.github.io/modelterms7/reference/term_value.md)
  : The Contribution of a Term at Its Current Coefficients
- [`term_converged()`](https://statmodels7.github.io/modelterms7/reference/term_converged.md)
  : Has a Term's Own Iteration Settled?
- [`term_stalled()`](https://statmodels7.github.io/modelterms7/reference/term_stalled.md)
  : Has a Term's Own Iteration Run Out of Step Control?

## Break-point terms

The break-points at which an effect changes slope, level or both. The
sharp constructions are estimated through a working linear model in the
break-point position, a smoothed construction treats each break-point as
an ordinary parameter, and the grid rule chooses where to start.

- [`seg()`](https://statmodels7.github.io/modelterms7/reference/seg.md)
  : Segmented Term: a Broken Line with Estimated Break-Points
- [`jump()`](https://statmodels7.github.io/modelterms7/reference/jump.md)
  : Step Term: a Level that Changes at Estimated Break-Points
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
- [`seg_reheat()`](https://statmodels7.github.io/modelterms7/reference/seg_reheat.md)
  : Reset a Break-Point Term's Scaling Schedule
- [`seg_relocate()`](https://statmodels7.github.io/modelterms7/reference/seg_relocate.md)
  : Move a Break-Point Term to Given Positions
- [`seg_polish()`](https://statmodels7.github.io/modelterms7/reference/seg_polish.md)
  : Polish a Break-Point Term's Positions on the Exact Profile
- [`seg_profile_rss()`](https://statmodels7.github.io/modelterms7/reference/seg_profile_rss.md)
  : The Exact Profile of a Break-Point Term at Its Positions
- [`seg_polish_exact()`](https://statmodels7.github.io/modelterms7/reference/seg_polish_exact.md)
  : Polish a Break-Point Term's Positions Over Every Interval
- [`seg_profile_intervals()`](https://statmodels7.github.io/modelterms7/reference/seg_profile_intervals.md)
  : The Profile of One Break-Point Over Every Interval
- [`seg_hold()`](https://statmodels7.github.io/modelterms7/reference/seg_hold.md)
  : Hold a Break-Point Term at Its Positions
- [`seg_pin()`](https://statmodels7.github.io/modelterms7/reference/seg_pin.md)
  : Confine the Free Break-Points of a Term Away From the Held Ones
- [`seg_limits()`](https://statmodels7.github.io/modelterms7/reference/seg_limits.md)
  : Confinement Intervals of the Break-Points of a Term

## Structural terms

Terms that rewrite the likelihood instead of adding to the predictor:
score-driven dynamics, a latent Markov chain and break-points integrated
out of the likelihood. Each of these terms carries its exact derivatives
through its recursion.

- [`gas()`](https://statmodels7.github.io/modelterms7/reference/gas.md)
  : Score-Driven Dynamics
- [`GasTerm()`](https://statmodels7.github.io/modelterms7/reference/GasTerm.md)
  : S7 Class for Score-Driven Dynamics
- [`gas_joint()`](https://statmodels7.github.io/modelterms7/reference/gas_joint.md)
  : Filter Several Score-Driven Terms as One Recursion
- [`GasJoint()`](https://statmodels7.github.io/modelterms7/reference/GasJoint.md)
  : A Joint Score-Driven Filter
- [`GasDriver()`](https://statmodels7.github.io/modelterms7/reference/GasDriver.md)
  : A Static Parameter Whose Score Drives a Joint Filter
- [`regime()`](https://statmodels7.github.io/modelterms7/reference/regime.md)
  : Markov Regime Switching
- [`RegimeTerm()`](https://statmodels7.github.io/modelterms7/reference/RegimeTerm.md)
  : S7 Class for Markov Regime Terms
- [`MarginalBreakTerm()`](https://statmodels7.github.io/modelterms7/reference/MarginalBreakTerm.md)
  : S7 Class for Marginal Break-Point Terms
- [`term_params()`](https://statmodels7.github.io/modelterms7/reference/term_params.md)
  : Parameters of a Structural Term
- [`term_links()`](https://statmodels7.github.io/modelterms7/reference/term_links.md)
  : Links of a Structural Term's Parameters
- [`term_charted()`](https://statmodels7.github.io/modelterms7/reference/term_charted.md)
  : Charted Coordinates of a Term
- [`term_start()`](https://statmodels7.github.io/modelterms7/reference/term_start.md)
  : Where a Term's Own Parameters Start
- [`term_starts()`](https://statmodels7.github.io/modelterms7/reference/term_starts.md)
  : Where a Fit of a Term's Own Parameters Starts, Several Times
- [`term_draw()`](https://statmodels7.github.io/modelterms7/reference/term_draw.md)
  : Drawing a Term's Own Parameters
- [`term_level_design()`](https://statmodels7.github.io/modelterms7/reference/term_level_design.md)
  : The Design of a Term's Level Development
- [`term_filter()`](https://statmodels7.github.io/modelterms7/reference/term_filter.md)
  : Apply a Structural Term to a Linear Predictor
- [`term_static_deriv()`](https://statmodels7.github.io/modelterms7/reference/term_static_deriv.md)
  : The Derivative of a Filtered Predictor in the Static One
- [`term_continue()`](https://statmodels7.github.io/modelterms7/reference/term_continue.md)
  : Continuing a Structural Term Past the Observed Series
- [`term_simulate()`](https://statmodels7.github.io/modelterms7/reference/term_simulate.md)
  : Drawing a Response From a Structural Term
- [`term_loglik()`](https://statmodels7.github.io/modelterms7/reference/term_loglik.md)
  : Log-Likelihood Contributions of a Structural Term
- [`term_adjoint()`](https://statmodels7.github.io/modelterms7/reference/term_adjoint.md)
  : Differentiate a Structural Term Backwards
- [`term_curvature()`](https://statmodels7.github.io/modelterms7/reference/term_curvature.md)
  : Second Derivatives of a Structural Term's Predictor
- [`term_third()`](https://statmodels7.github.io/modelterms7/reference/term_third.md)
  : Third Derivatives of a Structural Term's Predictor, in One Direction
- [`term_fourth()`](https://statmodels7.github.io/modelterms7/reference/term_fourth.md)
  : Fourth Derivatives of a Structural Term's Predictor, in Two
  Directions
- [`term_fifth()`](https://statmodels7.github.io/modelterms7/reference/term_fifth.md)
  : Fifth Derivatives of a Structural Term's Predictor, in Three
  Directions
- [`term_hessian()`](https://statmodels7.github.io/modelterms7/reference/term_hessian.md)
  : The Observed Hessian of a Likelihood Mixed Over States
- [`term_scores()`](https://statmodels7.github.io/modelterms7/reference/term_scores.md)
  : Per-Observation Scores of a Structural Term
- [`term_posterior()`](https://statmodels7.github.io/modelterms7/reference/term_posterior.md)
  : The Smoothed State Probabilities of a Latent Markov Term
- [`term_shift_jacobian()`](https://statmodels7.github.io/modelterms7/reference/term_shift_jacobian.md)
  : The Posterior Shift of a Term Mixed Over States, and Its Jacobian
- [`term_levels()`](https://statmodels7.github.io/modelterms7/reference/term_levels.md)
  : The Levels of a Likelihood-Shaped Structural Term
- [`term_latent()`](https://statmodels7.github.io/modelterms7/reference/term_latent.md)
  : The Posterior of a Structural Term's Latent Variable
- [`term_level_param()`](https://statmodels7.github.io/modelterms7/reference/term_level_param.md)
  : Which of a Term's Parameters Acts as an Intercept
- [`term_readable()`](https://statmodels7.github.io/modelterms7/reference/term_readable.md)
  : The Quantities That a Fitted Term Reports

## The formula layer

Reading a model formula into terms. A call is a term when its value
inherits from model_term, whatever the name of the function, so a term
class defined outside the package needs no registration.

- [`interpret_formula()`](https://statmodels7.github.io/modelterms7/reference/interpret_formula.md)
  : Interpret a Model Formula Into Terms

## The response

Marking a censored response on the left side of a formula and deriving
the status of each observation from the values that are given.
statmodels7::statmod() rejects a censored response.

- [`cens()`](https://statmodels7.github.io/modelterms7/reference/cens.md)
  : Mark a Response as Censored
- [`censored_response()`](https://statmodels7.github.io/modelterms7/reference/censored_response.md)
  : S7 Class for a Censored Response

## Fitted views

The effective degrees of freedom of a term, computed from the
coefficients, the curvature and the hyperparameters of a fit: one for
each unpenalized parameter, a count of the nonzero coefficients under a
kinked penalty, and one trace over the remaining parameters.

- [`edf()`](https://statmodels7.github.io/modelterms7/reference/edf.md)
  : Effective Degrees of Freedom of a Term

## Validation

The structural checks of an additive term against a data frame, among
them the subset check that a term fails when it rebuilds its block
instead of reapplying its `build_info`.

- [`check_term()`](https://statmodels7.github.io/modelterms7/reference/check_term.md)
  : Structural Checks on a Model Term
