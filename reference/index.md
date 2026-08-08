# Package index

## The classes and the generics

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
- [`term_npar()`](https://statmodels7.github.io/modelterms7/reference/term_npar.md)
  : Number of Coefficients of a Built Term
- [`term_coef_names()`](https://statmodels7.github.io/modelterms7/reference/term_coef_names.md)
  : Coefficient Names of a Built Term
- [`term_smooth()`](https://statmodels7.github.io/modelterms7/reference/term_smooth.md)
  : Whether a Term's Penalized Objective Is Smooth
- [`term_predict()`](https://statmodels7.github.io/modelterms7/reference/term_predict.md)
  : Design Block on New Data

## The parametric term

- [`linpar()`](https://statmodels7.github.io/modelterms7/reference/linpar.md)
  : Unpenalized Parametric Term
- [`LinparTerm()`](https://statmodels7.github.io/modelterms7/reference/LinparTerm.md)
  : S7 Class for the Unpenalized Parametric Term

## The penalized terms

- [`ridge()`](https://statmodels7.github.io/modelterms7/reference/ridge.md)
  [`lasso()`](https://statmodels7.github.io/modelterms7/reference/ridge.md)
  [`scad()`](https://statmodels7.github.io/modelterms7/reference/ridge.md)
  [`mcp()`](https://statmodels7.github.io/modelterms7/reference/ridge.md)
  : Penalized Parametric Terms
- [`PenalizedTerm()`](https://statmodels7.github.io/modelterms7/reference/PenalizedTerm.md)
  : S7 Class for Penalized Parametric Terms

## The random-effect term

- [`random()`](https://statmodels7.github.io/modelterms7/reference/random.md)
  : Grouped Random-Effect Term
- [`RandomTerm()`](https://statmodels7.github.io/modelterms7/reference/RandomTerm.md)
  : S7 Class for Grouped Random-Effect Terms

## Smooth terms

- [`s()`](https://statmodels7.github.io/modelterms7/reference/s.md) :
  Penalized Smooth of One Covariate
- [`te()`](https://statmodels7.github.io/modelterms7/reference/te.md) :
  Penalized Smooth of Several Covariates
- [`SmoothTerm()`](https://statmodels7.github.io/modelterms7/reference/SmoothTerm.md)
  : S7 Class for Smooth Terms

## Nonlinear terms

- [`nl()`](https://statmodels7.github.io/modelterms7/reference/nl.md) :
  Nonlinear Parametric Term
- [`NlTerm()`](https://statmodels7.github.io/modelterms7/reference/NlTerm.md)
  : S7 Class for Nonlinear Parametric Terms
- [`term_refresh()`](https://statmodels7.github.io/modelterms7/reference/term_refresh.md)
  : Refresh a Term at New Coefficients
- [`term_value()`](https://statmodels7.github.io/modelterms7/reference/term_value.md)
  : The Contribution of a Term at Its Current Coefficients

## Structural terms

- [`gas()`](https://statmodels7.github.io/modelterms7/reference/gas.md)
  : Score-Driven Dynamics
- [`GasTerm()`](https://statmodels7.github.io/modelterms7/reference/GasTerm.md)
  : S7 Class for Score-Driven Dynamics
- [`term_params()`](https://statmodels7.github.io/modelterms7/reference/term_params.md)
  : Parameters of a Structural Term
- [`term_links()`](https://statmodels7.github.io/modelterms7/reference/term_links.md)
  : Links of a Structural Term's Parameters
- [`term_filter()`](https://statmodels7.github.io/modelterms7/reference/term_filter.md)
  : Apply a Structural Term to a Linear Predictor

## The formula layer

- [`interpret_formula()`](https://statmodels7.github.io/modelterms7/reference/interpret_formula.md)
  : Interpret a Model Formula Into Terms

## The response

- [`cens()`](https://statmodels7.github.io/modelterms7/reference/cens.md)
  : Censored Response Constructor
- [`censored_response()`](https://statmodels7.github.io/modelterms7/reference/censored_response.md)
  : S7 Class for a Censored Response

## Fitted views

- [`edf()`](https://statmodels7.github.io/modelterms7/reference/edf.md)
  : Effective Degrees of Freedom of a Term

## Validation

- [`check_term()`](https://statmodels7.github.io/modelterms7/reference/check_term.md)
  : Numerical Validation of a Model Term
