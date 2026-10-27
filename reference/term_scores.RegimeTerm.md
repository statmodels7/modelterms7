# Per-Observation Scores of a Regime Term

Carries the first derivative, in forward mode, through the scaled
forward recursion and returns the derivative of each weighted increment.

## Arguments

- term:

  A built
  [`RegimeTerm()`](https://statmodels7.github.io/modelterms7/reference/RegimeTerm.md).

- eta:

  The static predictor.

- y:

  The response.

- logdens, grad:

  The log-density and its first derivatives in the predictors.

- psi:

  The term's parameters on the parameter scale.

- seed:

  The derivative of each predictor in the caller's unknowns.

- cols:

  The columns that the term's own parameters occupy.

- level:

  The distribution parameter that the regimes shift.

- weights:

  Observation weights: one number or one per observation, without
  missing values.

- ...:

  Unused.

## Value

A matrix with `n` rows and one column per unknown of the caller.

## Details

The recursion and the chart of the chain are those of
[`term_hessian()`](https://statmodels7.github.io/modelterms7/reference/term_hessian.md)
for a
[`RegimeTerm()`](https://statmodels7.github.io/modelterms7/reference/RegimeTerm.md),
with the derivative of order two omitted. The increment of an
observation is \\\log c_t\\, with \\c_t\\ the normalizing constant of
its step, and its derivative is \\\dot c_t / c_t\\. The callbacks are
called once per regime with the whole index vector.
