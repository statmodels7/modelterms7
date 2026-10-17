# Log-Likelihood of a Regime Term

Runs the forward recursion over each group in time order, normalizing at
every step, and returns the per-observation contributions with their
exact derivatives.

## Arguments

- term:

  A built `RegimeTerm`.

- eta:

  The static predictor.

- y:

  The response, reaching the recursion through the two callbacks.

- logdens, score:

  The log-density and its derivative in the predictor, as functions that
  each return one value per observation.

- psi:

  The parameters on the parameter scale, named as
  [`term_params()`](https://statmodels7.github.io/modelterms7/reference/term_params.md).

- ...:

  Unused.

## Value

A list with `loglik` and `jacobian`.

## Details

The density and the score of every observation under every regime are
computed first, by one call of each callback per regime with the whole
index vector, and the recursion and its first derivatives are then run
by the compiled kernel `regime_forward_cpp()`. The derivative in a gap
is taken in the gap itself. A gap that is not positive and a `psi` that
lacks a parameter signal an error, as does a callback that does not
return one value per observation.
