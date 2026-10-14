# Smoothed State Probabilities of a Regime Term

A forward pass with a backward pass beside it, both normalized and
computed in R by this method, giving the probability of each regime at
each observation given the whole series. The recursion runs over each
group in time order, from the chain's stationary distribution.

## Arguments

- term:

  A built
  [`RegimeTerm()`](https://statmodels7.github.io/modelterms7/reference/RegimeTerm.md).

- eta:

  The static predictor.

- y:

  The response.

- logdens:

  The log-density, returning one value per observation when given the
  whole index vector; any other length signals an error.

- psi:

  The term's parameters on the parameter scale, named as
  [`term_params()`](https://statmodels7.github.io/modelterms7/reference/term_params.md).

- ...:

  Unused.

## Value

A numeric matrix, one row per observation and one column per regime.
