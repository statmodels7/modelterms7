# The Smoothed States of a Regime Term

The probability of each regime at each observation given the whole
series, \\P(S_t = j \mid y_1, \dots, y_n)\\, as a data frame.

## Arguments

- term:

  A built
  [`RegimeTerm()`](https://statmodels7.github.io/modelterms7/reference/RegimeTerm.md).

- eta:

  The static predictor.

- y:

  The response.

- logdens:

  The log-density, as
  [`term_loglik()`](https://statmodels7.github.io/modelterms7/reference/term_loglik.md)
  takes it.

- psi:

  The term's parameters on the parameter scale.

- ...:

  Unused.

## Value

A data frame with one row per observation and one column per regime,
`state1`, `state2`, and so on; every row sums to one.

## Details

The probabilities are
[`term_posterior()`](https://statmodels7.github.io/modelterms7/reference/term_posterior.md)'s,
computed by the forward and backward recursions; this method only labels
them. The rows are in the order of the data the term was built on,
whatever `by` and `time` say, so the result can be bound to that data
directly.
