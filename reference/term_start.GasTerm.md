# Where a Score-Driven Term's Parameters Start

Zero on the unconstrained scale of every parameter except the score
loadings, which start at \\0.1\\ on the parameter scale, through
whatever link each one carries.

## Arguments

- term:

  A
  [`GasTerm`](https://statmodels7.github.io/modelterms7/reference/GasTerm.md).

- ...:

  Unused.

## Value

A named numeric vector on the unconstrained scale.

## Details

Zero is the natural point of every other chart – a level of zero, no
persistence, no deviation – so the term starts as near the model without
it as its charts allow. The loadings are the exception because zero on
the log scale is a loading of ONE, a response strong enough to
destabilize the recursion at ordinary curvatures; \\0.1\\ is a weak
response, and it is applied on the parameter scale so the start means
the same thing whatever chart a loading rides.
