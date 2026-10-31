# Fifth Derivatives of a Structural Term's Predictor, in Three Directions

Returns
[`term_fourth()`](https://statmodels7.github.io/modelterms7/reference/term_fourth.md)
differentiated once more along a third direction, with the contractions
of every lower order along the three directions.

## Usage

``` r
term_fifth(term, eta, y, score, psi, g, seed, blocks_data, directions, ...)
```

## Arguments

- term:

  A built term.

- eta:

  The static part of the predictor.

- y:

  The response.

- score:

  The score callback of
  [`term_filter()`](https://statmodels7.github.io/modelterms7/reference/term_filter.md),
  which gives the value of the driving score.

- psi:

  The term's parameters, named as
  [`term_params()`](https://statmodels7.github.io/modelterms7/reference/term_params.md).

- g:

  The weights against which the fifth derivative is contracted, one per
  observation.

- seed:

  The derivative of the static predictor in the caller's unknowns.

- blocks_data:

  The model's derivative arrays, as set out above.

- directions:

  A list of three numeric vectors, the directions \\v\\, \\w\\ and \\x\\
  in the caller's unknowns. The result is symmetric in the three.

- ...:

  Passed to methods.

## Value

A list with `jacobian`, the derivative of the predictor in the caller's
unknowns; `dphi`, a list of three, the second derivative contracted
against each direction; `dpsi`, a list of three, the third derivative
contracted against the pairs \\(v,w)\\, \\(v,x)\\ and \\(w,x)\\;
`domega`, the fourth derivative contracted against the three directions;
each of these with one row per observation; and `curvature`, the fifth
derivative contracted against `g` and the three directions.

## Details

The fifth derivative of the predictor is contracted against three
directions \\v\\, \\w\\ and \\x\\ and the caller's weights, which leaves
an \\m \times m\\ matrix, as at the orders below. The full derivative is
never formed: the recursion propagates one matrix per subset of the
three directions.

Each order of differentiation of the predictor involves one more order
of the family, so this order reads the family's sixth derivative. Unlike
the orders below, the model's derivatives are given as arrays rather
than as a callback, because the contractions that the fifth order forms
from them are many: `blocks_data` is a list with `H`, `D3`, `D4`, `D5`
and `D6`, the derivatives of the driving score in the predictors of
every equation, of orders one to five, one row per observation (`H` one
column per parameter, `D3` one column per pair of parameters at
`(a - 1) * npar + b`, and so on, as for the compiled routes of
[`term_fourth()`](https://statmodels7.github.io/modelterms7/reference/term_fourth.md));
`Vs`, a list of the Jacobian rows of every equation, one `n` by `m`
matrix per parameter (the filter's own entry is not read); and `ap`, the
index of the filter's own parameter.

The base method returns zeros, which is correct for a term whose
predictor is a block of columns times coefficients. A structural term
signals an error, for the reason that
[`term_third()`](https://statmodels7.github.io/modelterms7/reference/term_third.md)
gives.

## See also

[`term_fourth()`](https://statmodels7.github.io/modelterms7/reference/term_fourth.md)
for the order below.

## Examples

``` r
set.seed(1)
dd <- data.frame(t = 1:20, y = rnorm(20))
term <- term_build(linpar(~t), dd)
# an additive term bends no predictor, so every order above the first is
# zero, as the base method returns
out <- term_fifth(term, rep(0, 20), dd$y,
                  score = function(e, i) dd$y[i] - e, psi = list(),
                  g = rep(1, 20), seed = matrix(0, 20, 2),
                  blocks_data = NULL,
                  directions = list(c(1, 0), c(0, 1), c(1, 1)))
all(out$curvature == 0)
#> [1] TRUE
```
