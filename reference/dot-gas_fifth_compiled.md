# The Fifth Order of a Score-Driven Term, Compiled

Builds the groups that `gas_fifth_gen_cpp()` reads and runs it: for each
group its rows and active unknowns, and every value of the recursion
(the level, the loadings and the partial autocorrelations) as design
rows over the active unknowns with its link's derivatives of orders zero
to five per row. The kernel composes the autoregressive coefficients
from the partial autocorrelations by the Levinson-Durbin map, once per
group or, where a partial autocorrelation carries a subformula, once per
row.

## Usage

``` r
.gas_fifth_compiled(
  term,
  eta,
  score,
  psi,
  g,
  seed,
  blocks_data,
  dirs,
  score_values = NULL
)
```

## Arguments

- term:

  A built term.

- eta:

  The static part of the predictor.

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

- dirs:

  The three directions.

- score_values:

  `NULL`, or the driving score at every observation; with `NULL` the
  callback `score` is evaluated at each predictor.

## Value

As
[`term_fifth()`](https://statmodels7.github.io/modelterms7/reference/term_fifth.md).
