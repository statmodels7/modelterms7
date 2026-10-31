# The Score-Driven Recursion's Derivatives to the Fifth Order, R Route

Propagates the jets of the filtered predictor over three directions
through the recursion and returns the contractions of
[`term_fifth()`](https://statmodels7.github.io/modelterms7/reference/term_fifth.md).
`.gas_fifth_scalar()` serves a term without subformulas, whose chart is
the same at every observation, and `.gas_fifth_sub()` a term with
subformulas, whose chart is read per observation; both run
`.gas_fifth_run()`.

## Usage

``` r
.gas_fifth_run(term, eta, score, g, seed, blocks_data, dirs, chart_at)

.gas_start_jet(ch, m, dirs)

.gas_fifth_scalar(term, eta, y, score, psi, g, seed, blocks_data, dirs)

.gas_fifth_sub(term, eta, y, score, psi, g, seed, blocks_data, dirs)

.gas_sub_chart_fun(tm, psiv, cols, m, dirs)
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

- g:

  The weights against which the fifth derivative is contracted, one per
  observation.

- seed:

  The derivative of the static predictor in the caller's unknowns.

- blocks_data:

  The model's derivative arrays, as set out above.

- dirs:

  The three directions.

- chart_at:

  A function of a row returning the jets `omega`, `a` and `b` there.

- ch:

  The chart's jets at the group's first observation.

- m:

  The number of unknowns.

- y:

  The response.

- psi:

  The term's parameters, named as
  [`term_params()`](https://statmodels7.github.io/modelterms7/reference/term_params.md).

- tm:

  A built score-driven term, with or without subformulas.

- psiv:

  Its parameters in the order of
  [`term_params()`](https://statmodels7.github.io/modelterms7/reference/term_params.md):
  the parameter scale for a parameter without a subformula, the
  coefficients for a developed one.

- cols:

  The columns of those parameters among the unknowns.

## Value

`.gas_start_jet()` the jet of the starting level; the others as
[`term_fifth()`](https://statmodels7.github.io/modelterms7/reference/term_fifth.md).

## Details

At each observation the level is \\f_t = \omega + \sum_i a_i s\_{t-i} +
\sum_j b_j f\_{t-j}\\, and its jet is the sum of the products of jets
([`.gas_jet_times()`](https://statmodels7.github.io/modelterms7/reference/dot-gas_jet_times.md)).
The score's jet is composed from the family's derivative arrays by
[`.gas_score_jet()`](https://statmodels7.github.io/modelterms7/reference/dot-gas_compose.md).
The orders below the fifth are carried along every subset of the
directions, so this route is a second, independent computation of what
[`term_curvature()`](https://statmodels7.github.io/modelterms7/reference/term_curvature.md),
[`term_third()`](https://statmodels7.github.io/modelterms7/reference/term_third.md)
and
[`term_fourth()`](https://statmodels7.github.io/modelterms7/reference/term_fourth.md)
return. The starting level of a group is read off the fixed point at the
group's first observation, as at the lower orders.

`.gas_sub_chart_fun()` returns the chart at a row for a term whose
values are read per observation: each base parameter and each cross
loading is \\h(z^\top u)\\, with \\z\\ its design row among the unknowns
(a unit vector for a parameter without a subformula), and the
autoregressive coefficients come from
[`.gas_map_jets()`](https://statmodels7.github.io/modelterms7/reference/dot-gas_map_jets.md).
The result at a row carries `x`, the cross loadings' jets in the order
of the term's pairs (partner, lag).
