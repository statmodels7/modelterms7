# The Score-Driven Recursion's Fifth Derivative

Computes the fifth derivative of the filtered predictor in the caller's
unknowns, contracted against the caller's weights and three directions,
together with the lower orders' contractions along the same directions.

## Arguments

- score_values:

  `NULL`, or the driving score at every observation, read instead of
  evaluating `score`.

## Value

As
[`term_fifth()`](https://statmodels7.github.io/modelterms7/reference/term_fifth.md).

## Details

The jets of the level, of the score and of the chart are carried along
every subset of the three directions, and each product of the recursion
is differentiated by the product rule over those subsets. The chart's
derivatives come from
[`.gas_chart_contract()`](https://statmodels7.github.io/modelterms7/reference/dot-gas_chart_contract.md)
and
[`gas_levinson5()`](https://statmodels7.github.io/modelterms7/reference/gas_levinson5.md),
the score's from the family's arrays by
[`.gas_score_jet()`](https://statmodels7.github.io/modelterms7/reference/dot-gas_compose.md).

The compiled kernel
([`.gas_fifth_compiled()`](https://statmodels7.github.io/modelterms7/reference/dot-gas_fifth_compiled.md))
is the route; the R functions
[`.gas_fifth_scalar()`](https://statmodels7.github.io/modelterms7/reference/dot-gas_fifth_run.md)
and
[`.gas_fifth_sub()`](https://statmodels7.github.io/modelterms7/reference/dot-gas_fifth_run.md)
are an independent computation kept for the tests.
