# The Fifth Derivative of a Joint Score-Driven Filter

Computes the fifth derivative of the filtered predictors contracted
against the caller's weights and three directions, with the lower
orders' contractions along the same directions, through jets over every
subset of the directions.

## Usage

``` r
.gas_joint_fifth(
  term,
  psi,
  g,
  seed,
  blocks_data,
  directions,
  eta = NULL,
  sub = FALSE
)
```

## Arguments

- term:

  A built term.

- psi:

  The term's parameters, named as
  [`term_params()`](https://statmodels7.github.io/modelterms7/reference/term_params.md).

- g:

  The weights, one column per filtered parameter.

- seed:

  A list of one matrix of static rows per equation.

- blocks_data:

  The model's derivative arrays, as set out above.

- directions:

  A list of three numeric vectors, the directions \\v\\, \\w\\ and \\x\\
  in the caller's unknowns. The result is symmetric in the three.

- eta:

  The static part of the predictor.

- sub:

  Whether the members carry subformulas, whose charts are then read per
  observation.

## Value

A list with `jacobian`, `dphi`, `dpsi` and `domega`, each with one
element per filtered parameter (`dphi` and `dpsi` lists of three), and
`curvature`.

## Details

The levels' jets are sums of products of the chart's jets with the jets
of earlier driving scores and levels, as at the lower orders of
[`term_fourth()`](https://statmodels7.github.io/modelterms7/reference/term_fourth.md).
The jet of each driving score is composed from the jets of the
predictors of every equation by
[`.gas_compose()`](https://statmodels7.github.io/modelterms7/reference/dot-gas_compose.md),
with the driving score's derivatives in those predictors given as data
in `blocks_data`: `filtered`, `u`, `U1`, `U2`, `U3`, `U4` as for the
lower orders, and `U5`, an array with the observation, the driving score
and five predictor indices.
