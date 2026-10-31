# The Jet of a Function of Several Predictors

Composes the derivatives of a function \\u(e_1, \dots, e_P)\\ of the
predictors of every equation with the jets of those predictors, by Faa
di Bruno's formula over the set partitions of the slots: the directions
of a subset and the two free slots. `.gas_score_jet()` is the case of
one filtered equation, whose predictor alone has derivatives above the
first.

## Usage

``` r
.gas_compose(ell, s0, E)

.gas_score_jet(ell, s0, ej, V, ap)
```

## Arguments

- ell:

  A list whose element \\k\\ is the array of the derivatives of order
  \\k\\, with \\k\\ dimensions of `P` each.

- s0:

  The value.

- E:

  A list of `P` jets, one per equation's predictor.

- ej:

  The jet of the filtered predictor.

- V:

  A `P` by `m` matrix of the predictors' first derivatives (the row `ap`
  is replaced by the jet's).

- ap:

  The index of the filtered equation.

## Value

A jet.

## Details

A partition with \\k\\ blocks contributes \\\sum\_{r_1,\dots,r_k}
u\_{r_1\dots r_k}\prod\_\beta
\partial^{\|\beta\|}e\_{r\_\beta}\[\beta\]\\, where each factor is the
derivative of one predictor along the slots of its block: a number for a
block without a free slot, a row for a block with one, a matrix for a
block with both. The coefficients \\u\_{r_1\dots r_k}\\ are read from
`ell`, the arrays of the driving quantity's derivatives in the
predictors at the observation (`H`, `D3`, ..., `D6` of
[`term_fifth()`](https://statmodels7.github.io/modelterms7/reference/term_fifth.md),
or `U1` to `U5` of a joint filter).
