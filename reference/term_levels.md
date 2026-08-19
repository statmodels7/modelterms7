# The Levels of a Likelihood-Shaped Structural Term

The shifts a term of the likelihood shape adds to its equation's
predictor, one per mixture component, in the order the columns of
[`term_posterior`](https://statmodels7.github.io/modelterms7/reference/term_posterior.md)
carry the components.

## Usage

``` r
term_levels(term, psi, ...)
```

## Arguments

- term:

  A built structural term of the likelihood shape.

- psi:

  The term's parameters, named as
  [`term_params`](https://statmodels7.github.io/modelterms7/reference/term_params.md).

- ...:

  Passed to methods.

## Value

A numeric vector with one level per component, or a matrix with one row
per observation and one column per component.

## Details

By Fisher's identity the derivative of a likelihood mixed over latent
states, in any predictor the model carries, is the posterior-weighted
derivative of the ordinary one, each component read at the predictor
shifted by its own level.
[`term_posterior`](https://statmodels7.github.io/modelterms7/reference/term_posterior.md)
supplies the weights; this supplies the levels, so a fitting layer
assembles the identity without reading the term's internals. For
[`regime`](https://statmodels7.github.io/modelterms7/reference/regime.md)
the levels are the ordered regime means, one number per component; for a
marginal break-point term of the step kind they are the sums of the
changes of level over the active break-points, one number per side
pattern.

A component's shift may vary by observation – the quadrature nodes of a
marginal
[`seg`](https://statmodels7.github.io/modelterms7/reference/seg.md) or
[`jseg`](https://statmodels7.github.io/modelterms7/reference/jseg.md)
term shift each observation by its own hinge value – and the method then
returns a matrix with one row per observation and one column per
component, whose columns a caller reads in place of the constant levels.

## See also

[`term_posterior`](https://statmodels7.github.io/modelterms7/reference/term_posterior.md),
[`term_loglik`](https://statmodels7.github.io/modelterms7/reference/term_loglik.md)

## Examples

``` r
term_levels(regime(2), list(level1 = 0, gap2 = 3,
                            alr1.1 = 2, alr2.1 = -2))
#> [1] 0 3
```
