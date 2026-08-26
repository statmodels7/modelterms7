# The Levels of a Likelihood-Shaped Structural Term

The shifts a term of the likelihood shape adds to its equation's
predictor, one per mixture component, in the order the columns of
[`term_posterior()`](https://statmodels7.github.io/modelterms7/reference/term_posterior.md)
carry the components. A fitting layer reads the two together to assemble
Fisher's identity.

## Usage

``` r
term_levels(term, psi, ...)
```

## Arguments

- term:

  A built structural term of the likelihood shape.

- psi:

  The term's parameters on the parameter scale, named as
  [`term_params()`](https://statmodels7.github.io/modelterms7/reference/term_params.md).

- ...:

  Passed to methods.

## Value

A numeric vector of one level per component, or a numeric matrix of `n`
rows and one column per component where a shift varies by observation.
The order matches
[`term_posterior()`](https://statmodels7.github.io/modelterms7/reference/term_posterior.md)'s
columns.

## Details

By Fisher's identity the derivative of a likelihood mixed over latent
states, in **any** predictor the model carries, is the
posterior-weighted derivative of the ordinary one, each component read
at the predictor shifted by its own level:

\$\$\frac{\partial L}{\partial \eta\_{q,t}} = \sum_k \gamma_t(k)\\
\frac{\partial \ell(y_t; \theta_t(k))}{\partial \eta_q}.\$\$

[`term_posterior()`](https://statmodels7.github.io/modelterms7/reference/term_posterior.md)
supplies the \\\gamma_t(k)\\ and this supplies the shifts, so a fitting
layer assembles the identity without reading the term's internals.

For
[`regime()`](https://statmodels7.github.io/modelterms7/reference/regime.md)
the levels are the ordered regime means, one number per component. For a
marginal break-point term of the step kind they are the sums of the
changes of level over the active break-points, one number per side
pattern.

## A shift may vary by observation

The quadrature nodes of a marginal
[`seg()`](https://statmodels7.github.io/modelterms7/reference/seg.md) or
[`jseg()`](https://statmodels7.github.io/modelterms7/reference/jseg.md)
term shift each observation by its own hinge value, so the method may
return a **matrix** of one row per observation and one column per
component. A caller reads a column of it wherever it would read a
constant level, and must accept both shapes.

The method on
[`structural_term()`](https://statmodels7.github.io/modelterms7/reference/structural_term.md)
throws, naming the class: a term of the filter shape reports a predictor
instead of components, and answers
[`term_filter()`](https://statmodels7.github.io/modelterms7/reference/term_filter.md).

## See also

[`term_posterior()`](https://statmodels7.github.io/modelterms7/reference/term_posterior.md)
for the weights,
[`term_loglik()`](https://statmodels7.github.io/modelterms7/reference/term_loglik.md)
for the likelihood,
[`regime()`](https://statmodels7.github.io/modelterms7/reference/regime.md)
and
[`jump()`](https://statmodels7.github.io/modelterms7/reference/jump.md)
for the two implementations.

## Examples

``` r
# A two-regime chain: the first level and the cumulated gap.
term_levels(regime(2), list(level1 = 0, gap2 = 3, alr1.1 = 2, alr2.1 = -2))
#> [1] 0 3

# A one-break-point step term: no shift where no break-point is active,
# and the change of level where one is.
set.seed(1)
dd <- data.frame(id = rep(1:3, each = 8), x = rep(1:8, 3))
dd$y <- rnorm(24, 2 * (dd$x >= 4.5), 0.4)
tm <- term_build(jump(x, psi ~ random(~ 1 | id), marginal = TRUE), dd)
term_levels(tm, list(m1 = 4.5, tau1 = 0.5, delta1 = 2))
#> [1] 0 2
```
