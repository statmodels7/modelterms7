# The Posterior of a Structural Term's Latent Variable

A summary of the latent variable a structural term integrates over,
given the whole sample. For a marginal break-point term it is the
posterior mean and standard deviation of each group's break-points. That
is what a reader wants from such a fit: where each group's change
happened, and how sure the data are about it.

## Usage

``` r
term_latent(term, eta, y, logdens, psi, ...)
```

## Arguments

- term:

  A built structural term.

- eta:

  The static predictor of the equation the term sits in, one value per
  observation.

- y:

  The response.

- logdens:

  A function `(e, i)` returning the log-density of observation `i` at
  predictor `e`, as
  [`term_loglik()`](https://statmodels7.github.io/modelterms7/reference/term_loglik.md)
  takes it.

- psi:

  The term's parameters on the parameter scale, named as
  [`term_params()`](https://statmodels7.github.io/modelterms7/reference/term_params.md).

- ...:

  Passed to methods.

## Value

A data frame with one row per group and break-point and four columns:
`group`, the grouping level; `psi`, which break-point; `mean` and `sd`,
the posterior moments of its position. `NA` in a moment the prior does
not possess.

## Details

[`term_posterior()`](https://statmodels7.github.io/modelterms7/reference/term_posterior.md)
answers the fitting layer's question, the component weights Fisher's
identity needs at every observation. This one answers the reader's. For
the marginal break-point term the two come from the same decomposition:
the mean and variance within an interval are those of the prior
truncated to it, and under quadrature the moments of the node posterior.

## What a heavy-tailed prior can refuse

The moments are the prior's, so a prior without them has none to report.
A Student t below one degree of freedom has no mean on an edge interval
and below two no variance, and the quadrature returns `NA` there instead
of a number. That is a property of the prior the caller chose.

The method on
[`structural_term()`](https://statmodels7.github.io/modelterms7/reference/structural_term.md)
throws, naming the class:
[`gas()`](https://statmodels7.github.io/modelterms7/reference/gas.md)
and
[`regime()`](https://statmodels7.github.io/modelterms7/reference/regime.md)
have no continuous latent to summarize this way, a regime's latent being
the discrete state
[`term_posterior()`](https://statmodels7.github.io/modelterms7/reference/term_posterior.md)
already reports.

## See also

[`term_posterior()`](https://statmodels7.github.io/modelterms7/reference/term_posterior.md)
for the component weights a fit reads,
[`jump()`](https://statmodels7.github.io/modelterms7/reference/jump.md)
and
[`seg()`](https://statmodels7.github.io/modelterms7/reference/seg.md)
for the terms that implement it.

## Examples

``` r
set.seed(1)
dd <- data.frame(id = rep(1:3, each = 8), x = rep(1:8, 3))
dd$y <- rnorm(24, 2 * (dd$x >= 4.5), 0.4)
tm <- term_build(jump(x, psi ~ random(~ 1 | id), marginal = TRUE), dd)

# Every group's step is at 4.5, and the posterior finds it there,
# with a spread well inside the prior's own 0.5.
term_latent(tm, rep(0, 24), dd$y,
            logdens = function(e, i) dnorm(dd$y[i], e, 0.4, log = TRUE),
            psi = list(m1 = 4.5, tau1 = 0.5, delta1 = 2))
#>   group psi     mean        sd
#> 1     1   1 4.498117 0.2722611
#> 2     2   1 4.500010 0.2698055
#> 3     3   1 4.499987 0.2697968
```
