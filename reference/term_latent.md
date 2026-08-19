# The Posterior of a Structural Term's Latent Variable

A summary of the latent variable a structural term integrates over,
given the whole sample: for a marginal break-point term, the posterior
mean and standard deviation of each group's break-points.

## Usage

``` r
term_latent(term, eta, y, logdens, psi, ...)
```

## Arguments

- term:

  A built structural term.

- eta:

  The static predictor of the equation the term sits in.

- y:

  The response.

- logdens:

  The log-density as a function of a predictor value and a row index, as
  [`term_loglik`](https://statmodels7.github.io/modelterms7/reference/term_loglik.md)
  takes it.

- psi:

  The term's parameters, named as
  [`term_params`](https://statmodels7.github.io/modelterms7/reference/term_params.md).

- ...:

  Passed to methods.

## Value

A data frame with one row per group and break-point: `group`, `psi`
(which break-point), `mean` and `sd`.

## Details

[`term_posterior`](https://statmodels7.github.io/modelterms7/reference/term_posterior.md)
answers the fitting layer's question, the component weights Fisher's
identity needs at every observation. This one answers the reader's:
where each group's latent positions sit once the data have been seen.
For the marginal break-point term the two come from the same
decomposition; the mean and variance within an interval are those of the
prior truncated to it, and under quadrature the moments of the node
posterior.

## See also

[`term_posterior`](https://statmodels7.github.io/modelterms7/reference/term_posterior.md),
[`jump`](https://statmodels7.github.io/modelterms7/reference/jump.md)

## Examples

``` r
set.seed(1)
dd <- data.frame(id = rep(1:3, each = 8), x = rep(1:8, 3))
dd$y <- rnorm(24, 2 * (dd$x >= 4.5), 0.4)
tm <- term_build(jump(x, psi ~ random(~1 | id), marginal = TRUE), dd)
term_latent(tm, rep(0, 24), dd$y,
            logdens = function(e, i) dnorm(dd$y[i], e, 0.4, log = TRUE),
            psi = list(m1 = 4.5, tau1 = 0.5, delta1 = 2))
#>   group psi     mean        sd
#> 1     1   1 4.498117 0.2722611
#> 2     2   1 4.500010 0.2698055
#> 3     3   1 4.499987 0.2697968
```
