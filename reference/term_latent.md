# The Posterior of a Structural Term's Latent Variable

Returns a summary of the latent variable over which a structural term
integrates, given the whole sample. For a marginal break-point term it
consists of the posterior mean and standard deviation of each group's
break-points, which locate each group's change and measure the
uncertainty about its position. For a
[`regime()`](https://statmodels7.github.io/modelterms7/reference/regime.md)
term it is the smoothed probability of each regime at each observation.

## Usage

``` r
term_latent(term, eta, y, logdens, psi, ...)
```

## Arguments

- term:

  A built structural term.

- eta:

  The static predictor of the equation in which the term sits, one value
  per observation.

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

For a marginal break-point term, a data frame with one row per group and
break-point and four columns: `group`, the grouping level; `psi`, which
break-point; `mean` and `sd`, the posterior moments of its position. A
moment that the prior does not possess is `NA`. For a regime term, a
data frame with one row per observation and one column per regime,
`state1`, `state2`, and so on.

## Details

[`term_posterior()`](https://statmodels7.github.io/modelterms7/reference/term_posterior.md)
returns the component weights that Fisher's identity needs at every
observation, for use by a fitting layer. This function returns a summary
of the latent variable for the reader of a fit. For the marginal
break-point term the two come from the same decomposition: the mean and
variance within an interval are those of the prior truncated to it, and
under quadrature the moments of the node posterior.

## Priors without moments

The moments are those of the prior, so a prior that lacks them cannot
give them. A Student t prior with fewer than one degree of freedom has
no mean on an edge interval, and with fewer than two it has no variance.
Where such an interval carries posterior weight of at least 1e-9, the
corresponding column holds `NA`.

The method on
[`structural_term()`](https://statmodels7.github.io/modelterms7/reference/structural_term.md)
signals an error that names the class. A
[`gas()`](https://statmodels7.github.io/modelterms7/reference/gas.md)
term has no latent variable, because its level is a deterministic
function of the data. The latent variable of a regime term is the
discrete state, and its method returns the probabilities that
[`term_posterior()`](https://statmodels7.github.io/modelterms7/reference/term_posterior.md)
computes, labeled by regime.

## See also

[`term_posterior()`](https://statmodels7.github.io/modelterms7/reference/term_posterior.md)
for the component weights that a fit reads,
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
            psi = list(psi1.mean = 4.5, psi1.sd = 0.5, delta1 = 2))
#>   group psi     mean        sd
#> 1     1   1 4.498117 0.2722611
#> 2     2   1 4.500010 0.2698055
#> 3     3   1 4.499987 0.2697968
```
