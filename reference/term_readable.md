# The Quantities a Fitted Term Reports

What a reader reads, with the Jacobian from the term's own parameters on
the unconstrained scale, so that a caller holding their variance matrix
can carry it across by the delta method.

## Usage

``` r
term_readable(term, zeta, ...)
```

## Arguments

- term:

  A built term.

- zeta:

  The term's parameters on the unconstrained scale, named as
  [`term_params`](https://statmodels7.github.io/modelterms7/reference/term_params.md).

- ...:

  Passed to methods.

## Value

A list with `name`, `value`, `jacobian` (one row per quantity and one
column per parameter) and `scale`, the link an interval for each
quantity is built on.

## Details

A term's parameters are the coordinates it is ESTIMATED on, chosen so
that a search runs unconstrained, and they are not always the quantities
the model is about. The clearest case is a score-driven persistence,
which rides a partial autocorrelation because the stationary region is
not a box: what the literature writes as \\\beta_1\\ is the
autoregressive coefficient, which is a function of the whole chart and
coincides with it only at \\q = 1\\. Reporting the coordinate under the
coefficient's name would promise one quantity and print another.

Each row gives a value and the row of \\\partial(\text{value}) /
\partial\zeta\\ at the current parameters, so a standard error is
\\\sqrt{J V J^\top}\\ and an interval is built on whichever scale keeps
the quantity in its own set, exactly as parameters7's `param_readable()`
does for a matrix parameter.

The base method reports the parameters themselves on the PARAMETER
scale, with the diagonal Jacobian of their links, which is what every
term whose coordinates are already its quantities wants.

## See also

[`term_params`](https://statmodels7.github.io/modelterms7/reference/term_params.md),
[`term_links`](https://statmodels7.github.io/modelterms7/reference/term_links.md)

## Examples

``` r
term_readable(gas(p = 1, q = 1), c(omega = 0.3, alpha1 = 0.4, pacf1 = 0.8))
#> $name
#> [1] "omega"  "alpha1" "beta1" 
#> 
#> $value
#> [1] 0.3000000 1.4918247 0.6640368
#> 
#> $jacobian
#>        omega   alpha1     pacf1
#> omega      1 0.000000 0.0000000
#> alpha1     0 1.491825 0.0000000
#> beta1      0 0.000000 0.5590552
#> 
#> $scale
#> $scale$omega
#> S7 Link Object: identity
#>   - Parameter domain (theta): (-Inf, Inf)
#> 
#> $scale$alpha1
#> S7 Link Object: log
#>   - Parameter domain (theta): (0, Inf)
#> 
#> $scale$beta1
#> S7 Link Object: identity
#>   - Parameter domain (theta): (-Inf, Inf)
#> 
#> 
```
