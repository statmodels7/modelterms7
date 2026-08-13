# Parameters of a Structural Term

The names of a structural term's own parameters, in the order its filter
expects them. A structural term contributes no design block, so its
parameters are not coefficients: they are estimated alongside the
distribution's, on the unconstrained scale its links define.

## Usage

``` r
term_params(term, ...)
```

## Arguments

- term:

  An object inheriting from
  [`structural_term`](https://statmodels7.github.io/modelterms7/reference/structural_term.md).

- ...:

  Passed to methods.

## Value

A character vector.

## See also

[`term_links`](https://statmodels7.github.io/modelterms7/reference/term_links.md),
[`term_filter`](https://statmodels7.github.io/modelterms7/reference/term_filter.md)

## Examples

``` r
term_params(gas(p = 1, q = 1))
#> [1] "omega"  "alpha1" "pacf1" 
```
