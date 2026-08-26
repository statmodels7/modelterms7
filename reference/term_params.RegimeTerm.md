# The Parameters of a Regime Term

`"level1"`, then `"gap2"` ... `"gapK"`, then the free names of the
transition matrix, `"alr1.1"` and the rest. A chain of \\K\\ regimes has
\\K + K(K-1)\\ parameters in all: one free level, \\K-1\\ positive gaps
that order the others, and \\K(K-1)\\ additive log-ratios, \\K-1\\ per
row of the matrix.

## Arguments

- term:

  A
  [`RegimeTerm()`](https://statmodels7.github.io/modelterms7/reference/RegimeTerm.md),
  built or not: the names depend on `k` alone.

- ...:

  Unused.

## Value

A character vector of length \\K + K(K-1)\\.

## Details

The order matters: it is the order
[`term_links()`](https://statmodels7.github.io/modelterms7/reference/term_links.md),
[`term_start()`](https://statmodels7.github.io/modelterms7/reference/term_start.md),
[`term_loglik()`](https://statmodels7.github.io/modelterms7/reference/term_loglik.md)'s
Jacobian columns and the joint variance matrix are all indexed by, and
`psi` must supply the names whatever order they come in.

The log-ratio names come from
[`parameters7::transition_matrix()`](https://statmodels7.github.io/parameters7/reference/transition_matrix.html)'s
own `free_names`, so the spelling is that package's.

## See also

[`term_links()`](https://statmodels7.github.io/modelterms7/reference/term_links.md)
for the chart each rides,
[`regime()`](https://statmodels7.github.io/modelterms7/reference/regime.md)
for what they mean.

## Examples

``` r
term_params(regime(2))
#> [1] "level1" "gap2"   "alr1.1" "alr2.1"
term_params(regime(3))
#> [1] "level1" "gap2"   "gap3"   "alr1.1" "alr1.2" "alr2.1" "alr2.2" "alr3.1"
#> [9] "alr3.2"
vapply(2:4, function(k) length(term_params(regime(k))), integer(1))
#> [1]  4  9 16
```
