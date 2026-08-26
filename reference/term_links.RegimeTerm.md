# The Charts of a Regime Term's Parameters

The **log** link on every gap and the identity on everything else. A gap
must be positive, that being what orders the levels; a level and an
additive log-ratio are already unconstrained.

## Arguments

- term:

  A
  [`RegimeTerm()`](https://statmodels7.github.io/modelterms7/reference/RegimeTerm.md).

- ...:

  Unused.

## Value

A named list of linkfunctions7 links, one per entry of
[`term_params()`](https://statmodels7.github.io/modelterms7/reference/term_params.md).

## Details

The log-ratios need no chart because
[`parameters7::transition_matrix()`](https://statmodels7.github.io/parameters7/reference/transition_matrix.html)
has already provided one: they are the additive log-ratios of each row,
so any real values at all give a row of probabilities summing to one.

The consequence for a caller is that `psi` and `zeta` differ only in the
gaps.
[`term_loglik()`](https://statmodels7.github.io/modelterms7/reference/term_loglik.md)
takes the parameter scale, where `gap2 = 3` is a gap of three;
[`term_readable()`](https://statmodels7.github.io/modelterms7/reference/term_readable.md)
takes the unconstrained scale, where the same number is a gap of
\\e^3\\.

## See also

[`term_params()`](https://statmodels7.github.io/modelterms7/reference/term_params.md),
[`regime()`](https://statmodels7.github.io/modelterms7/reference/regime.md).

## Examples

``` r
vapply(term_links(regime(3)), function(l) l@link_name, character(1))
#>     level1       gap2       gap3     alr1.1     alr1.2     alr2.1     alr2.2 
#> "identity"      "log"      "log" "identity" "identity" "identity" "identity" 
#>     alr3.1     alr3.2 
#> "identity" "identity" 
```
