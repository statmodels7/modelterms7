# The Charts of a Marginal Break-Point Term's Parameters

The **log** link on every `tauk`, a prior's own link on any parameter it
contributes, and the identity on everything else. A prior scale must be
positive; a position, a change of level and a change of slope are
already unconstrained.

## Arguments

- term:

  A
  [`MarginalBreakTerm()`](https://statmodels7.github.io/modelterms7/reference/MarginalBreakTerm.md),
  built or not.

- ...:

  Unused.

## Value

A named list of linkfunctions7 links, one per entry of
[`term_params()`](https://statmodels7.github.io/modelterms7/reference/term_params.md).

## Details

Where `random(distrib = )` named a family other than the Gaussian, that
family's parameters keep the links it declares in `link_params`, so a
Student t prior's degrees of freedom ride whatever chart distributions7
gives them. Nothing about the chart is restated here.

## See also

[`term_params()`](https://statmodels7.github.io/modelterms7/reference/term_params.md),
[`term_start()`](https://statmodels7.github.io/modelterms7/reference/term_start.md).

## Examples

``` r
set.seed(1)
dd <- data.frame(id = rep(1:3, each = 8), x = rep(1:8, 3))
dd$y <- rnorm(24, 2 * (dd$x >= 4.5), 0.4)
tm <- term_build(jump(x, psi ~ random(~ 1 | id), marginal = TRUE), dd)
vapply(term_links(tm), function(l) l@link_name, character(1))
#>         m1       tau1     delta1 
#> "identity"      "log" "identity" 
```
