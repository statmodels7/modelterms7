# Links of a Structural Term's Parameters

One linkfunctions7 link per parameter of
[`term_params`](https://statmodels7.github.io/modelterms7/reference/term_params.md),
carrying it to the unconstrained scale the model layer optimizes on.

## Usage

``` r
term_links(term, ...)
```

## Arguments

- term:

  An object inheriting from
  [`structural_term`](https://statmodels7.github.io/modelterms7/reference/structural_term.md).

- ...:

  Passed to methods.

## Value

A named list of link objects.

## See also

[`term_params`](https://statmodels7.github.io/modelterms7/reference/term_params.md)

## Examples

``` r
vapply(term_links(gas(p = 1, q = 1)), function(l) l@link_name, character(1))
#>      omega     alpha1      pacf1 
#> "identity" "identity"   "rhobit" 
```
