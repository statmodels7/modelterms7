# Where a Term's Own Parameters Start

The starting values of a structural term's parameters, on the
unconstrained scale of
[`term_links`](https://statmodels7.github.io/modelterms7/reference/term_links.md),
one per parameter of
[`term_params`](https://statmodels7.github.io/modelterms7/reference/term_params.md).

## Usage

``` r
term_start(term, ...)
```

## Arguments

- term:

  A built structural term.

- ...:

  Passed to methods.

## Value

A named numeric vector on the unconstrained scale.

## Details

The start belongs to the term because only the term knows what a
coordinate of zero means on each of its charts. The base method returns
zero everywhere, which is each link's own natural point; a term whose
chart makes zero mean something other than "the model without the term"
overrides it, as
[`gas`](https://statmodels7.github.io/modelterms7/reference/gas.md) does
for its score loadings, whose log chart puts a loading of one at zero.

## See also

[`term_params`](https://statmodels7.github.io/modelterms7/reference/term_params.md),
[`term_links`](https://statmodels7.github.io/modelterms7/reference/term_links.md)

## Examples

``` r
term_start(gas(p = 1, q = 1))
#>     omega    alpha1     pacf1 
#>  0.000000 -2.302585  0.000000 
```
