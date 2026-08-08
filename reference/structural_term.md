# S7 Class for Structural Terms

The branch of
[`model_term`](https://statmodels7.github.io/modelterms7/reference/model_term.md)
for terms that rewrite the likelihood contributions rather than adding a
design block, such as score-driven dynamics. The class exists so the
formula interpreter can route such terms; no structural term is
implemented yet, and
[`term_build`](https://statmodels7.github.io/modelterms7/reference/term_build.md)
on this branch signals an error saying so.

## Usage

``` r
structural_term(label = character(0))
```

## Arguments

- label:

  A character string naming the term.

## Value

An object inheriting from class `structural_term`.

## See also

[`model_term`](https://statmodels7.github.io/modelterms7/reference/model_term.md)

## Examples

``` r
S7::S7_inherits(linpar(~1), structural_term)
#> [1] FALSE
```
