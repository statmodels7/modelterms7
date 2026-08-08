# Number of Coefficients of a Built Term

The number of columns of the term's design block.

## Usage

``` r
term_npar(term, ...)
```

## Arguments

- term:

  A built term (see
  [`term_build`](https://statmodels7.github.io/modelterms7/reference/term_build.md)).

- ...:

  Passed to methods.

## Value

An integer.

## Examples

``` r
term_npar(term_build(linpar(~x), data.frame(x = 1:4)))
#> [1] 2
```
