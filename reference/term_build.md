# Build a Term on Data

Turns a term specification into a built term: the design block is
computed from the data, the coefficient names are assigned, and the
blueprint that reproduces the mapping on new data is recorded. The
returned object is a copy of the specification with those properties
filled; the specification itself is unchanged.

## Usage

``` r
term_build(term, data, ...)
```

## Arguments

- term:

  An object inheriting from class
  [`model_term`](https://statmodels7.github.io/modelterms7/reference/model_term.md).

- data:

  A data frame.

- ...:

  Passed to methods.

## Value

A built term of the same class as `term`.

## Examples

``` r
built <- term_build(linpar(~x), data.frame(x = 1:4))
term_matrix(built)
#>   (Intercept) x
#> 1           1 1
#> 2           1 2
#> 3           1 3
#> 4           1 4
```
