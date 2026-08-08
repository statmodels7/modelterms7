# Penalty of a Term

The penalty attached to the term's coefficients, or `NULL` for an
unpenalized term. The hyperparameters, their bounds and links, and every
derivative in the coefficients and the hyperparameters are the penalty
object's, not the term's.

## Usage

``` r
term_penalty(term, ...)
```

## Arguments

- term:

  An object inheriting from class
  [`additive_term`](https://statmodels7.github.io/modelterms7/reference/additive_term.md).

- ...:

  Passed to methods.

## Value

A penalty object, or `NULL`.

## Examples

``` r
term_penalty(linpar(~x))
#> NULL
```
