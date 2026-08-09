# Whether a Term's Penalized Objective Is Smooth

`TRUE` when the term's contribution to the penalized objective is
differentiable in the coefficients. The answer is read from the penalty
rather than declared by the term: an unpenalized term is smooth, and a
penalized one is smooth exactly when its penalty declares no kinks, so a
term cannot disagree with its own penalty. The model layer uses this
flag to split the coefficient vector into the block the classical
optimizers handle and the block that needs non-smooth strategies.

## Usage

``` r
term_smooth(term, ...)
```

## Arguments

- term:

  An object inheriting from class
  [`additive_term`](https://statmodels7.github.io/modelterms7/reference/additive_term.md).

- ...:

  Passed to methods.

## Value

A logical scalar.

## See also

[`term_penalty`](https://statmodels7.github.io/modelterms7/reference/term_penalty.md),
[`edf`](https://statmodels7.github.io/modelterms7/reference/edf.md)

## Examples

``` r
term_smooth(linpar(~x))
#> [1] TRUE
```
