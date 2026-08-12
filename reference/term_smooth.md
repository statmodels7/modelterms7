# Whether a Term's Penalized Objective Is Smooth

`TRUE` when the term's contribution to the penalized objective is
differentiable in the coefficients. The answer is read from the
penalties rather than declared by the term: an unpenalized term is
smooth, and a penalized one is smooth exactly when no penalty it carries
declares a kink, so a term cannot disagree with its own penalties. The
model layer uses this flag to split the coefficient vector into the
block the classical optimizers handle and the block that needs
non-smooth strategies.

## Usage

``` r
term_smooth(term, ...)
```

## Arguments

- term:

  An object inheriting from class
  [`model_term`](https://statmodels7.github.io/modelterms7/reference/model_term.md).

- ...:

  Passed to methods.

## Value

A logical scalar.

## Details

The enumeration is
[`term_penalties`](https://statmodels7.github.io/modelterms7/reference/term_penalties.md),
so a term carrying one penalty over part of its parameters and none over
the rest answers for the part: `seg(x, penalty = "lasso")` is not
smooth, its slope changes sitting at a kink, although its linear effect
and its break-points are unpenalized.

## See also

[`term_penalties`](https://statmodels7.github.io/modelterms7/reference/term_penalties.md),
[`term_penalty`](https://statmodels7.github.io/modelterms7/reference/term_penalty.md),
[`edf`](https://statmodels7.github.io/modelterms7/reference/edf.md)

## Examples

``` r
term_smooth(linpar(~x))
#> [1] TRUE
```
