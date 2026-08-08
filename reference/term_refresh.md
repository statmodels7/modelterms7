# Refresh a Term at New Coefficients

Recomputes whatever a term's design block depends on when the
coefficients move. For every ordinary term this is the identity, the
block being a function of the data alone; for a nonlinear term the block
is the Jacobian of its contribution, which is a function of where the
parameters currently are.

## Usage

``` r
term_refresh(term, coef, ...)
```

## Arguments

- term:

  A built term.

- coef:

  The current coefficients of the term's block.

- ...:

  Passed to methods.

## Value

A built term, refreshed.

## See also

[`nl`](https://statmodels7.github.io/modelterms7/reference/nl.md),
[`term_value`](https://statmodels7.github.io/modelterms7/reference/term_value.md)

## Examples

``` r
dd <- data.frame(x = seq(0, 2, length.out = 20))
built <- term_build(nl(~ a * exp(-r * x), start = list(a = 1, r = 1)), dd)
r1 <- term_refresh(built, c(2, 0.5))
max(abs(term_matrix(r1) - term_matrix(built))) > 0
#> [1] TRUE
```
