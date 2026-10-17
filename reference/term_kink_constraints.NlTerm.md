# The Constraints of a Kinked Sub-Term, Carried to a Nonlinear Term

Returns the rows of
[`term_kink_constraints()`](https://statmodels7.github.io/modelterms7/reference/term_kink_constraints.md)
of every sharp
[`seg()`](https://statmodels7.github.io/modelterms7/reference/seg.md)
sub-term that develops a parameter, placed in the columns of the term's
coefficients that the sub-term owns. Every other nonlinear term has no
rows.

## Arguments

- term:

  A built
  [`NlTerm()`](https://statmodels7.github.io/modelterms7/reference/NlTerm.md).

- coef:

  The term's coefficients. `NULL` reads the stored ones.

- ...:

  Unused.

## Value

A numeric matrix with one column per coefficient of the term and one row
per constraint.

## Examples

``` r
d <- data.frame(t = seq(0, 10, length.out = 41), x = rep(1:2, length.out = 41))
b <- term_build(nl(~ a * exp(-r * x), r ~ seg(t),
                   start = list(a = 2, r = 0.5)), d)
term_kink_constraints(b)
#>      [,1] [,2] [,3] [,4] [,5]
#> [1,]    0    0    0    0    1
```
