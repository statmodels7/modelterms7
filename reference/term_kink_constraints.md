# The Linear Constraints a Kink Places on a Term's Coefficients

The rows of the linear constraints that hold the term's objective on its
kinks at the coefficients given: one row per constraint, one column per
coefficient of the term. A perturbation \\d\\ of the coefficients keeps
the kinks where they are if and only if \\A d = 0\\.

## Usage

``` r
term_kink_constraints(term, coef = NULL, ...)
```

## Arguments

- term:

  A built term.

- coef:

  The term's coefficients, `NULL` for the stored ones.

- ...:

  Passed to methods. The shipped methods do not use it.

## Value

A numeric matrix with one column per coefficient of the term and one row
per constraint, with no rows where the term has no kink at `coef`.

## Details

Where
[`term_kinks()`](https://statmodels7.github.io/modelterms7/reference/term_kinks.md)
names a coordinate, the kink is held by holding that coordinate, and the
row is the corresponding unit vector; that is the base method's answer.
A break-point whose position is developed over a design, \\\psi_i =
z_i^\top\kappa\\, sits on an observation \\j\\ along the hyperplane
\\z_j^\top\kappa = x_j\\, so the row of that constraint is \\z_j\\ and
the coefficients can still move along the hyperplane. The
[`seg()`](https://statmodels7.github.io/modelterms7/reference/seg.md)
method returns those rows and the
[`nl()`](https://statmodels7.github.io/modelterms7/reference/nl.md)
method carries the rows of a
[`seg()`](https://statmodels7.github.io/modelterms7/reference/seg.md)
sub-term to the term's coefficients. Rows that repeat are given once.

## See also

[`term_kinks()`](https://statmodels7.github.io/modelterms7/reference/term_kinks.md),
[`seg()`](https://statmodels7.github.io/modelterms7/reference/seg.md).

## Examples

``` r
d <- data.frame(x = seq(0, 10, length.out = 41))
term_kink_constraints(term_build(seg(x), d), c(0.5, 2, 6))
#>      [,1] [,2] [,3]
#> [1,]    0    0    1
term_kink_constraints(term_build(seg(x), d), c(0.5, 2, 6.1))
#>      [,1] [,2] [,3]
```
