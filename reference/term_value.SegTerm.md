# The Contribution of a Break-Point Term

The values the term contributes to the linear predictor at its current
break-points: the broken line for
[`seg()`](https://statmodels7.github.io/modelterms7/reference/seg.md),
the step for
[`jump()`](https://statmodels7.github.io/modelterms7/reference/jump.md),
both for
[`jseg()`](https://statmodels7.github.io/modelterms7/reference/jseg.md).
It is what a Gauss-Newton step needs beside the block, and for the
continuous construction it is **not** the block times the coefficients.

## Arguments

- term:

  A built
  [`SegTerm()`](https://statmodels7.github.io/modelterms7/reference/SegTerm.md).

- coef:

  Optional coefficients to evaluate at; `NULL`, the default, uses the
  ones
  [`term_refresh()`](https://statmodels7.github.io/modelterms7/reference/term_refresh.md)
  last committed.

- newdata:

  An optional data frame; the contribution comes back on its rows
  instead of the fitting ones.

- ...:

  Unused.

## Value

A numeric vector of one value per observation.

## Details

For `seg` the block is the Jacobian, so \\X\beta\\ is the linearization
and the two differ by whatever that drops: a step at the break-point, in
a construction that is continuous. For `jump` the columns satisfy
\\X\beta = \\ the contribution exactly, the identity behind the read-off
making them so, and the two agree to the last bit.

With `newdata` the contribution is returned on those rows, at the
break-points the term **carries**. They are not re-derived from the new
rows, exactly as
[`term_predict()`](https://statmodels7.github.io/modelterms7/reference/term_predict.md)
reapplies rather than rebuilds.

## See also

[`term_value()`](https://statmodels7.github.io/modelterms7/reference/term_value.md)
for the generic,
[`term_refresh.SegTerm()`](https://statmodels7.github.io/modelterms7/reference/term_refresh.SegTerm.md)
for the block at the same coefficients,
[`seg_psi()`](https://statmodels7.github.io/modelterms7/reference/seg_psi.md)
for the positions.

## Examples

``` r
set.seed(1)
d <- data.frame(x = sort(runif(120, 0, 10)))
d$y <- 1 + 0.5 * d$x + 2 * pmax(d$x - 6, 0) + rnorm(120, sd = 0.4)

# seg: the block is the Jacobian, so X beta is the linearization and
# differs from the contribution.
b <- term_refresh(term_build(seg(x, npsi = 1), d), c(0.5, 2, 6))
max(abs(as.numeric(term_matrix(b) %*% b@blueprint$coef) - term_value(b)))
#> [1] 12

# jump: the columns reproduce the contribution exactly.
bj <- term_build(jump(x, npsi = 1), d)
max(abs(as.numeric(term_matrix(bj) %*% bj@blueprint$coef) - term_value(bj)))
#> [1] 4.440892e-16

# And the seg contribution really is the broken line at psi.
cf <- b@blueprint$coef
max(abs(term_value(b) - (cf[1] * d$x + cf[2] * pmax(d$x - cf[3], 0))))
#> [1] 0
```
