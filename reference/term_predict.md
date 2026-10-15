# Design Block on New Data

Applies a built term's recorded mapping to new rows, returning the block
the term would have produced had those rows been in the data on which it
was built. Factor levels, contrasts, spline knots and a Demmler-Reinsch
reparametrization come from the blueprint and are reused; the spreads of
a standardization enter the penalty only and are not applied to the
block. A factor level that the blueprint does not contain is rejected,
except for the grouping levels of a
[`random()`](https://statmodels7.github.io/modelterms7/reference/random.md)
term with `unseen = "zero"`.

## Usage

``` r
term_predict(term, newdata, ...)
```

## Arguments

- term:

  A built additive term (see
  [`term_build()`](https://statmodels7.github.io/modelterms7/reference/term_build.md)).
  A specification signals the error
  `"the term has not been built; call term_build(term, data) first."`.

- newdata:

  A data frame carrying every variable the term names, with any number
  of rows. Any other object is rejected by the generic, before dispatch,
  with the error `"'newdata' must be a data frame."`. A factor here need
  carry only the levels that its own rows use; the rest come from the
  blueprint.

- ...:

  Passed to methods. The method for
  [`random()`](https://statmodels7.github.io/modelterms7/reference/random.md)
  reads `unseen`, which is `"error"` by default and rejects a grouping
  level that was absent at the build and a missing grouping value; with
  `unseen = "zero"` such a row receives a row of zeros. The other
  shipped methods do not use it.

## Value

A block of `nrow(newdata)` rows and
[`term_npar()`](https://statmodels7.github.io/modelterms7/reference/term_npar.md)
columns, in the same storage the term built: a numeric matrix, or a
Matrix object where the block is sparse.

## The identity that makes it useful

The block returned is \\\tilde{X}\_t\\ such that \\\tilde{\eta} =
\tilde{X}\_t \beta_t\\ is the term's contribution at the new rows, at
the coefficients the model already carries. The identity holds because
the mapping is reused. A rebuilt encoding gives a block of the same
shape multiplying the same coefficients and meaning something else: a
factor whose new rows omit a level loses a column, and a basis rebuilt
on a narrower range is a different set of functions.

For a smooth, predicting on a subset of the fitting rows reproduces the
corresponding rows of the original block exactly, whereas rebuilding the
term on the subset gives a block that multiplies the coefficients by
different functions. The subset check of
[`check_term()`](https://statmodels7.github.io/modelterms7/reference/check_term.md)
makes this comparison.

Predicting on the fitting data returns the block itself, so
`term_predict(b, data)` and `term_matrix(b)` agree to the last bit.

## Which terms have a method

The methods are registered for
[`linpar()`](https://statmodels7.github.io/modelterms7/reference/linpar.md),
the penalized terms,
[`random()`](https://statmodels7.github.io/modelterms7/reference/random.md),
the smooths,
[`nl()`](https://statmodels7.github.io/modelterms7/reference/nl.md) and
the break-point terms. A structural term contributes no block and has no
method, so `term_predict()` on one signals S7's method-not-found error;
[`term_continue()`](https://statmodels7.github.io/modelterms7/reference/term_continue.md)
is the corresponding operation there.

## See also

[`term_matrix()`](https://statmodels7.github.io/modelterms7/reference/term_matrix.md)
for the block on the fitting data,
[`term_build()`](https://statmodels7.github.io/modelterms7/reference/term_build.md)
for what records the blueprint,
[`check_term()`](https://statmodels7.github.io/modelterms7/reference/check_term.md)
for the check this identity is the subject of, and
[`term_continue()`](https://statmodels7.github.io/modelterms7/reference/term_continue.md)
for a structural term.

## Examples

``` r
d <- data.frame(x = 1:6, g = factor(rep(c("a", "b", "c"), 2)))
b <- term_build(linpar(~ x + g), d)

# New rows, the same mapping.
term_predict(b, data.frame(x = c(0.5, 2.5), g = factor(c("a", "c"))))
#>   (Intercept)   x gb gc
#> 1           1 0.5  0  0
#> 2           1 2.5  0  1

# On the fitting data it returns the block itself.
all.equal(term_predict(b, d), term_matrix(b))
#> [1] TRUE

# A subset that drops a level keeps the blueprint's columns, where a
# rebuild would lose one.
nd <- droplevels(d[d$g != "c", ])
dim(term_predict(b, nd))
#> [1] 4 4
dim(model.matrix(~ x + g, nd))
#> [1] 4 3

# A basis is not replaced on the narrower range: reapplying agrees with
# the original rows exactly, rebuilding does not.
d2  <- data.frame(x = seq(0, 1, length.out = 60))
bs  <- term_build(s(x, basis7::bspline_smooth(k = 6)), d2)
X   <- term_matrix(bs)
sub <- 1:20
max(abs(term_predict(bs, d2[sub, , drop = FALSE]) - X[sub, ]))
#> [1] 0
max(abs(term_matrix(term_build(s(x, basis7::bspline_smooth(k = 6)), d2[sub, , drop = FALSE])) -
        X[sub, ]))
#> [1] 3.487998
```
