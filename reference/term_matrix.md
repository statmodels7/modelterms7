# Design Block of a Built Term

Returns the \\n \times k\\ block a built additive term contributes to
the linear predictor, with the term's coefficient names as column names
and one row per observation of the data it was built on. The block is
the term's `X` property, returned as it is stored.

## Usage

``` r
term_matrix(term, ...)
```

## Arguments

- term:

  A built additive term (see
  [`term_build()`](https://statmodels7.github.io/modelterms7/reference/term_build.md)).
  A specification throws
  `"the term has not been built; call term_build(term, data) first."`.

- ...:

  Passed to methods. No shipped method reads anything here.

## Value

The design block: a numeric matrix, or a two-dimensional Matrix object
where the term built one, with `nrow` the number of observations the
term was built on and `ncol` equal to
[`term_npar()`](https://statmodels7.github.io/modelterms7/reference/term_npar.md).

## Details

The block is **not necessarily a base matrix**.
[`random()`](https://statmodels7.github.io/modelterms7/reference/random.md)
builds a grouping indicator as a `dgCMatrix`, since a row belongs to one
group and the density is \\1/m\\, and
[`linpar()`](https://statmodels7.github.io/modelterms7/reference/linpar.md),
the penalized terms and a smooth with a factor `by` build sparse when
asked. Code that reads a block therefore tests
`is.matrix(x) && is.numeric(x)` or a two-dimensional S4 object;
[`is.matrix()`](https://rdrr.io/r/base/matrix.html) alone is `FALSE` for
every Matrix class and would reject a block for being economical.

For a term whose block moves with its own coefficients,
[`nl()`](https://statmodels7.github.io/modelterms7/reference/nl.md) and
[`seg()`](https://statmodels7.github.io/modelterms7/reference/seg.md),
this returns the block at the coefficients last committed by
[`term_refresh()`](https://statmodels7.github.io/modelterms7/reference/term_refresh.md).
[`term_jacobian_block()`](https://statmodels7.github.io/modelterms7/reference/term_jacobian_block.md)
says whether that block is a Jacobian or a frozen working linearization.

A structural term contributes no block and registers no method, so
`term_matrix()` on one stops with S7's method-not-found error.

## See also

[`term_predict()`](https://statmodels7.github.io/modelterms7/reference/term_predict.md)
for the same mapping at other rows,
[`term_coef_names()`](https://statmodels7.github.io/modelterms7/reference/term_coef_names.md)
for the column names,
[`term_npar()`](https://statmodels7.github.io/modelterms7/reference/term_npar.md)
for the count, and
[`term_refresh()`](https://statmodels7.github.io/modelterms7/reference/term_refresh.md)
for a block that moves.

## Examples

``` r
d <- data.frame(x = 1:4, g = factor(c("a", "b", "a", "b")))

term_matrix(term_build(linpar(~ x), d))
#>   (Intercept) x
#> 1           1 1
#> 2           1 2
#> 3           1 3
#> 4           1 4

# The column names are the term's coefficient names.
X <- term_matrix(term_build(linpar(~ x + g, label = "lin"), d))
colnames(X)
#> [1] "lin.(Intercept)" "lin.x"           "lin.gb"         

# A grouping indicator comes back sparse, not as a base matrix.
R <- term_matrix(term_build(random(~ 1 | g), d))
c(class = class(R)[1], is.matrix = is.matrix(R))
#>       class   is.matrix 
#> "dgCMatrix"     "FALSE" 

# A specification has nothing to return.
try(term_matrix(linpar(~ x)))
#> Error : the term has not been built; call term_build(term, data) first.
```
