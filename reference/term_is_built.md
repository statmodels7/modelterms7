# Whether a Term Has Been Built

Returns `TRUE` for a term that
[`term_build()`](https://statmodels7.github.io/modelterms7/reference/term_build.md)
has filled and `FALSE` for a bare specification.
[`term_matrix()`](https://statmodels7.github.io/modelterms7/reference/term_matrix.md),
[`term_npar()`](https://statmodels7.github.io/modelterms7/reference/term_npar.md),
[`term_coef_names()`](https://statmodels7.github.io/modelterms7/reference/term_coef_names.md),
[`term_predict()`](https://statmodels7.github.io/modelterms7/reference/term_predict.md)
and [`plot()`](https://rdrr.io/r/graphics/plot.default.html) apply this
test before reading a design block, so the result determines which error
a caller receives from those functions.

## Usage

``` r
term_is_built(term)
```

## Arguments

- term:

  An object inheriting from
  [`model_term()`](https://statmodels7.github.io/modelterms7/reference/model_term.md).
  Any other object is rejected with the error
  `"'term' must inherit from 'model_term'."`.

## Value

A single logical, never `NA`.

## Details

The two branches record being built in different places, so the
predicate reads the record of each branch. An additive term is built
when it has coefficient names, `length(term@coef_names) > 0L`; a
structural term contributes no design columns and so has none to count,
and is built when its blueprint is filled,
`length(term@blueprint) > 0L`.

A built additive term with no columns would give `FALSE`. The shipped
constructors do not produce one: a formula that gives no columns, such
as `~ 0`, is rejected when
[`term_build()`](https://statmodels7.github.io/modelterms7/reference/term_build.md)
runs.

The `blueprint` property is read directly. It is declared on
[`additive_term()`](https://statmodels7.github.io/modelterms7/reference/additive_term.md)
and on the abstract
[`structural_term()`](https://statmodels7.github.io/modelterms7/reference/structural_term.md),
which the three shipped structural classes inherit, so every term
inheriting from
[`model_term()`](https://statmodels7.github.io/modelterms7/reference/model_term.md)
has one. A structural class written outside the package inherits an
empty blueprint, which reads as not built until
[`term_build()`](https://statmodels7.github.io/modelterms7/reference/term_build.md)
fills it, and the predicate returns a logical for such a term as well.

## See also

[`term_build()`](https://statmodels7.github.io/modelterms7/reference/term_build.md),
which makes it `TRUE`;
[`term_matrix()`](https://statmodels7.github.io/modelterms7/reference/term_matrix.md)
and
[`term_npar()`](https://statmodels7.github.io/modelterms7/reference/term_npar.md),
which reject a term for which it is `FALSE`.

[`term_build()`](https://statmodels7.github.io/modelterms7/reference/term_build.md),
[`term_matrix()`](https://statmodels7.github.io/modelterms7/reference/term_matrix.md),
[`term_coef_names()`](https://statmodels7.github.io/modelterms7/reference/term_coef_names.md),
[`term_npar()`](https://statmodels7.github.io/modelterms7/reference/term_npar.md)

## Examples

``` r
d <- data.frame(x = 1:4)
term_is_built(linpar(~ x))
#> [1] FALSE
term_is_built(term_build(linpar(~ x), d))
#> [1] TRUE

# The accessors apply this test, so it predicts their error.
try(term_matrix(linpar(~ x)))
#> Error : the term has not been built; call term_build(term, data) first.

# A structural term carries no design block, so the predicate reads its
# blueprint instead.
gs <- gas(p = 1, q = 1)
gb <- term_build(gs, data.frame(y = rnorm(20)))
c(spec = term_is_built(gs), built = term_is_built(gb))
#>  spec built 
#> FALSE  TRUE 
```
