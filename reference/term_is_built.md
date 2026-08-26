# Whether a Term Has Been Built

`TRUE` for a term that
[`term_build()`](https://statmodels7.github.io/modelterms7/reference/term_build.md)
has filled, `FALSE` for a bare specification. It is the test
[`term_matrix()`](https://statmodels7.github.io/modelterms7/reference/term_matrix.md),
[`term_npar()`](https://statmodels7.github.io/modelterms7/reference/term_npar.md),
[`term_coef_names()`](https://statmodels7.github.io/modelterms7/reference/term_coef_names.md),
[`term_predict()`](https://statmodels7.github.io/modelterms7/reference/term_predict.md)
and [`plot()`](https://rdrr.io/r/graphics/plot.default.html) apply
before reading a design block, so it decides which error a caller gets
from those.

## Usage

``` r
term_is_built(term)
```

## Arguments

- term:

  An object inheriting from
  [`model_term()`](https://statmodels7.github.io/modelterms7/reference/model_term.md).
  Anything else throws `"'term' must inherit from 'model_term'."`.

## Value

A single logical, never `NA`.

## Details

The two branches record being built in different places, so the
predicate asks each about its own. An additive term is built when it has
coefficient names, `length(term@coef_names) > 0L`; a structural term
contributes no design columns and so has none to count, and is built
when its blueprint is filled, `length(term@blueprint) > 0L`.

**A built additive term with no columns would answer `FALSE`.** No
shipped constructor produces one:
[`linpar()`](https://statmodels7.github.io/modelterms7/reference/linpar.md)
with an empty formula fails in the class validator before it gets here.

**`blueprint` is asked for rather than assumed.** It is declared on
[`additive_term()`](https://statmodels7.github.io/modelterms7/reference/additive_term.md)
and on each of the three shipped structural classes, and not on the
abstract
[`structural_term()`](https://statmodels7.github.io/modelterms7/reference/structural_term.md),
so a structural class written outside the package need not carry one.
Where it does not, the answer is `FALSE` rather than an error, this
predicate promising a logical for anything inheriting from
[`model_term()`](https://statmodels7.github.io/modelterms7/reference/model_term.md).

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

# It is what the accessors test, so it predicts the error.
try(term_matrix(linpar(~ x)))
#> Error : the term has not been built; call term_build(term, data) first.

# A structural term carries no design block, so it is asked about the
# blueprint instead, and the answer is the same question either way.
gs <- gas(p = 1, q = 1)
gb <- term_build(gs, data.frame(y = rnorm(20)))
c(spec = term_is_built(gs), built = term_is_built(gb))
#>  spec built 
#> FALSE  TRUE 
```
