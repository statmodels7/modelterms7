# Build a Random-Effect Term

Builds the within-group design from the left of the bar, interacts it
with the group indicators, and attaches the effects' distribution as the
penalty on the resulting coefficients. The levels of the grouping
variable are recorded, so
[`term_predict()`](https://statmodels7.github.io/modelterms7/reference/term_predict.md)
maps new rows onto the same ones.

## Arguments

- term:

  An unbuilt or built
  [`RandomTerm()`](https://statmodels7.github.io/modelterms7/reference/RandomTerm.md).

- data:

  A data frame carrying the grouping variable and the within-group
  covariates.

- ...:

  Unused.

## Value

The term with `X` (a `dgCMatrix` of \\md\\ columns), `coef_names`,
`blueprint` and `penalty` filled.

## The block

With \\m\\ levels and a within-group design \\Z_i\\ of \\d\\ columns,
the block is \\\mathrm{diag}(Z_1, \dots, Z_m)\\, ordered
**group-major**, so the \\d\\ coefficients of one group are adjacent. It
is built as a `dgCMatrix`: a row belongs to one group, so the density is
\\1/m\\.

The coefficient names are `label.level.column`, so `random(~ x | g)`
over three levels gives `random.a.(Intercept)`, `random.a.x`,
`random.b.(Intercept)` and so on, which is the group-major order read
off.

## The penalty, and what the build checks

Where `distrib` is `NULL` the default is chosen here: a centered
`gaussian1_distrib` at one column or under `correlated = FALSE`, and a
centered multivariate Gaussian on an unstructured covariance for several
correlated ones.

Whatever the distribution, its location parameters must be **held**. A
free location is confounded with the intercept of the equation the term
sits in, and the build rejects it with a message naming the parameter
and the fix. A multivariate distribution must also match the
within-group dimension, and one that carries no matrix parameter cannot
express correlated effects at all.

Any value named in `hyper` is checked here too, against the penalty's
own names, this being the first point at which the penalty exists.

## See also

[`random()`](https://statmodels7.github.io/modelterms7/reference/random.md),
[`term_predict.RandomTerm()`](https://statmodels7.github.io/modelterms7/reference/term_predict.RandomTerm.md),
[`term_penalties()`](https://statmodels7.github.io/modelterms7/reference/term_penalties.md).

## Examples

``` r
set.seed(1)
dd <- data.frame(x = rnorm(9), g = factor(rep(c("a", "b", "c"), 3)))

# Group-major: the two coefficients of one level are adjacent.
b <- term_build(random(~ x | g), dd)
term_coef_names(b)
#> [1] "random.a.(Intercept)" "random.a.x"           "random.b.(Intercept)"
#> [4] "random.b.x"           "random.c.(Intercept)" "random.c.x"          
as.matrix(term_matrix(b))
#>       random.a.(Intercept) random.a.x random.b.(Intercept) random.b.x
#>  [1,]                    1 -0.6264538                    0  0.0000000
#>  [2,]                    0  0.0000000                    1  0.1836433
#>  [3,]                    0  0.0000000                    0  0.0000000
#>  [4,]                    1  1.5952808                    0  0.0000000
#>  [5,]                    0  0.0000000                    1  0.3295078
#>  [6,]                    0  0.0000000                    0  0.0000000
#>  [7,]                    1  0.4874291                    0  0.0000000
#>  [8,]                    0  0.0000000                    1  0.7383247
#>  [9,]                    0  0.0000000                    0  0.0000000
#>       random.c.(Intercept) random.c.x
#>  [1,]                    0  0.0000000
#>  [2,]                    0  0.0000000
#>  [3,]                    1 -0.8356286
#>  [4,]                    0  0.0000000
#>  [5,]                    0  0.0000000
#>  [6,]                    1 -0.8204684
#>  [7,]                    0  0.0000000
#>  [8,]                    0  0.0000000
#>  [9,]                    1  0.5757814

# A free location is refused, naming the parameter.
try(term_build(random(~ 1 | g,
                      distrib = distributions7::gaussian1_distrib()), dd))
#> Error : 'distrib' has a free location ('mu'), which is confounded with the
#>   intercept of the equation the term sits in. Hold it with
#>   distributions7::fixed(): at zero for a family on the line, and at
#>   whatever centers the effects where the prior is a transformation of
#>   another family, the parameter being the mean on the original scale.
```
