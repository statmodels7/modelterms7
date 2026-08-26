# Build a Term on Data

Turns a term specification into a built term: an additive term computes
its design block from `data`, assigns the coefficient names and records
the blueprint that will reproduce the mapping on other rows; a
structural term records whatever its recursion needs, its grouping and
its ordering. The returned object is a copy of the specification with
those properties filled, and the specification is unchanged.

## Usage

``` r
term_build(term, data, ...)
```

## Arguments

- term:

  An object inheriting from
  [`model_term()`](https://statmodels7.github.io/modelterms7/reference/model_term.md),
  built or not. A built term is rebuilt.

- data:

  A data frame carrying every variable the term names. Anything else
  throws `"'data' must be a data frame."` from the generic, before
  dispatch.

- ...:

  Passed to methods. No shipped method reads anything here.

## Value

A built term of the same class as `term`. For an additive term the `X`,
`coef_names` and `blueprint` properties are filled; for a structural
term the `blueprint` is.
[`term_is_built()`](https://statmodels7.github.io/modelterms7/reference/term_is_built.md)
is `TRUE` either way, reading whichever of the two the branch fills.

## What building produces

An additive term contributes to the linear predictor through a block,
and through a penalty on that block's coefficients when it is penalized:

\$\$\eta = \sum\_{t} X_t \beta_t, \qquad \text{penalized objective}
\quad -\ell(\beta) + \sum\_{t} \rho_t(\beta_t; \theta_t).\$\$

Building is what produces \\X_t\\ from the data and attaches \\\rho_t\\.
[`term_matrix()`](https://statmodels7.github.io/modelterms7/reference/term_matrix.md)
then reads the block,
[`term_penalties()`](https://statmodels7.github.io/modelterms7/reference/term_penalties.md)
the penalties, and
[`term_predict()`](https://statmodels7.github.io/modelterms7/reference/term_predict.md)
reproduces \\X_t\\ at other rows through the blueprint.

A structural term has no such block.
[`gas()`](https://statmodels7.github.io/modelterms7/reference/gas.md)
records the group each row belongs to and its place in that group's
series;
[`regime()`](https://statmodels7.github.io/modelterms7/reference/regime.md)
and the marginal break-point terms record what their forward recursion
reads. They report themselves through
[`term_filter()`](https://statmodels7.github.io/modelterms7/reference/term_filter.md)
or
[`term_loglik()`](https://statmodels7.github.io/modelterms7/reference/term_loglik.md),
and
[`term_matrix()`](https://statmodels7.github.io/modelterms7/reference/term_matrix.md)
has no method for them.

## Building twice, and building on other data

Building is not idempotent in general: a term built again on new data
re-derives its factor levels and its knots from those rows. That is what
[`term_predict()`](https://statmodels7.github.io/modelterms7/reference/term_predict.md)
exists to avoid, and what
[`check_term()`](https://statmodels7.github.io/modelterms7/reference/check_term.md)'s
subset check tests for. Build once, predict thereafter.

## The default

There is one, on
[`model_term()`](https://statmodels7.github.io/modelterms7/reference/model_term.md),
and it throws `"the term class 'X' does not implement term_build()."`,
naming the class, so a term class that supplies nothing else says so
clearly. It covers both branches: every shipped structural term
registers a method of its own, and a structural class written elsewhere
that does not is named the same way.

## See also

[`term_predict()`](https://statmodels7.github.io/modelterms7/reference/term_predict.md)
for the block at other rows,
[`term_matrix()`](https://statmodels7.github.io/modelterms7/reference/term_matrix.md)
and
[`term_coef_names()`](https://statmodels7.github.io/modelterms7/reference/term_coef_names.md)
for what a build filled,
[`term_refresh()`](https://statmodels7.github.io/modelterms7/reference/term_refresh.md)
for a block that moves with its coefficients, and
[`check_term()`](https://statmodels7.github.io/modelterms7/reference/check_term.md)
for validating the result.

## Examples

``` r
d <- data.frame(x = rnorm(20), g = factor(rep(letters[1:4], 5)))

# A specification carries no block; building fills it.
spec <- linpar(~ x + g)
built <- term_build(spec, d)
c(spec = term_is_built(spec), built = term_is_built(built))
#>  spec built 
#> FALSE  TRUE 
dim(term_matrix(built))
#> [1] 20  5
term_coef_names(built)
#> [1] "(Intercept)" "x"           "gb"          "gc"          "gd"         

# The specification is untouched: building returns a copy.
term_is_built(spec)
#> [1] FALSE

# A structural term records its recursion's bookkeeping and no block.
g <- term_build(gas(p = 1, q = 1), data.frame(y = rnorm(30)))
term_params(g)
#> [1] "omega"  "alpha1" "pacf1" 
try(term_matrix(g))
#> Error : Can't find method for `term_matrix(<modelterms7::GasTerm>)`.

# A class that implements nothing is told which class it is.
Foo <- S7::new_class("Foo", parent = additive_term)
try(term_build(Foo(), d))
#> Error : the term class 'Foo' does not implement term_build().
```
