# The Hyperparameters a Term Holds

Reports the hyperparameter values the caller fixed in the constructor,
one entry per penalty the term carries. A hyperparameter left `NULL`, or
written out as a grid, is **estimated** by the fit and does not appear
here, so an empty result means that the term holds none.

## Usage

``` r
term_hyper(term, ...)
```

## Arguments

- term:

  A term, built or not.

- ...:

  Passed to methods. The shipped methods do not use it.

## Value

A named list, one entry per penalty of the term and named as
[`term_penalties()`](https://statmodels7.github.io/modelterms7/reference/term_penalties.md)
names it, each entry a named list of single numbers keyed by the
penalty's own hyperparameter names. Empty where the term holds nothing.

## The three states of one argument

A constructor's hyperparameter argument carries three meanings, read for
one hyperparameter at a time. `NULL` leaves the hyperparameter to be
estimated: a path builds a grid for a penalty with a kink, and a
marginal criterion estimates the value at the mode for a smooth one.
`NULL` is the default for `lambda` and `alpha`, whereas the shape of
[`scad()`](https://statmodels7.github.io/modelterms7/reference/scad.md)
and
[`mcp()`](https://statmodels7.github.io/modelterms7/reference/mcp.md)
defaults to a held number, and `NULL` there has the shape estimated.
**One number** holds the hyperparameter, and that is what this function
reports. **Several numbers** are a grid that the path visits as they
stand, which
[`term_values()`](https://statmodels7.github.io/modelterms7/reference/term_values.md)
reports; the hyperparameter is still estimated, the values to search
being given and the estimate being left to the fit.

## How the entries are keyed

The names are the ones
[`term_penalties()`](https://statmodels7.github.io/modelterms7/reference/term_penalties.md)
gives its entries: `""` for a penalty covering the whole term, and
`parameter::label` for a penalty brought in by a subformula. A term
whose penalties come from sub-terms, such as
[`nl()`](https://statmodels7.github.io/modelterms7/reference/nl.md) or
[`seg()`](https://statmodels7.github.io/modelterms7/reference/seg.md)
with subformulas, propagates their held values without a method of its
own, every entry of that enumeration carrying its own in the field
`fixed`.

A built term reports those entries, and an unbuilt term reports its own
`hyper` property with the names as written. The two answers are the same
for the penalized terms, and for
[`s()`](https://statmodels7.github.io/modelterms7/reference/s.md) and
[`random()`](https://statmodels7.github.io/modelterms7/reference/random.md)
where the names are not qualified. A qualified name, such as `lambda.a`
for the level `a` of a factor, is reported as written before the build
and under its entry after it. For a term whose held values sit in
subformulas, such as
[`nl()`](https://statmodels7.github.io/modelterms7/reference/nl.md) and
[`seg()`](https://statmodels7.github.io/modelterms7/reference/seg.md),
the sub-terms exist only after the build, so an unbuilt term reports an
empty list.

## See also

[`term_values()`](https://statmodels7.github.io/modelterms7/reference/term_values.md)
for a written-out grid,
[`term_grid()`](https://statmodels7.github.io/modelterms7/reference/term_grid.md)
for how many values a built grid visits,
[`term_penalties()`](https://statmodels7.github.io/modelterms7/reference/term_penalties.md)
for the entries this is keyed by.

## Examples

``` r
# Held, and free.
term_hyper(lasso(~ x, lambda = 3))
#> [[1]]
#> [[1]]$lambda
#> [1] 3
#> 
#> 
term_hyper(lasso(~ x))
#> list()

# Several values are a grid, not a held value, so they go elsewhere.
term_hyper(lasso(~ x, lambda = c(0.1, 1, 10)))
#> list()
term_values(lasso(~ x, lambda = c(0.1, 1, 10)))
#> [[1]]
#> [[1]]$lambda
#> [1]  0.1  1.0 10.0
#> 
#> 

# One of two held, the other estimated.
term_hyper(enet(~ x, alpha = 0.5))
#> [[1]]
#> [[1]]$alpha
#> [1] 0.5
#> 
#> 

# The shape of SCAD is held by default; NULL has it estimated.
term_hyper(scad(~ x))
#> [[1]]
#> [[1]]$a
#> [1] 3.7
#> 
#> 
term_hyper(scad(~ x, a = NULL))
#> list()

# A held value on a sub-term is reported under that entry's name once the
# term is built.
set.seed(4)
d <- data.frame(x = runif(40, 0, 5), g = factor(rep(1:4, 10)))
nu <- nl(~ a * exp(-r * x), a ~ 0 + lasso(~ g, lambda = 2), r ~ 1)
term_hyper(nu)
#> list()
term_hyper(term_build(nu, d))
#> $`a::lasso(~g, lambda = 2)`
#> $`a::lasso(~g, lambda = 2)`$lambda
#> [1] 2
#> 
#> 
```
