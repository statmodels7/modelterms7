# Interpret a Model Formula Into Terms

Splits a model formula into a response specification and a list of term
specifications. Term constructors are recognized by what they return: a
call on the right-hand side whose value inherits from
[`model_term`](https://statmodels7.github.io/modelterms7/reference/model_term.md)
becomes a term of its own, and everything else (bare covariates,
transformations such as `log(x)`, interactions) is collected into one
[`linpar`](https://statmodels7.github.io/modelterms7/reference/linpar.md)
block with the usual
[`model.matrix`](https://rdrr.io/r/stats/model.matrix.html) conventions.

## Usage

``` r
interpret_formula(formula, data)
```

## Arguments

- formula:

  A model formula.

- data:

  A data frame in which the formula's symbols are evaluated.

## Value

A list with elements `response` (the evaluated left-hand side, or `NULL`
for a one-sided formula), `terms` (a named list of term specifications,
the collected parametric block first under the name `"linpar"`),
`intercept` (logical) and `formula` (the input).

## Details

Recognition by evaluation is what makes the interpreter extensible: a
term class defined outside the package works in a formula the day it is
written, with no list of special names to amend. `log(x)` evaluates to a
numeric vector and stays a covariate; a constructor call evaluates to a
term specification and is routed as one. Interaction labels and bare
symbols are never evaluated directly.

The left-hand side, when present, is evaluated in the data: a plain
expression gives a numeric response, and a response constructor such as
[`cens`](https://statmodels7.github.io/modelterms7/reference/cens.md)
gives its response object. The intercept convention is the formula's
own, carried into the collected parametric block, so `y ~ ridge_like(R)`
still produces an intercept-only `linpar` block and
`y ~ ridge_like(R) - 1` produces none.

## Examples

``` r
dd <- data.frame(y = rnorm(6), x1 = 1:6, x2 = runif(6))
out <- interpret_formula(y ~ x1 + log(x2), dd)
names(out$terms)
#> [1] "linpar"
```
