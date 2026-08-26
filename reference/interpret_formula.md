# Interpret a Model Formula Into Terms

Splits a model formula into a response and a named list of term
specifications. A call on the right-hand side is evaluated, and if its
value inherits from
[`model_term()`](https://statmodels7.github.io/modelterms7/reference/model_term.md)
it becomes a term of its own; everything else, meaning bare covariates,
transformations such as `log(x)` and interactions, is collected into a
single
[`linpar()`](https://statmodels7.github.io/modelterms7/reference/linpar.md)
block with the usual
[`stats::model.matrix()`](https://rdrr.io/r/stats/model.matrix.html)
conventions. The terms come back unbuilt, so the caller decides which
data each is built against.

## Usage

``` r
interpret_formula(formula, data, linpar = list())
```

## Arguments

- formula:

  A two-sided or one-sided model formula. Anything else throws
  `"'formula' must be a formula."`. Its environment is where a term
  call's symbols are looked up when `data` does not carry them, and it
  is carried onto the formula of the implicit `linpar` block.

- data:

  A data frame in which the formula's symbols are evaluated. Anything
  else throws `"'data' must be a data frame."`.

- linpar:

  A named list of arguments for the **implicit**
  [`linpar()`](https://statmodels7.github.io/modelterms7/reference/linpar.md)
  term, the one the bare covariates collapse into: `sparse` and
  `contrasts`. This is the only place they can be given, that term never
  being written by the caller. Empty by default, which leaves
  [`linpar()`](https://statmodels7.github.io/modelterms7/reference/linpar.md)'s
  own defaults. A
  [`linpar()`](https://statmodels7.github.io/modelterms7/reference/linpar.md)
  term the caller writes out takes its own arguments and ignores this
  one. Anything that is not a list throws.

## Value

A list of four elements:

- `response`:

  The evaluated left-hand side: a numeric vector, a
  [`censored_response()`](https://statmodels7.github.io/modelterms7/reference/censored_response.md),
  or `NULL` for a one-sided formula.

- `terms`:

  A named list of unbuilt term specifications. The collected parametric
  block comes first under the name `"linpar"` and is absent when the
  formula has no bare covariates and no intercept. Every other name is
  the term's label as it appears in the formula, deparsed, such as
  `"s(x2, k = 5)"`.

- `intercept`:

  `TRUE` unless the formula removes the intercept.

- `formula`:

  The input, unchanged.

## Recognition by evaluation

A term constructor is identified by what its call returns, so a term
class defined outside the package works in a formula the day it is
written, with no list of special names to amend. `log(x)` evaluates to a
numeric vector and stays a covariate; `s(x, k = 5)` evaluates to a
`SmoothTerm` and is routed as a term.

Some labels are never evaluated: `:`, `*`, `^`, `%in%`, `+`, `-`, `(`
and `I`. A bare interaction `x1:x2` parses as a call to `:`, which on
numeric vectors is the sequence operator and would return something
unrelated to the column
[`stats::model.matrix()`](https://rdrr.io/r/stats/model.matrix.html)
builds. Those labels go straight to the parametric block.

A call whose value is neither a term nor something a model matrix can
hold throws, naming the label and the class it produced. The commonest
way to arrive there is a masked name. `mgcv` also exports
[`s()`](https://statmodels7.github.io/modelterms7/reference/s.md) and
[`te()`](https://statmodels7.github.io/modelterms7/reference/te.md) and
`segmented` exports
[`seg()`](https://statmodels7.github.io/modelterms7/reference/seg.md),
so where the package exports a term of that name and the name resolves
elsewhere, the message says where it was found and suggests
`modelterms7::`.

## The response and the intercept

The left-hand side, when there is one, is evaluated in `data`: a plain
expression gives a numeric vector, and a response constructor such as
[`cens()`](https://statmodels7.github.io/modelterms7/reference/cens.md)
gives its object. A one-sided formula gives `response = NULL`.

The intercept convention is the formula's own, carried into the
collected parametric block. `y ~ ridge(~ g)` still produces an
intercept-only `linpar` block, and `y ~ ridge(~ g) - 1` produces no
`linpar` block at all.

## One covariate is removed

A [`seg()`](https://statmodels7.github.io/modelterms7/reference/seg.md)
or
[`jseg()`](https://statmodels7.github.io/modelterms7/reference/jseg.md)
term built with `linear = TRUE`, which is the default, contributes the
same column the bare covariate would, so `y ~ x + seg(x)` is rank
deficient by one. The term owns that effect, so the covariate is dropped
from the parametric block and a warning names both. Writing
`seg(x, linear = FALSE)` keeps the linear effect outside the term.

Only the bare main effect is removed, matched on the deparsed
expression: an interaction spans no main effect and is left alone. A
term that spans the same direction without being that one column is
reported by a warning and left unmodified. A spline basis is the case:
it contains the line, and its penalty leaves the line unpenalized, so
`s(x) + seg(x)` is confounded too. There is no single column to remove
there, and reshaping another term is not this function's business.

## What it does not do

Nothing is built: every element of `terms` is a specification, and
[`term_is_built()`](https://statmodels7.github.io/modelterms7/reference/term_is_built.md)
is `FALSE` for all of them. `data` is read for
[`stats::terms()`](https://rdrr.io/r/stats/terms.html)'s variable
classification and for evaluating the term calls and the response, and
the term constructors themselves do not touch it.

Labels are deduplicated by
[`stats::terms()`](https://rdrr.io/r/stats/terms.html) before any of
this, so `y ~ ridge(~ x) + ridge(~ x)` gives one term, as `y ~ x + x`
gives one column.

## See also

[`linpar()`](https://statmodels7.github.io/modelterms7/reference/linpar.md)
for the block the covariates collapse into,
[`cens()`](https://statmodels7.github.io/modelterms7/reference/cens.md)
for the response constructor,
[`check_term()`](https://statmodels7.github.io/modelterms7/reference/check_term.md)
for validating one of the returned specifications, and
[`term_build()`](https://statmodels7.github.io/modelterms7/reference/term_build.md)
for building it.

## Examples

``` r
dd <- data.frame(y = rnorm(20), x1 = 1:20, x2 = runif(20),
                 g = factor(rep(letters[1:4], 5)))

# Covariates and transformations collapse into one parametric block.
out <- interpret_formula(y ~ x1 + log(x2), dd)
names(out$terms)
#> [1] "linpar"
out$terms$linpar@formula
#> ~x1 + log(x2)
#> <environment: 0x55d8f1736d70>

# A constructor call becomes a term, keyed by its label in the formula.
out2 <- interpret_formula(y ~ x1 + s(x2, k = 5) + ridge(~ g), dd)
names(out2$terms)
#> [1] "linpar"       "s(x2, k = 5)" "ridge(~g)"   
vapply(out2$terms, function(t) class(t)[1], character(1))
#>                       linpar                 s(x2, k = 5) 
#>    "modelterms7::LinparTerm"    "modelterms7::SmoothTerm" 
#>                    ridge(~g) 
#> "modelterms7::PenalizedTerm" 

# Nothing is built yet.
vapply(out2$terms, term_is_built, logical(1))
#>       linpar s(x2, k = 5)    ridge(~g) 
#>        FALSE        FALSE        FALSE 

# The intercept convention is the formula's own.
names(interpret_formula(y ~ ridge(~ g), dd)$terms)
#> [1] "linpar"    "ridge(~g)"
names(interpret_formula(y ~ ridge(~ g) - 1, dd)$terms)
#> [1] "ridge(~g)"

# An interaction is a covariate: `:` is never evaluated.
interpret_formula(y ~ x1:x2 + g, dd)$terms$linpar@formula
#> ~g + x1:x2
#> <environment: 0x55d8f1736d70>

# seg() carries the linear effect, so the bare covariate is removed.
w <- interpret_formula(y ~ x1 + seg(x1), dd)
#> Warning: the covariate 'x1' is exactly collinear with the linear effect that 'seg' carries, and has been removed from the parametric part. Write seg(x1, linear = FALSE) to keep the linear effect outside the term instead.
w$terms$linpar@formula          # x1 is gone; the term carries it
#> ~1
#> <environment: 0x55d8f1736d70>

# Unless the term is told not to own it.
interpret_formula(y ~ x1 + seg(x1, linear = FALSE), dd)$terms$linpar@formula
#> ~x1
#> <environment: 0x55d8f1736d70>

# Arguments for the implicit block go through `linpar`.
sp <- interpret_formula(y ~ g, dd, linpar = list(sparse = TRUE))
class(term_matrix(term_build(sp$terms$linpar, dd)))
#> [1] "dgCMatrix"
#> attr(,"package")
#> [1] "Matrix"

# A call that returns neither a term nor a covariate is refused by name.
e <- new.env(parent = globalenv())
assign("s", function(x, ...) structure(list(), class = "gamObject"), envir = e)
f <- y ~ s(x1)
environment(f) <- e
try(interpret_formula(f, dd))
#> Error : the term 's(x1)' evaluated to an object of class 'gamObject', which is neither a model term nor a covariate. The name 's' is masked here, and modelterms7 exports a term of that name: write modelterms7::s().
```
