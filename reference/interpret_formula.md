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
conventions. An [`offset()`](https://rdrr.io/r/stats/offset.html) term
is not among the labels that
[`stats::terms()`](https://rdrr.io/r/stats/terms.html) returns, so it is
dropped without a message;
[`statmodels7::statmod()`](https://statmodels7.github.io/statmodels7/reference/statmod.html)
separates offsets before calling this function. The terms are returned
unbuilt, so the caller chooses the data against which each term is
built.

## Usage

``` r
interpret_formula(formula, data, linpar = list())
```

## Arguments

- formula:

  A two-sided or one-sided model formula. Anything else signals the
  error `"'formula' must be a formula."`. The symbols of a term call are
  looked up in the formula's environment when `data` does not carry
  them, and that environment is carried onto the formula of the implicit
  `linpar` block.

- data:

  A data frame in which the formula's symbols are evaluated. Anything
  else signals the error `"'data' must be a data frame."`.

- linpar:

  A named list of arguments for the implicit
  [`linpar()`](https://statmodels7.github.io/modelterms7/reference/linpar.md)
  term, the term into which the bare covariates collapse; the arguments
  are `sparse` and `contrasts`. This is the only place where they can be
  given, because the caller never writes that term. The default is an
  empty list, which leaves
  [`linpar()`](https://statmodels7.github.io/modelterms7/reference/linpar.md)'s
  own defaults in place. A
  [`linpar()`](https://statmodels7.github.io/modelterms7/reference/linpar.md)
  term that the caller writes out takes its own arguments and ignores
  this one. Anything that is not a list signals an error.

## Value

A list with the elements:

- `response`:

  The evaluated left-hand side: a numeric vector, a
  [`censored_response()`](https://statmodels7.github.io/modelterms7/reference/censored_response.md),
  or `NULL` for a one-sided formula.

- `terms`:

  A named list of unbuilt term specifications. The collected parametric
  block comes first under the name `"linpar"` and is absent when the
  formula has no bare covariates and no intercept. Every other name is
  the term's label as it appears in the formula, deparsed, such as
  `"s(x2, basis7::bspline_smooth(k = 5))"`.

- `intercept`:

  `TRUE` unless the formula removes the intercept.

- `formula`:

  The input, unchanged.

## Recognition by evaluation

A term constructor is identified by what its call returns, so a term
class defined outside the package works in a formula without a list of
special names to amend. `log(x)` evaluates to a numeric vector and stays
a covariate; `s(x, basis7::bspline_smooth(k = 5))` evaluates to a
`SmoothTerm` and is routed as a term. A call whose evaluation signals an
error, such as a term constructor called with an argument that it does
not have, signals the error `"evaluating '<label>' failed: <message>"`,
so the call is not treated as a covariate.

Some labels are never evaluated: `:`, `*`, `^`, `%in%`, `+`, `-`, `(`
and `I`. A bare interaction `x1:x2` parses as a call to `:`, which on
numeric vectors is the sequence operator and would return something
unrelated to the column that
[`stats::model.matrix()`](https://rdrr.io/r/stats/model.matrix.html)
builds. Those labels go straight to the parametric block.

A call whose value is neither a term nor something that a model matrix
can hold signals an error that names the label and the class of the
value. The commonest cause is a masked name. The package `mgcv` also
exports
[`s()`](https://statmodels7.github.io/modelterms7/reference/s.md) and
[`te()`](https://statmodels7.github.io/modelterms7/reference/te.md), and
`segmented` exports
[`seg()`](https://statmodels7.github.io/modelterms7/reference/seg.md),
so where this package exports a term of that name and the name resolves
elsewhere, the message states where the function was found and suggests
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

A
[`linpar()`](https://statmodels7.github.io/modelterms7/reference/linpar.md)
written in a formula that carries an intercept drops its own column
`(Intercept)` after its model matrix is built, so its factors keep the
coding that they have with an intercept and `y ~ linpar(~ x + g)` gives
the columns of `lm(y ~ x + g)`. Under `0 +` it keeps the column.

## One covariate is removed

A [`seg()`](https://statmodels7.github.io/modelterms7/reference/seg.md)
or
[`jseg()`](https://statmodels7.github.io/modelterms7/reference/jseg.md)
term built with `linear = TRUE`, which is the default, contributes the
same column as the bare covariate, so `y ~ x + seg(x)` is rank deficient
by one. The term carries that effect, so the covariate is dropped from
the parametric block and a warning names both. The same holds for a
marginal
[`seg()`](https://statmodels7.github.io/modelterms7/reference/seg.md) or
[`jseg()`](https://statmodels7.github.io/modelterms7/reference/jseg.md)
built with `marginal = TRUE`. Writing `seg(x, linear = FALSE)` keeps the
linear effect outside the term.

Only the bare main effect is removed, matched on the deparsed
expression: an interaction spans no main effect and is left alone. A
term that spans the same direction without being that one column is
reported by a warning and left unmodified. A spline basis is the case:
it contains the line, and its penalty leaves the line unpenalized, so
`s(x) + seg(x)` is confounded too. There is no single column to remove
there, and this function does not modify another term.

## What it does not do

The function does not build any term: every element of `terms` is a
specification, and
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
for the block into which the covariates collapse,
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
#> <environment: 0x55bddcbf6060>

# A constructor call becomes a term, keyed by its label in the formula.
out2 <- interpret_formula(y ~ x1 + s(x2, basis7::bspline_smooth(k = 5)) + ridge(~ g), dd)
names(out2$terms)
#> [1] "linpar"                              
#> [2] "s(x2, basis7::bspline_smooth(k = 5))"
#> [3] "ridge(~g)"                           
vapply(out2$terms, function(t) class(t)[1], character(1))
#>                               linpar s(x2, basis7::bspline_smooth(k = 5)) 
#>            "modelterms7::LinparTerm"            "modelterms7::SmoothTerm" 
#>                            ridge(~g) 
#>         "modelterms7::PenalizedTerm" 

# The terms are not built yet.
vapply(out2$terms, term_is_built, logical(1))
#>                               linpar s(x2, basis7::bspline_smooth(k = 5)) 
#>                                FALSE                                FALSE 
#>                            ridge(~g) 
#>                                FALSE 

# The intercept convention is the formula's own.
names(interpret_formula(y ~ ridge(~ g), dd)$terms)
#> [1] "linpar"    "ridge(~g)"
names(interpret_formula(y ~ ridge(~ g) - 1, dd)$terms)
#> [1] "ridge(~g)"

# A linpar() written beside the formula's intercept drops its own.
li <- interpret_formula(y ~ linpar(~ x1 + g), dd)
term_coef_names(term_build(li$terms[["linpar(~x1 + g)"]], dd))
#> [1] "x1" "gb" "gc" "gd"

# An interaction is a covariate: `:` is never evaluated.
interpret_formula(y ~ x1:x2 + g, dd)$terms$linpar@formula
#> ~g + x1:x2
#> <environment: 0x55bddcbf6060>

# seg() carries the linear effect, so the bare covariate is removed.
w <- interpret_formula(y ~ x1 + seg(x1), dd)
#> Warning: the covariate 'x1' is exactly collinear with the linear effect that 'seg' carries, and has been removed from the parametric part. Write seg(x1, linear = FALSE) to keep the linear effect outside the term instead.
w$terms$linpar@formula          # x1 is gone; the term carries it
#> ~1
#> <environment: 0x55bddcbf6060>

# The bare covariate is kept when the term is told not to own the linear
# effect.
interpret_formula(y ~ x1 + seg(x1, linear = FALSE), dd)$terms$linpar@formula
#> ~x1
#> <environment: 0x55bddcbf6060>

# Arguments for the implicit block go through `linpar`.
sp <- interpret_formula(y ~ g, dd, linpar = list(sparse = TRUE))
class(term_matrix(term_build(sp$terms$linpar, dd)))
#> [1] "dgCMatrix"
#> attr(,"package")
#> [1] "Matrix"

# A call that returns neither a term nor a covariate is rejected by name.
e <- new.env(parent = globalenv())
assign("s", function(x, ...) structure(list(), class = "gamObject"), envir = e)
f <- y ~ s(x1)
environment(f) <- e
try(interpret_formula(f, dd))
#> Error : the term 's(x1)' evaluated to an object of class 'gamObject', which is neither a model term nor a covariate. The name 's' is masked here, and modelterms7 exports a term of that name: write modelterms7::s().
```
