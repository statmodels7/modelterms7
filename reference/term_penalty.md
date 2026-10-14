# Penalty of a Term

Returns the penalties7 penalty attached to the whole of a term's
coefficients, or `NULL` where there is none. The hyperparameters, their
bounds and links, and every derivative in the coefficients and in the
hyperparameters belong to the penalty object; the term only carries it.

## Usage

``` r
term_penalty(term, ...)
```

## Arguments

- term:

  An object inheriting from
  [`additive_term()`](https://statmodels7.github.io/modelterms7/reference/additive_term.md),
  built or not.

- ...:

  Passed to methods. The shipped methods do not use it.

## Value

A penalties7 penalty object, or `NULL` when the term is unpenalized,
when its penalty covers only part of its parameters, or when it has not
been built.

## A penalty over the whole block

A term whose penalty reaches part of its parameters returns `NULL` here
and declares that penalty through
[`term_penalties()`](https://statmodels7.github.io/modelterms7/reference/term_penalties.md),
which names the parameters it covers. A
[`seg()`](https://statmodels7.github.io/modelterms7/reference/seg.md)
term and an
[`nl()`](https://statmodels7.github.io/modelterms7/reference/nl.md) term
carry a penalty only where a subformula on one of their parameters
declares one, as in `gamma ~ 0 + lasso(~ 1)` over the changes of slope
of a
[`seg()`](https://statmodels7.github.io/modelterms7/reference/seg.md),
and their other parameters stay free;
[`gas()`](https://statmodels7.github.io/modelterms7/reference/gas.md)
declares its penalties through
[`term_penalties()`](https://statmodels7.github.io/modelterms7/reference/term_penalties.md)
in the same way. Reporting one of those penalties here would state that
it covers the whole block, so this generic answers only where that
holds.

[`term_penalties()`](https://statmodels7.github.io/modelterms7/reference/term_penalties.md)
is therefore the general function, and a fitting layer enumerates
penalties through it. This one is the shortcut for the common case.

## A specification carries no penalty

The penalty is attached at
[`term_build()`](https://statmodels7.github.io/modelterms7/reference/term_build.md),
its width being the number of columns the data produce, so
`term_penalty(ridge(~ x))` is `NULL` and
`term_penalty(term_build(ridge(~ x), d))` is the quadratic penalty. The
same holds for
[`term_penalties()`](https://statmodels7.github.io/modelterms7/reference/term_penalties.md)
and, through it, for
[`term_smooth()`](https://statmodels7.github.io/modelterms7/reference/term_smooth.md).

The one method is registered on
[`additive_term()`](https://statmodels7.github.io/modelterms7/reference/additive_term.md)
and reads the `penalty` property. A structural term has no such
property, so `term_penalty()` on one signals S7's method-not-found
error.

## See also

[`term_penalties()`](https://statmodels7.github.io/modelterms7/reference/term_penalties.md)
for the general form,
[`term_smooth()`](https://statmodels7.github.io/modelterms7/reference/term_smooth.md)
for whether the result has a kink,
[`edf()`](https://statmodels7.github.io/modelterms7/reference/edf.md)
for what it costs, and
[`penalties7::penalty_value()`](https://statmodels7.github.io/penalties7/reference/penalty_value.html)
for what the returned object computes.

## Examples

``` r
d <- data.frame(x = rnorm(20), g = factor(rep(c("a", "b"), 10)))

# Unpenalized, and unbuilt: both give NULL.
term_penalty(linpar(~ x))
#> NULL
term_penalty(ridge(~ x))
#> NULL

# Built, it is the penalty object itself.
p <- term_penalty(term_build(ridge(~ x), d))
c(name = p@penalty_name, params = p@params, n_coef = p@n_coef)
#>        name      params      n_coef 
#> "quadratic"    "lambda"         "1" 

# A term whose penalty covers part of its parameters answers NULL here
# and declares it through term_penalties() instead: a lasso on the
# change of slope of a break-point term.
sb <- term_build(seg(x, npsi = 1, gamma ~ 0 + lasso(~ 1)),
                 data.frame(x = sort(runif(50, 0, 10))))
term_penalty(sb)
#> NULL
vapply(term_penalties(sb), function(z) z$name, character(1))
#> [1] "gamma1::lasso(~1)"
```
