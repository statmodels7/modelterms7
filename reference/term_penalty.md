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

  Passed to methods. No shipped method reads anything here.

## Value

A penalties7 penalty object, or `NULL` when the term is unpenalized,
when its penalty covers only part of its parameters, or when it has not
been built.

## It answers only for a penalty over the whole block

A term whose penalty reaches part of its parameters returns `NULL` here
and declares that penalty through
[`term_penalties()`](https://statmodels7.github.io/modelterms7/reference/term_penalties.md),
which names the parameters it covers.
[`seg()`](https://statmodels7.github.io/modelterms7/reference/seg.md)
penalizes its changes of slope and leaves the linear effect and the
break-points free;
[`nl()`](https://statmodels7.github.io/modelterms7/reference/nl.md) and
[`gas()`](https://statmodels7.github.io/modelterms7/reference/gas.md)
carry the penalties of whatever sub-terms develop their own parameters.
Reporting one of those here would say it covers the block, so this
generic answers only where that is true.

[`term_penalties()`](https://statmodels7.github.io/modelterms7/reference/term_penalties.md)
is therefore the general question, and it is the one a fitting layer
asks. This one is the convenience for the common case.

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
property, so `term_penalty()` on one stops with S7's method-not-found
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
# and declares it through term_penalties() instead.
sb <- term_build(seg(x, npsi = 1), data.frame(x = sort(runif(50, 0, 10))))
term_penalty(sb)
#> NULL
```
