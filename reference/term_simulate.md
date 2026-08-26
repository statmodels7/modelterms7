# Drawing a Response From a Structural Term

The predictor a structural term produces when the response is being
generated rather than read, together with whatever latent quantity the
term drew on the way.

## Usage

``` r
term_simulate(term, psi, eta, draw, ...)
```

## Arguments

- term:

  A built structural term.

- psi:

  The term's parameters, named as
  [`term_params()`](https://statmodels7.github.io/modelterms7/reference/term_params.md).

- eta:

  The static part of the predictor, one value per observation.

- draw:

  A function `(e, i)` returning one response value drawn at predictor
  `e` for observation `i`.

- ...:

  Passed to methods.

## Value

A list with `eta`, the predictor; `y`, the responses drawn or `NULL`;
and `latent`, whatever the term drew.

## Details

Simulating from a model that carries state is not the same operation as
fitting one, and the difference is which direction the response moves
in. A term whose contribution does not read the response can report that
contribution and leave the drawing to the caller: a latent chain's
levels and a group's break-point drawn from its prior are both like
that. A score-driven term cannot. Its level at one time is driven by the
score of the response at the time before, so the response has to be
drawn AS the recursion runs.

One contract covers both. The caller supplies `draw`, a function of a
predictor and a row index returning one response value, and the method
returns the predictor it produced; a method that drew returns the
responses as well and one that did not returns `NULL` there, leaving the
caller to draw at the predictor.

## See also

[`term_filter()`](https://statmodels7.github.io/modelterms7/reference/term_filter.md),
[`term_continue()`](https://statmodels7.github.io/modelterms7/reference/term_continue.md)

## Examples

``` r
set.seed(1)
dd <- data.frame(t = 1:30)
term <- term_build(gas(p = 1, q = 1, time = t), dd)
out <- term_simulate(term, list(omega = 0.5, alpha1 = 0.3, pacf1 = 0.6),
                     rep(0, 30),
                     draw = function(e, i) stats::rnorm(1, e, 1))
head(out$y, 3)
#> [1] 0.6235462 1.2457072 0.3567027
```
