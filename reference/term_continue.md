# Continuing a Structural Term Past the Observed Series

The contribution a structural term makes at rows that come after the
ones it was built on, continuing its recursion rather than restarting
it.

## Usage

``` r
term_continue(term, psi, f_past, s_past, newdata, ...)
```

## Arguments

- term:

  A built structural term.

- psi:

  The term's parameters, named as
  [`term_params`](https://statmodels7.github.io/modelterms7/reference/term_params.md).

- f_past:

  The term's contribution at each observed row.

- s_past:

  The driving quantity at each observed row.

- newdata:

  The rows to continue onto.

- ...:

  Passed to methods.

## Value

A numeric vector of `nrow(newdata)` contributions.

## Details

A structural term's contribution at one observation is not a function of
that observation: it is the state a recursion has reached, so predicting
past the series means carrying the state forward. What makes it possible
without simulation is that the quantity driving the recursion has zero
conditional mean – for a score-driven term the score itself – so beyond
the data the recursion is deterministic.

The base method signals an error rather than returning zero: a term with
state that cannot say what its state does next has nothing to offer a
prediction, and a zero would read as a term with no effect.

## See also

[`term_filter`](https://statmodels7.github.io/modelterms7/reference/term_filter.md)

## Examples

``` r
set.seed(1)
dd <- data.frame(t = 1:20, y = rnorm(20))
term <- term_build(gas(p = 1, q = 1, time = t), dd)
psi <- list(omega = 0.1, alpha1 = 0.3, pacf1 = 0.5)
out <- term_filter(term, eta = rep(0, 20), y = dd$y,
                   score = function(e, i) dd$y[i] - e,
                   curvature = function(e, i) -1, psi = psi)
sc <- dd$y - out$eta
term_continue(term, psi, out$eta, sc, data.frame(t = 21:23))
#> [1] 0.3637729 0.2818864 0.2409432
```
