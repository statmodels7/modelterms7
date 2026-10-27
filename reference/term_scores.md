# Per-Observation Scores of a Structural Term

Returns, for each observation, the derivative of the weighted increment
of the term's likelihood in the caller's whole unknown vector. The sum
of the rows is the `gradient` that
[`term_hessian()`](https://statmodels7.github.io/modelterms7/reference/term_hessian.md)
returns.

## Usage

``` r
term_scores(
  term,
  eta,
  y,
  logdens,
  grad,
  psi,
  seed,
  cols,
  level,
  weights = NULL,
  ...
)
```

## Arguments

- term:

  A built structural term.

- eta:

  The static predictor.

- y:

  The response.

- logdens, grad:

  The log-density and its first derivatives in the predictors, as in
  [`term_hessian()`](https://statmodels7.github.io/modelterms7/reference/term_hessian.md).

- psi:

  The term's parameters on the parameter scale.

- seed:

  The derivative of each predictor in the caller's unknowns.

- cols:

  The columns that the term's own parameters occupy.

- level:

  The distribution parameter that the regimes shift.

- weights:

  Optional observation weights: one number or one number per
  observation, without missing values. `NULL` gives every observation
  the weight one.

- ...:

  Passed to methods.

## Value

A matrix with one row per observation and one column per unknown of the
caller, in the order of the rows of the data. Each row is multiplied by
the weight of its observation.

## Details

For a
[`regime()`](https://statmodels7.github.io/modelterms7/reference/regime.md)
term the increment of observation \\t\\ is \\w_t \log p(y_t \mid y_1,
\dots, y\_{t-1})\\, with \\p\\ the conditional density given by the
scaled forward recursion, and the row is its derivative. The derivative
is carried in forward mode through the recursion, which is the
first-order part of the recursion of
[`term_hessian()`](https://statmodels7.github.io/modelterms7/reference/term_hessian.md).
A term class that does not implement the generic signals an error.

## See also

[`term_hessian()`](https://statmodels7.github.io/modelterms7/reference/term_hessian.md)

## Examples

``` r
set.seed(1)
dd <- data.frame(t = 1:30, y = c(rnorm(15), rnorm(15, 3)))
term <- term_build(regime(2, time = t), dd)
psi <- list(level1 = 0, gap2 = 3, alr1.1 = 2, alr2.1 = -2)
seed <- list(matrix(c(rep(1, 30), rep(0, 30 * 4)), 30, 5L))
ld <- function(e, i) dnorm(dd$y[i], e, log = TRUE)
gr <- function(e, i) matrix(dd$y[i] - e, ncol = 1L)
U <- term_scores(term, eta = rep(0, 30), y = dd$y, logdens = ld,
                 grad = gr, psi = psi, seed = seed, cols = 2:5,
                 level = 1L)
dim(U)
#> [1] 30  5
```
