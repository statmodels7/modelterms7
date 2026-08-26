# Log-Likelihood Contributions of a Structural Term

The per-observation log-likelihood contributions a structural term
produces, with their derivatives in the term's own parameters. This is
the second shape the structural branch takes, beside
[`term_filter()`](https://statmodels7.github.io/modelterms7/reference/term_filter.md):
a term that shifts the predictor implements the filter, and one that
rewrites the likelihood itself implements this, its contribution being
no predictor at all.

## Usage

``` r
term_loglik(term, eta, y, logdens, score, psi, ...)
```

## Arguments

- term:

  A built structural term.

- eta:

  The static part of the linear predictor, one value per observation.

- y:

  The response. It reaches the recursion through the callbacks; the
  argument is passed for methods that need it directly.

- logdens:

  A function `(e, i)` returning the log-density of observation `i` at
  predictor `e`.

- score:

  A function `(e, i)` returning the derivative of that log-density with
  respect to the predictor.

- psi:

  The term's parameters on the **parameter** scale, named as
  [`term_params()`](https://statmodels7.github.io/modelterms7/reference/term_params.md).

- ...:

  Passed to methods.

## Value

A list of two: `loglik`, a numeric vector of `n` contributions summing
to the term's log-likelihood, and `jacobian`, an `n` by `length(psi)`
matrix of their derivatives in the term's parameters, exact and
propagated through the recursion.

## What the contributions are

The contribution of observation \\t\\ is the logarithm of its one-step
predictive density given everything before it,

\$\$\ell_t(\psi) = \log f(y_t \mid y_1, \dots, y\_{t-1}; \psi), \qquad
\sum\_{t=1}^{n} \ell_t(\psi) = \log f(y_1, \dots, y_n; \psi),\$\$

so the vector returned sums to the term's log-likelihood by the chain
rule of probability, whatever the dependence between observations.

For
[`regime()`](https://statmodels7.github.io/modelterms7/reference/regime.md)
it is the normalizing constant of the forward recursion, \\\ell_t = \log
\sum\_{k} \pi\_{t \mid t-1, k} f(y_t \mid S_t = k)\\, and the Jacobian
\\\partial \ell_t / \partial \psi_j\\ is propagated beside the filtered
distribution, never differenced.

## How the model reaches it

The term knows the chain and the levels and nothing about the family.
The two callbacks are how the model's own density enters: `logdens`
returns the log-density of one observation at a given predictor and
`score` its derivative in that predictor. Both are called with the
predictor shifted by each regime's level, so a term of \\K\\ regimes
evaluates them \\K\\ times per observation.

Unlike
[`term_filter()`](https://statmodels7.github.io/modelterms7/reference/term_filter.md)'s
callbacks these do **not** depend on the state the recursion has
reached: a regime shifts a predictor known in advance. That is why the
density and the score of every observation under every regime can be
computed once, vectorized, before the recursion starts.

The method on
[`structural_term()`](https://statmodels7.github.io/modelterms7/reference/structural_term.md)
throws, naming the class: a term of the filter shape implements
[`term_filter()`](https://statmodels7.github.io/modelterms7/reference/term_filter.md)
instead, and a fitting layer tells the two apart by which one answers.

## See also

[`regime()`](https://statmodels7.github.io/modelterms7/reference/regime.md)
and the marginal break-point terms for the two implementations,
[`term_filter()`](https://statmodels7.github.io/modelterms7/reference/term_filter.md)
for the other structural shape,
[`term_posterior()`](https://statmodels7.github.io/modelterms7/reference/term_posterior.md)
for the smoothed states,
[`term_hessian()`](https://statmodels7.github.io/modelterms7/reference/term_hessian.md)
for the observed information.

## Examples

``` r
set.seed(1)
dd <- data.frame(t = 1:40, y = c(rnorm(20), rnorm(20, 3)))
term <- term_build(regime(2, time = t), dd)
psi <- list(level1 = 0, gap2 = 3, alr1.1 = 2, alr2.1 = -2)

out <- term_loglik(term, rep(0, 40), dd$y,
                   logdens = function(e, i) dnorm(dd$y[i], e, log = TRUE),
                   score = function(e, i) dd$y[i] - e,
                   psi = psi)
sum(out$loglik)
#> [1] -59.65148
dim(out$jacobian)
#> [1] 40  4

# The forward recursion is the sum over every state path. On eight
# observations that sum can be taken in full, and it agrees exactly.
d2 <- data.frame(t = 1:8, y = c(rnorm(4), rnorm(4, 3)))
t2 <- term_build(regime(2, time = t), d2)
fwd <- sum(term_loglik(t2, rep(0, 8), d2$y,
                       logdens = function(e, i) dnorm(d2$y[i], e, log = TRUE),
                       score = function(e, i) d2$y[i] - e,
                       psi = psi)$loglik)

P <- as.matrix(parameters7::param_value(
  parameters7::transition_matrix(2), c(alr1.1 = 2, alr2.1 = -2)))
p0 <- Re(eigen(t(P))$vectors[, 1]); p0 <- p0 / sum(p0)
lev <- c(0, 3)                       # psi is on the parameter scale
paths <- as.matrix(expand.grid(rep(list(1:2), 8)))
tot <- sum(apply(paths, 1, function(s) {
  pr <- p0[s[1]]
  for (j in 2:8) pr <- pr * P[s[j - 1], s[j]]
  pr * prod(dnorm(d2$y, lev[s]))
}))
c(forward = fwd, all_256_paths = log(tot))
#>       forward all_256_paths 
#>     -12.07278     -12.07278 
```
