# The Posterior Shift of a Term Mixed Over States, and Its Jacobian

For a structural term whose contribution is a likelihood mixed over
latent states, the shift that the term adds to its equation's predictor
at each observation, averaged over the smoothed posterior of the states,
together with its exact derivative in the whole of a caller's unknown
vector.

## Usage

``` r
term_shift_jacobian(term, eta, y, logdens, grad, psi, seed, cols, level, ...)
```

## Arguments

- term:

  A built term of the likelihood shape.

- eta:

  The static part of the predictor the term shifts.

- y:

  The response.

- logdens:

  A function of a predictor value and a row index returning the
  log-density there, as
  [`term_loglik()`](https://statmodels7.github.io/modelterms7/reference/term_loglik.md)
  takes it.

- grad:

  A function of the same two arguments returning a matrix with one row
  per observation and one column per distribution parameter, the
  derivative of the log-density in each predictor.

- psi:

  The term's parameters on the parameter scale, named as
  [`term_params()`](https://statmodels7.github.io/modelterms7/reference/term_params.md).

- seed:

  A list with one matrix per distribution parameter, each with one row
  per observation and one column per unknown, giving the derivative of
  that predictor in the caller's unknowns, with zeros in the columns
  that the term's own parameters occupy.

- cols:

  The columns of the unknown vector that the term's own parameters
  occupy, in the order of
  [`term_params()`](https://statmodels7.github.io/modelterms7/reference/term_params.md).

- level:

  The index, among the distribution parameters, of the one the term
  shifts.

- ...:

  Passed to methods.

## Value

A list with `shift`, the `n` posterior shifts, and `jacobian`, an `n` by
`m` matrix with `m` the length of the unknown vector.

## Details

A term of the likelihood shape reports no predictor of its own. The
predictor of its equation at observation \\i\\ is the static part plus
the posterior mean of the shift, \$\$\bar s_i = \sum_c \gamma\_{ic}\\
s\_{ic},\$\$ with \\\gamma\_{ic}\\ the posterior probability of
component \\c\\ given the whole series
([`term_posterior()`](https://statmodels7.github.io/modelterms7/reference/term_posterior.md))
and \\s\_{ic}\\ its shift
([`term_levels()`](https://statmodels7.github.io/modelterms7/reference/term_levels.md)).
The derivative in the unknown vector \\u\\ is \$\$\frac{\partial \bar
s_i}{\partial u} = \sum_c \frac{\partial \gamma\_{ic}}{\partial u}\\
s\_{ic} + \gamma\_{ic}\\ \frac{\partial s\_{ic}}{\partial u},\$\$ and
the posterior moves with every unknown: with the coefficients of every
equation, through the emissions, and with the term's own parameters,
through the emissions, the shifts and the prior of the states. The delta
method for a prediction reads this derivative.

Where the posterior is the smoothed marginal of a hidden Markov chain (a
regime chain, and the side chain of a marginal step), the derivatives
are carried by forward mode through the normalized forward and backward
recursions and combined by the quotient rule; the scale of each step
multiplies a whole vector and cancels in the normalization of the
posterior, so its own derivative is not needed. Where the posterior is a
group's distribution over quadrature nodes (the continuous marginal
kinds), \\\partial \gamma_c = \gamma_c (\partial a_c - \sum\_{c'}
\gamma\_{c'} \partial a\_{c'})\\ with \\a_c\\ the log of the node's
weight times the group's emissions.

The derivatives in the term's own parameters are taken on the
unconstrained scale of
[`term_links()`](https://statmodels7.github.io/modelterms7/reference/term_links.md),
as
[`term_hessian()`](https://statmodels7.github.io/modelterms7/reference/term_hessian.md)
takes them.

## See also

[`term_posterior()`](https://statmodels7.github.io/modelterms7/reference/term_posterior.md),
[`term_levels()`](https://statmodels7.github.io/modelterms7/reference/term_levels.md),
[`term_hessian()`](https://statmodels7.github.io/modelterms7/reference/term_hessian.md)

## Examples

``` r
set.seed(1)
dd <- data.frame(t = 1:30, y = c(rnorm(15), rnorm(15, 3)))
term <- term_build(regime(2, time = t), dd)
psi <- list(level1 = 0, gap2 = 3, alr1.1 = 2, alr2.1 = -2)
out <- term_shift_jacobian(
  term, eta = rep(0, 30), y = dd$y,
  logdens = function(e, i) dnorm(dd$y[i], e, log = TRUE),
  grad = function(e, i) matrix(dd$y[i] - e, ncol = 1L),
  psi = psi, seed = list(matrix(c(rep(1, 30), rep(0, 30 * 4)), 30, 5)),
  cols = 2:5, level = 1L)
dim(out$jacobian)
#> [1] 30  5
```
