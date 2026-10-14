# The Observed Hessian of a Likelihood Mixed Over States

The log-likelihood of a term that mixes over latent states, with its
exact gradient and its exact Hessian in the whole of a caller's unknown
vector: the coefficients of every equation together with the term's own
parameters.

## Usage

``` r
term_hessian(
  term,
  eta,
  y,
  logdens,
  grad,
  hess,
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

  A built
  [`RegimeTerm()`](https://statmodels7.github.io/modelterms7/reference/RegimeTerm.md).

- eta:

  The static part of the predictor the regimes shift.

- y:

  The response.

- logdens:

  A function of a predictor value and a row index returning the
  log-density there, one value per observation, as
  [`term_loglik()`](https://statmodels7.github.io/modelterms7/reference/term_loglik.md)
  takes it.

- grad:

  A function of the same two arguments returning a matrix with one row
  per observation and one column per distribution parameter, the
  derivative of the log-density in each predictor.

- hess:

  A function of the same two arguments returning an array of dimension
  `c(n, p, p)`, observation first, holding the second derivatives of the
  log-density in the predictors, where `p` is the number of distribution
  parameters.

- psi:

  The term's parameters on the parameter scale, named as
  [`term_params()`](https://statmodels7.github.io/modelterms7/reference/term_params.md).
  The derivatives in the columns `cols` are taken in the unconstrained
  coordinates, so the derivative in a gap is taken in its logarithm, as
  it is for the scale of
  [`term_links()`](https://statmodels7.github.io/modelterms7/reference/term_links.md).

- seed:

  A list with one matrix per distribution parameter, each with one row
  per observation and one column per unknown, giving the derivative of
  that predictor in the caller's unknowns with zeros in the columns the
  term's own parameters occupy.

- cols:

  The columns of the unknown vector the term's own parameters occupy, in
  the order of
  [`term_params()`](https://statmodels7.github.io/modelterms7/reference/term_params.md).

- level:

  The index, among the distribution parameters, of the one the regimes
  shift.

- weights:

  Optional observation weights: one number or one number per
  observation, none missing. `NULL` gives every observation the weight
  one.

- ...:

  Passed to methods.

## Value

A list of three: `loglik`, the `n` per-observation contributions,
without the weights; `gradient`, the derivative of their weighted sum in
the caller's whole unknown vector; and `hessian`, the observed Hessian
of the weighted sum in that same vector, an `m` by `m` matrix with `m`
its length. The Hessian is the mixture's own. Its negative, the observed
information, is smaller in the Loewner order than the complete-data
information, the difference being the information that the unobserved
states cost.

## Details

The smoothed probabilities of
[`term_posterior()`](https://statmodels7.github.io/modelterms7/reference/term_posterior.md)
give the gradient by Fisher's identity. A caller can also assemble from
them the complete-data information, which is the ordinary information
averaged over the states and is the matrix that an EM step inverts. The
observed information of the mixture is a different matrix, smaller than
the complete-data information by the information that the unobserved
states cost, so a standard error computed from the complete-data matrix
is too small.

Louis's identity expresses the difference as the conditional variance of
the complete-data score, which needs the pairwise smoothed probabilities
and a recursion carrying a second moment. This method uses a shorter
route. The scaled forward recursion computes the observed log-likelihood
exactly, as \\\log L = \sum_t \log c_t\\ with \\c_t\\ the normalizing
constant of one step, so differentiating that arithmetic twice gives the
observed Hessian with no identity involved. Carrying \\a_t(k)\\ together
with its first and second derivatives in the unknown vector \\u\\,

\$\$\frac{\partial \log c_t}{\partial u} = \frac{\dot c_t}{c_t}, \qquad
\frac{\partial^2 \log c_t}{\partial u \partial u^\top} = \frac{\ddot
c_t}{c_t} - \frac{\dot c_t}{c_t}\frac{\dot c_t^\top}{c_t},\$\$

and the state is renormalized by the quotient rule,

\$\$\ddot a = \bigl(\ddot{\tilde a} - \dot a \otimes \dot c - \dot c
\otimes \dot a - a\\\ddot c\bigr)/c.\$\$

The emission enters through \\\ddot w = w(gg^\top + G)\\, with \\g\\ and
\\G\\ the score and the Hessian of the log-density in \\u\\ at that
state, which is where the caller's derivatives are used. Groups are
independent series and their contributions add.

The cost is \\O(nK^2m^2)\\ in time. The second derivatives of the state
take \\O(Km^2)\\ in storage, beside \\O(nKm)\\ for the emission scores.
The computation is meant to be run once, at a fitted point.

## See also

[`term_posterior()`](https://statmodels7.github.io/modelterms7/reference/term_posterior.md),
[`term_loglik()`](https://statmodels7.github.io/modelterms7/reference/term_loglik.md)

## Examples

``` r
set.seed(1)
dd <- data.frame(t = 1:30, y = c(rnorm(15), rnorm(15, 3)))
term <- term_build(regime(2, time = t), dd)
# a gaussian mean of unit variance: one predictor, one unknown besides
# the term's own four
m <- 5L
out <- term_hessian(
  term, eta = rep(0, 30), y = dd$y,
  logdens = function(e, i) dnorm(dd$y[i], e, log = TRUE),
  grad = function(e, i) matrix(dd$y[i] - e, ncol = 1L),
  hess = function(e, i) array(-1, c(length(i), 1L, 1L)),
  psi = list(level1 = 0, gap2 = 3, alr1.1 = 2, alr2.1 = -2),
  seed = list(matrix(c(rep(1, 30), rep(0, 30 * 4)), 30, m)),
  cols = 2:5, level = 1L)
dim(out$hessian)
#> [1] 5 5
```
