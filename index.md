# modelterms7

Model terms as S7 objects. A term is the recipe a formula names: the
mapping from a data frame to a design-matrix block, the penalty attached
to the block’s coefficients, and the metadata a fit reads – coefficient
names, hyperparameters with their links, whether the penalized objective
is differentiable, and how effective degrees of freedom are counted.

Part of the [statmodels7](https://statmodels7.github.io) toolkit.

## Installation

``` r

# install.packages("pak")
pak::pak("statmodels7/modelterms7")
```

## The terms

An **additive** term contributes a block of columns, so the predictor is
$`\eta = \sum_t X_t \beta_t`$ and a penalized fit minimizes
$`-\ell(\beta) + \sum_t \rho_t(\beta_t; \theta_t)`$:

- [`linpar()`](https://statmodels7.github.io/modelterms7/reference/linpar.md)
  is the unpenalized parametric block;
- [`ridge()`](https://statmodels7.github.io/modelterms7/reference/ridge.md),
  [`lasso()`](https://statmodels7.github.io/modelterms7/reference/lasso.md),
  [`scad()`](https://statmodels7.github.io/modelterms7/reference/scad.md),
  [`mcp()`](https://statmodels7.github.io/modelterms7/reference/mcp.md)
  and
  [`enet()`](https://statmodels7.github.io/modelterms7/reference/enet.md)
  attach the corresponding
  [penalties7](https://statmodels7.github.io/penalties7/) object to
  their block;
- [`random()`](https://statmodels7.github.io/modelterms7/reference/random.md)
  builds grouped intercepts and slopes with the effect distribution as
  the penalty: a multivariate one lets the effects of a group depend on
  each other, its matrix parameter carrying the dependence on whichever
  [parameters7](https://statmodels7.github.io/parameters7/) chart is
  given it, and a second bar labels effects that share a covariance
  block with those of other terms;
- [`s()`](https://statmodels7.github.io/modelterms7/reference/s.md) and
  [`te()`](https://statmodels7.github.io/modelterms7/reference/te.md)
  are the penalized smooths of one and of several covariates, built on
  [basis7](https://statmodels7.github.io/basis7/);
- [`nl()`](https://statmodels7.github.io/modelterms7/reference/nl.md)
  carries a contribution nonlinear in its own parameters, its block the
  Jacobian refreshed as they move;
- [`seg()`](https://statmodels7.github.io/modelterms7/reference/seg.md),
  [`jump()`](https://statmodels7.github.io/modelterms7/reference/jump.md)
  and
  [`jseg()`](https://statmodels7.github.io/modelterms7/reference/jseg.md)
  estimate the break-points at which an effect changes slope, level, or
  both, with
  [`seg_start()`](https://statmodels7.github.io/modelterms7/reference/seg_start.md)
  choosing where the iteration begins.

A **structural** term rewrites the likelihood rather than adding to it:
[`gas()`](https://statmodels7.github.io/modelterms7/reference/gas.md)
gives score-driven dynamics and
[`regime()`](https://statmodels7.github.io/modelterms7/reference/regime.md)
a latent Markov chain, both returning their exact derivative alongside
the state because the recursion is the only place it can be computed.
[`cens()`](https://statmodels7.github.io/modelterms7/reference/cens.md)
marks a censored response on the left side of the formula.

A structural term has parameters of its own rather than coefficients,
and declares them: they are estimated on the unconstrained scale of each
one’s own chart, so a persistence stays stationary and a loading stays
positive whatever value the optimizer proposes.

``` r

g <- gas(p = 1, q = 1)
term_params(g)
#> [1] "omega"  "alpha1" "pacf1"
vapply(term_links(g), function(l) l@link_name, character(1))
#>      omega     alpha1      pacf1 
#> "identity"      "log"   "rhobit"
```

## A term’s own parameters take formulas

Anything a term carries can be developed over covariates by writing a
two-sided formula whose left side names the parameter. The development
acts on the unconstrained scale of that parameter’s chart, so the
constraint holds at every observation and not merely on average: a
break-point given a random effect per group, a loading given a
covariate.

``` r

pan <- data.frame(x = runif(60), id = factor(rep(1:6, each = 10)))
pan$y <- rnorm(60)

sg <- term_build(seg(x, psi ~ random(~ 1 | id)), pan)
term_coef_names(sg)[1:4]
#> [1] "seg.beta"             "seg.gamma1"           "seg.psi1.(Intercept)"
#> [4] "seg.psi1.random.1"
```

That is the random-changepoint model, and it is one line rather than a
separate term. The same route carries a penalty:
`nl(~ a * exp(-r * x), a ~ 0 + ridge(~ g))` puts a ridge on the
amplitudes across the levels of `g`, and
[`term_penalties()`](https://statmodels7.github.io/modelterms7/reference/term_penalties.md)
reports every entry a term declares with the subset of its parameters
each one covers.

## Effects that share a covariance block

A second bar labels a random effect, following **brms**: in
`~ 1 + x | u | id` the last position is the grouping variable and the
middle one says which effects are correlated with each other. R nests
bars to the left, so that reading comes out of the parser. The label is
a name and is never evaluated, so a column of the data carrying it is
not looked at.

``` r

tm <- random(~ 1 + x | u | id)
term_tag(tm)
#> [1] "u"
term_penalties(term_build(tm, pan))
#> list()
```

The empty list is the answer rather than a gap: the coefficients of a
labelled term share a block with those of every other term carrying the
label, and that block may span another equation entirely, so the prior
over it belongs to the class and is built where the model is fitted.
`statmodels7` collects them, so
`y ~ random(~ 1 | u | id) | sigma ~ random(~ 1 | u | id)` is a random
intercept in the mean correlated with a random intercept in the scale.

## Break-points three ways

The break-point terms come in three constructions, and which to use is a
statement about the model rather than about the numerics. The default
estimates the positions by the working fits of Muggeo and of Fasola,
with the position read off the fitted coefficients.

`smoothed =` replaces the absolute value inside the term with a
[penalties7](https://statmodels7.github.io/penalties7/) `abs_smoother`,
and then a break-point is an ordinary parameter of a differentiable
model.
[`term_jacobian_block()`](https://statmodels7.github.io/modelterms7/reference/term_jacobian_block.md)
is where the difference shows.
[`seg()`](https://statmodels7.github.io/modelterms7/reference/seg.md) is
continuous and its block is already a Jacobian; the two discontinuous
constructions read their position off a product of the fitted
coefficients, so their block is a working linearization with a frozen
weight until the term is smoothed.

``` r

jb <- function(...) term_jacobian_block(term_build(jump(x, ...), pan))
c(sharp    = jb(npsi = 1),
  smoothed = jb(npsi = 1, smoothed = penalties7::smooth_probit()))
#>    sharp smoothed 
#>    FALSE     TRUE
```

That is also what makes the developments above available on a
discontinuous term: the read-off cannot carry them and a Jacobian can.

`marginal = TRUE` integrates a latent break-point per group out of the
likelihood instead of estimating one, which is the exact treatment
rather than a mode-based one. For a jump the integral is exact: the
process of active break-points is a hidden Markov chain over the side
patterns, so the forward recursion costs $`nK2^K`$ where the cell sum
costs $`(n+1)^K`$.

## Simulating from a term

[`term_simulate()`](https://statmodels7.github.io/modelterms7/reference/term_simulate.md)
is how
[`statmodels7::rstatmod()`](https://statmodels7.github.io/statmodels7/reference/rstatmod.html)
draws a response from a term that carries state. A score-driven filter
needs no second recursion for it: the filter evaluates its score once
per observation, in time order, at the predictor it has just produced,
so a callback that draws the response there and returns the score of
what it drew turns the filter into a generator. A latent chain draws its
path from the stationary law the likelihood is written with, since any
other start would simulate a different model. A term without state has
nothing of its own to draw, and the base method signals an error.

## A formula, interpreted

The interpreter recognizes a term by what its call evaluates to, not by
its name, so a term class defined outside the package works in a formula
without registration.

``` r

dd <- data.frame(y = rnorm(50), x1 = rnorm(50), x2 = rnorm(50),
                 g = factor(sample(letters[1:5], 50, TRUE)))
dd$L <- matrix(rnorm(50 * 3), 50, 3)

out <- interpret_formula(y ~ x1 + log(abs(x2)) + lasso(L) +
                           random(~ 1 | g), dd)
names(out$terms)
#> [1] "linpar"         "lasso(L)"       "random(~1 | g)"
```

`x1` and `log(abs(x2))` collapse into one `linpar`: a bare covariate is
a covariate, and so is a call whose value is not a term.

``` r

built <- lapply(out$terms, term_build, data = dd)
vapply(built, term_npar, integer(1))
#>         linpar       lasso(L) random(~1 | g) 
#>              3              3              5
vapply(built, term_smooth, logical(1))
#>         linpar       lasso(L) random(~1 | g) 
#>           TRUE          FALSE           TRUE
```

[`term_smooth()`](https://statmodels7.github.io/modelterms7/reference/term_smooth.md)
is read from the penalty rather than declared by the term, so the lasso
block reports `FALSE` and the model layer knows which coefficients need
a non-smooth method.

## What a built term carries

``` r

sm <- term_build(s(x1, k = 8), dd)
dim(term_matrix(sm))
#> [1] 50  7
term_coef_names(sm)
#> [1] "s(x1).lin" "s(x1).z1"  "s(x1).z2"  "s(x1).z3"  "s(x1).z4"  "s(x1).z5" 
#> [7] "s(x1).z6"
term_penalty(sm)@params
#> [1] "lambda"
```

Every built term also carries a blueprint that reproduces the mapping on
new data, reapplying the factor levels, contrasts and knot placement
recorded at build time rather than deriving them again.

``` r

newdd <- dd[1:5, ]
identical(term_predict(sm, newdd), term_matrix(sm)[1:5, ])
#> [1] TRUE
```

[`check_term()`](https://statmodels7.github.io/modelterms7/reference/check_term.md)
validates a term against a data frame – the block’s shape and names,
that prediction on the same data reproduces the block, and that
prediction on a subset equals the corresponding rows, which is the check
a term that rebuilds instead of reapplying fails.

``` r

res <- check_term(linpar(~ x1 + g), dd, verbose = FALSE)
all(res$status == "OK")
#> [1] TRUE
```

[`edf()`](https://statmodels7.github.io/modelterms7/reference/edf.md)
computes effective degrees of freedom from the pieces a fit supplies:
exactly the column count for an unpenalized block,
$`\operatorname{tr}[(H + S)^{-1} H]`$ for a smooth penalty, and the
number of nonzero coefficients for the lasso, SCAD and MCP.
