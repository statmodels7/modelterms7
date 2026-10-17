# modelterms7

Model terms as S7 objects. A term is the recipe that a formula names:
the mapping from a data frame to a design-matrix block, the penalty
attached to the block’s coefficients, and the metadata that a fit uses:
coefficient names, hyperparameters with their links, whether the
penalized objective is differentiable, and how effective degrees of
freedom are counted.

The package is part of the [statmodels7](https://statmodels7.github.io)
toolkit.

## Installation

``` r

# install.packages("pak")
pak::pak("statmodels7/modelterms7")
```

## The terms

An additive term contributes a block of columns, so the predictor is
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
  each other, and its matrix parameter carries the dependence on
  whichever [parameters7](https://statmodels7.github.io/parameters7/)
  chart is given to it; a second bar labels effects that share a
  covariance block with those of other terms;
- [`s()`](https://statmodels7.github.io/modelterms7/reference/s.md) and
  [`te()`](https://statmodels7.github.io/modelterms7/reference/te.md)
  are the penalized smooths of one and of several covariates, built on
  [basis7](https://statmodels7.github.io/basis7/);
- [`nl()`](https://statmodels7.github.io/modelterms7/reference/nl.md)
  carries a contribution that is nonlinear in its own parameters, and
  its block is the Jacobian, which is refreshed as the parameters move;
- [`seg()`](https://statmodels7.github.io/modelterms7/reference/seg.md),
  [`jump()`](https://statmodels7.github.io/modelterms7/reference/jump.md)
  and
  [`jseg()`](https://statmodels7.github.io/modelterms7/reference/jseg.md)
  estimate the break-points at which an effect changes slope, level, or
  both, and
  [`seg_start()`](https://statmodels7.github.io/modelterms7/reference/seg_start.md)
  chooses where the iteration begins.

A structural term rewrites the likelihood instead of adding to the
predictor:
[`gas()`](https://statmodels7.github.io/modelterms7/reference/gas.md)
gives score-driven dynamics and
[`regime()`](https://statmodels7.github.io/modelterms7/reference/regime.md)
a latent Markov chain. Both return their exact derivative alongside the
state, because the recursion is the only place where that derivative can
be computed.

A structural term has parameters of its own instead of coefficients, and
it declares them. They are estimated on the unconstrained scale of each
one’s own chart, so a persistence stays stationary and a loading stays
positive for every value that the optimizer proposes.

``` r

g <- gas(p = 1, q = 1)
term_params(g)
#> [1] "omega"  "kappa1" "pacf1"
vapply(term_links(g), function(l) l@link_name, character(1))
#>      omega     kappa1      pacf1 
#> "identity"      "log"   "rhobit"
```

[`cens()`](https://statmodels7.github.io/modelterms7/reference/cens.md)
is not a term. It marks a censored response on the left side of a
formula, and
[`interpret_formula()`](https://statmodels7.github.io/modelterms7/reference/interpret_formula.md)
returns it as the response.
[`statmodels7::statmod()`](https://statmodels7.github.io/statmodels7/reference/statmod.html)
and `rstatmod()` reject a censored response, because a censored
likelihood for it is not yet implemented.

## A term’s own parameters take formulas

A term that has parameters of its own
([`nl()`](https://statmodels7.github.io/modelterms7/reference/nl.md),
[`seg()`](https://statmodels7.github.io/modelterms7/reference/seg.md),
[`jump()`](https://statmodels7.github.io/modelterms7/reference/jump.md),
[`jseg()`](https://statmodels7.github.io/modelterms7/reference/jseg.md)
and
[`gas()`](https://statmodels7.github.io/modelterms7/reference/gas.md))
can develop each of them over covariates by a two-sided formula whose
left side names the parameter. The development acts on the unconstrained
scale of that parameter’s chart, so the constraint holds at every
observation and not merely on average. Examples are a break-point with a
random effect for each group and a loading that depends on a covariate.

``` r

pan <- data.frame(x = runif(60), id = factor(rep(1:6, each = 10)))
pan$y <- rnorm(60)

sg <- term_build(seg(x, psi ~ random(~ 1 | id)), pan)
term_coef_names(sg)[1:4]
#> [1] "seg.beta"             "seg.gamma1"           "seg.psi1.(Intercept)"
#> [4] "seg.psi1.random.1"
```

This is the random-changepoint model, and it needs no separate term. The
same route carries a penalty:
`nl(~ a * exp(-r * x), a ~ 0 + ridge(~ g))` puts a ridge on the
amplitudes across the levels of `g`, and
[`term_penalties()`](https://statmodels7.github.io/modelterms7/reference/term_penalties.md)
reports every entry that a term declares, with the subset of its
parameters that each entry covers.

## Effects that share a covariance block

A second bar labels a random effect, following brms: in
`~ 1 + x | u | id` the last position is the grouping variable, and the
middle one is a label shared by the terms whose effects are correlated
with each other. R nests bars to the left, so the parser yields that
reading. The label is a name and is never evaluated, so a column of the
data with the same name is not read.

``` r

tm <- random(~ 1 + x | u | id)
term_tag(tm)
#> [1] "u"
term_penalties(term_build(tm, pan))
#> list()
```

The empty list is the intended result. The coefficients of a labelled
term share a block with those of every other term that carries the
label, and that block may span another equation, so the fitting layer
builds the prior over the whole block. `statmodels7` collects the terms
with the same label, so
`y ~ random(~ 1 | u | id) | sigma ~ random(~ 1 | u | id)` is a random
intercept in the mean correlated with a random intercept in the scale.

## Break-points three ways

A break-point term estimates the positions in one of three ways. The
default uses the working fits of Muggeo for
[`seg()`](https://statmodels7.github.io/modelterms7/reference/seg.md)
and of Fasola, Muggeo and Kuchenhoff for
[`jump()`](https://statmodels7.github.io/modelterms7/reference/jump.md)
and
[`jseg()`](https://statmodels7.github.io/modelterms7/reference/jseg.md).

`smoothed =` replaces the absolute value inside the term with a
[numericals7](https://statmodels7.github.io/numericals7/)
`abs_smoother`, and each break-point then becomes an ordinary parameter
of a differentiable model.
[`term_jacobian_block()`](https://statmodels7.github.io/modelterms7/reference/term_jacobian_block.md)
shows the difference. The block of
[`seg()`](https://statmodels7.github.io/modelterms7/reference/seg.md) is
already a Jacobian, because the construction is continuous. The two
discontinuous constructions read their position off a product of the
fitted coefficients, so their block is a working linearization with a
frozen weight until the term is smoothed.

``` r

jb <- function(...) term_jacobian_block(term_build(jump(x, ...), pan))
c(sharp    = jb(npsi = 1),
  smoothed = jb(npsi = 1, smoothed = numericals7::smooth_probit()))
#>    sharp smoothed 
#>    FALSE     TRUE
```

The same difference determines which developments a discontinuous term
accepts. On a sharp
[`jump()`](https://statmodels7.github.io/modelterms7/reference/jump.md)
or
[`jseg()`](https://statmodels7.github.io/modelterms7/reference/jseg.md)
the estimated coefficients are products of the change of level and the
development, so a penalty on the development would act on the product,
and the term rejects a development that carries a penalty. A development
without a penalty is accepted, and a smoothed term accepts both kinds,
because its block is a Jacobian.

`marginal = TRUE`, together with a subformula `psi ~ random(~ 1 | g)`,
integrates the latent break-point of each group out of the likelihood.
The prior’s parameters are then estimated by maximum likelihood, and the
positions are not parameters of the fit. For
[`jump()`](https://statmodels7.github.io/modelterms7/reference/jump.md)
the sum over the positions is exact: the process of active break-points
is a hidden Markov chain over the side patterns, so for a group of $`n`$
observations and $`K`$ break-points the forward recursion costs
$`nK2^K`$ where the sum over the cells costs $`(n+1)^K`$. The continuous
kinds,
[`seg()`](https://statmodels7.github.io/modelterms7/reference/seg.md)
and
[`jseg()`](https://statmodels7.github.io/modelterms7/reference/jseg.md),
integrate one break-point per group by Gauss-Kronrod quadrature on each
interval between the group’s ordered observations.

## Simulating from a term

[`term_simulate()`](https://statmodels7.github.io/modelterms7/reference/term_simulate.md)
is the function through which
[`statmodels7::rstatmod()`](https://statmodels7.github.io/statmodels7/reference/rstatmod.html)
draws a response from a term that carries state. A score-driven filter
does not need a second recursion for this. The filter evaluates its
score once per observation, in time order, at the predictor that it has
just produced, so a callback that draws the response at that predictor
and returns the score of the drawn value turns the filter into a
generator. A latent chain draws its path from the stationary
distribution that the likelihood uses for the initial state, since any
other start would simulate a different model. A term without state has
nothing of its own to draw, and the base method signals an error.

## A formula, interpreted

The interpreter evaluates each call in the formula and treats it as a
term when the value inherits from `model_term`. A term class defined
outside the package therefore works in a formula without registration.

``` r

dd <- data.frame(y = rnorm(50), x1 = rnorm(50), x2 = rnorm(50),
                 g = factor(sample(letters[1:5], 50, TRUE)))
dd$L <- matrix(rnorm(50 * 3), 50, 3)

out <- interpret_formula(y ~ x1 + log(abs(x2)) + lasso(L) +
                           random(~ 1 | g), dd)
names(out$terms)
#> [1] "linpar"         "lasso(L)"       "random(~1 | g)"
```

`x1` and `log(abs(x2))` collapse into one `linpar` block: a bare
covariate is a covariate of that block, and so is a call whose value is
not a term. A call whose evaluation fails signals an error that names
the call.

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
is computed from the kink set of the term’s penalty and is not a
property of the term. The lasso block therefore reports `FALSE`, and a
fitting layer can use the answer to select the terms that need a
non-smooth method.

## What a built term carries

``` r

sm <- term_build(s(x1, basis7::bspline_smooth(k = 8)), dd)
dim(term_matrix(sm))
#> [1] 50  7
term_coef_names(sm)
#> [1] "s(x1).lin" "s(x1).z1"  "s(x1).z2"  "s(x1).z3"  "s(x1).z4"  "s(x1).z5" 
#> [7] "s(x1).z6"
term_penalty(sm)@params
#> [1] "lambda"
```

Every built term also carries a blueprint that reproduces the mapping on
new data. Prediction reapplies the factor levels, contrasts and knot
placement recorded at build time and does not derive them again from the
new rows.

``` r

newdd <- dd[1:5, ]
identical(term_predict(sm, newdd), term_matrix(sm)[1:5, ])
#> [1] TRUE
```

[`check_term()`](https://statmodels7.github.io/modelterms7/reference/check_term.md)
validates an additive term against a data frame: the shape and names of
the block, that prediction on the same data reproduces the block, and
that prediction on a subset equals the corresponding rows. A term that
rebuilds its encoding on new rows instead of reapplying its blueprint
fails the last check.

``` r

res <- check_term(linpar(~ x1 + g), dd, verbose = FALSE)
all(res$status == "OK")
#> [1] TRUE
```

[`edf()`](https://statmodels7.github.io/modelterms7/reference/edf.md)
computes effective degrees of freedom from the quantities that a fit
supplies. An unpenalized block counts its columns, a smooth penalty
counts $`\operatorname{tr}[(H + S)^{-1} H]`$, with $`H`$ the unpenalized
curvature and $`S`$ the Hessian of the penalty, and a kinked penalty
(the lasso, the elastic net, SCAD and MCP) counts the nonzero
coefficients.
