# Segmented Term: a Broken Line with Estimated Break-Points

`seg()` specifies a covariate whose slope changes at \\K\\ break-points
estimated with everything else, and the relationship stays continuous
(Muggeo, 2003). Written in the equation of any parameter of any
distribution, the term contributes

\$\$\beta x_i + \sum\_{k=1}^{K} \gamma_k\\(x_i - \psi_k)\_{+}, \qquad
(u)\_{+} = \max(0, u),\$\$

to that parameter's linear predictor, so that in a model \\g(\theta_i) =
\eta_i\\ in which the rest of the equation supplies \\z_i'\alpha\\,

\$\$\eta_i = z_i'\alpha + \beta x_i + \sum\_{k=1}^{K} \gamma_k\\(x_i -
\psi_k)\_{+}.\$\$

\\\beta\\ is the slope before the first break-point and \\\gamma_k\\ the
change of slope at \\\psi_k\\, so the slope over the \\k\\-th segment is
\\\beta + \gamma_1 + \cdots + \gamma_k\\. The gaussian additive model is
the case \\g = \mathrm{id}\\ on the mean of a `gaussian1_distrib`, and
the term applies to any other distribution and parameter in the same
way.

## Usage

``` r
seg(
  x,
  ...,
  npsi = 1,
  psi = NULL,
  by = NULL,
  linear = TRUE,
  smoothed = NULL,
  marginal = FALSE,
  n_boot = 10,
  label = "seg"
)
```

## Arguments

- x:

  The covariate, an expression evaluated in the data.

- ...:

  Two-sided formulas developing the term's own coefficients; see the
  section above. Cannot be combined with `by`.

- npsi:

  The number of break-points. Defaults to 1.

- psi:

  Optional starting positions; defaults to evenly spaced quantiles of
  the covariate. Where a break-point carries a development, they seed
  it, and each starting vector solves \\Wp_k \approx \psi_k^{0}\\.

- by:

  A one-sided formula giving every coefficient of the term the same
  development, e.g. `by = ~0 + g` for an independent set per level of a
  factor. A bare variable is rejected; write the formula.

- linear:

  Whether the block carries the linear effect \\\beta x\\. Defaults to
  `TRUE`.

- smoothed:

  `NULL` (the default: the construction exactly as documented above) or
  a numericals7
  [`numericals7::abs_smoother()`](https://statmodels7.github.io/numericals7/reference/abs_smoother.html),
  e.g.
  [`numericals7::smooth_probit()`](https://statmodels7.github.io/numericals7/reference/smooth_probit.html).
  The smoother replaces the step and the hinge by their smooth versions,
  \\(1 + s'(u))/2\\ and \\(u + s(u))/2\\, so every break-point becomes
  an ordinary parameter of a model as smooth as the smoother
  (\\C^\infty\\ for the probit and the hyperbolic, \\C^3\\ for the
  quintic; the page of the quintic gives its cost under an outer
  criterion, where the probit is the recommended smoother). There is no
  working parametrization, no auxiliary coefficient and no scaling
  schedule (for
  [`jump()`](https://statmodels7.github.io/modelterms7/reference/jump.md)
  and
  [`jseg()`](https://statmodels7.github.io/modelterms7/reference/jseg.md),
  a `c0` that is passed is ignored with a message). The block is the
  true Jacobian and the term is fitted by Gauss-Newton like
  [`nl()`](https://statmodels7.github.io/modelterms7/reference/nl.md). A
  smoothed term accepts a development of any coefficient, including a
  break-point (`psi ~ random(~1 | id)`) and penalized developments, for
  every kind, since a smoothed term has no read-off that would constrain
  the development. The smoother's width is resolved at build from the
  covariate's spacing (the median gap between distinct values, within
  groups where a break-point development supplies a partition) unless
  the object carries one. For the quintic, which is exact outside its
  width, the width is raised to at least 0.55 times the largest gap over
  the range that a break-point may take, so that a break-point cannot
  sit in a gap with no observation within its width. The width is
  reported: it is the width of the transition (the bent-cable reading).
  The smoothing bias is zero outside the transition for the quintic,
  decays with the Gaussian tails for the probit, and decays as
  \\c/(4\|u\|)\\ for the hyperbolic, whose parameter \\c\\ is a squared
  length. The objective remains multimodal in the positions: smoothing
  rounds the local optima and does not remove them, so the profile start
  and the `n_boot` restarts remain necessary.

- marginal:

  Whether the break-point is a latent variable per group, integrated out
  of the likelihood exactly instead of being an estimated position.
  `FALSE`, the default, is the construction documented above. `TRUE`
  requires the subformula `psi ~ random(~1 | g)` and returns a
  structural term of the likelihood shape
  ([`MarginalBreakTerm()`](https://statmodels7.github.io/modelterms7/reference/MarginalBreakTerm.md));
  see the section on the marginal construction in
  [`jump()`](https://statmodels7.github.io/modelterms7/reference/jump.md).
  For a `seg` term the native random-changepoint fit
  (`psi ~ random(~1 | g)` without `marginal`) estimates the break-points
  as the mode of the latent distribution and is the recommended route,
  while the marginal construction is the exact-likelihood alternative.
  The marginal `seg` has one break-point: the conditional is smooth in
  the position, so the integral runs on Gauss-Kronrod panels over the
  intervals between a group's ordered observations, and several latent
  break-points would need a product quadrature with a number of
  components that a fitting layer cannot carry.

- n_boot:

  How many bootstrap restarts the fitting layer runs after the iteration
  first converges (Wood 2001, the device that `segmented` runs by
  default): each restart re-estimates on a bootstrap resample from the
  current break-points and then on the data again from where the
  resample ended, keeping the better fit. The objective has local optima
  in the break-points, and with several of them a grid start does not
  reach the right basin from every placement. Defaults to 10,
  `segmented`'s own default; 0 disables. The term itself only declares
  the value; running the restarts is the fitting layer's, as with a
  penalty's hyperparameters.

- label:

  A single non-empty string prefixed to the coefficient names.

## Value

An object of class
[`SegTerm()`](https://statmodels7.github.io/modelterms7/reference/SegTerm.md)
(a specification; see
[`term_build()`](https://statmodels7.github.io/modelterms7/reference/term_build.md)),
or of class
[`MarginalBreakTerm()`](https://statmodels7.github.io/modelterms7/reference/MarginalBreakTerm.md)
with `marginal = TRUE`.

## Details

### How the break-points are estimated

The contribution is differentiable in \\\psi_k\\ away from the
break-point, with \\\partial/\partial\psi_k =
-\gamma_k\\\mathbb{1}(x\>\psi_k)\\, so the design block is the ordinary
Jacobian and the break-point is an ordinary coefficient. The working
variables \\U = (x-\psi)\_{+}\\ and \\V = -\mathbb{1}(x\>\psi)\\ of
Muggeo (2003) are, up to the factor \\\gamma_k\\ on \\V\\, the columns
of the block for \\\gamma_k\\ and \\\psi_k\\, so a Gauss-Newton step on
the block is an iteration of his linearization. An iteration of this
scheme refreshes the term at the current coefficients
([`term_refresh()`](https://statmodels7.github.io/modelterms7/reference/term_refresh.md))
and then fits, and
[`nl()`](https://statmodels7.github.io/modelterms7/reference/nl.md) uses
the same scheme.

The objective has local optima in the break-points and the iteration
converges from within a basin around where it starts, so
[`seg_start()`](https://statmodels7.github.io/modelterms7/reference/seg_start.md)
chooses the starting positions on a grid and not at a conventional
point. A break-point is confined to the interval between the 5th and the
95th percentile of the covariate, so that both sides of a break-point
carry observations. A position that reaches a limit is clamped there,
and the column of a clamped break-point is zero, since the contribution
no longer moves with it. A run ending against the limit has not located
a break-point, and
[`seg_psi()`](https://statmodels7.github.io/modelterms7/reference/seg_psi.md)
then reports the limit itself.

The iteration can settle into a cycle of period two in the break-point,
in which case the rule of
[`seg_converged()`](https://statmodels7.github.io/modelterms7/reference/seg_step.md)
is never met while the objective has long since stopped moving. A caller
that can evaluate the objective should stop on its relative change, as
`segmented` does.

## The coefficients

The coefficients are `beta` (present when `linear = TRUE`), `gamma1` ...
`gammaK` and `psi1` ... `psiK`, in that order. The break-points are
coefficients of the working block here, and they are not in the sharp
discontinuous constructions, so
[`seg_psi()`](https://statmodels7.github.io/modelterms7/reference/seg_psi.md)
reports them under every construction.

## Developing a coefficient with covariates

Every coefficient of the term is a parameter like any other and may be
developed on covariates, through a two-sided formula in `...` whose left
side names it:

    seg(x, psi ~ id)                   a break-point per subject
    seg(x, psi ~ random(~1 | id))      the random-changepoint model
    seg(x, gamma1 ~ z)                 a first slope change varying with z
    seg(x, by = ~0 + g)                every coefficient, per level of g

The right side goes through
[`interpret_formula()`](https://statmodels7.github.io/modelterms7/reference/interpret_formula.md)
and takes any term with a fixed design block, penalized ones included. A
structural term
([`gas()`](https://statmodels7.github.io/modelterms7/reference/gas.md),
[`regime()`](https://statmodels7.github.io/modelterms7/reference/regime.md))
and a term whose block moves with its coefficients (`seg()`,
[`nl()`](https://statmodels7.github.io/modelterms7/reference/nl.md)) are
rejected. The hyperparameters that a sub-term declares are reported by
[`term_penalties()`](https://statmodels7.github.io/modelterms7/reference/term_penalties.md)
and estimated by the fitting layer. The left side names either one
coefficient (`gamma2`) or a whole kind (`gamma`, meaning every
\\\gamma_k\\, and each of them has coefficients of its own over the
shared design). `by = ~f` is the shorthand giving every coefficient the
same development, and does not combine with per-coefficient formulas.

A developed coefficient carries one value per observation:
\\\gamma\_{k,i} = w_i'a_k\\ and \\\psi\_{k,i} = w_i'p_k\\, with \\w_i\\
the row of the sub-design. Both enter the block linearly: the truncated
line becomes \\(x_i-\psi\_{k,i})\_{+}w_i'\\ and the Jacobian column
becomes \\-\gamma\_{k,i}\\\mathbb{1}(x_i\>\psi\_{k,i})\\w_i'\\, so the
continuous construction accepts a development of any coefficient and any
combination of them.

A penalty on the changes themselves, the lasso that selects the number
of break-points, is the development on an intercept alone,
`seg(x, npsi = 4, gamma ~ 0 + lasso(~1))`: the \\K\\ changes are then
one penalized block under one hyperparameter, because the entries of a
subformula shared by every coefficient of a kind are pooled into one. A
subformula follows the intercept convention of a formula, so
`gamma ~ lasso(~1)` carries the subformula's own unpenalized intercept
beside the penalized column, which for an intercept-only development
repeats the same constant column. That design is rejected when the term
is built, and the `0 +` removes the subformula's own intercept.

## The linear effect

With `linear = TRUE`, the default, the term carries \\\beta x\\ itself,
so `y ~ seg(x)` is the whole relationship. Writing the covariate again
in the same equation is then exactly collinear with that column, and
[`interpret_formula()`](https://statmodels7.github.io/modelterms7/reference/interpret_formula.md)
removes it with a warning. `linear = FALSE` is the explicit way to keep
the linear effect outside the term, and `y ~ x + seg(x, linear = FALSE)`
is the same model written the other way round.

## References

Muggeo, V. M. R. (2003). Estimating regression models with unknown
break-points. *Statistics in Medicine*, 22(19), 3055–3071.

Muggeo, V. M. R., Atkins, D. C., Gallop, R. J. and Dimidjian, S. (2014).
Segmented mixed models with random changepoints: a maximum likelihood
approach with application to treatment for depression study.
*Statistical Modelling*, 14(4), 293–313.

Wood, S. N. (2001). Minimizing model fitting objectives that contain
spurious local minima by bootstrap restarting. *Biometrics*, 57(1),
240–244.

## See also

[`jump()`](https://statmodels7.github.io/modelterms7/reference/jump.md),
[`jseg()`](https://statmodels7.github.io/modelterms7/reference/jseg.md),
[`seg_psi()`](https://statmodels7.github.io/modelterms7/reference/seg_psi.md),
[`seg_start()`](https://statmodels7.github.io/modelterms7/reference/seg_start.md),
[`seg_step()`](https://statmodels7.github.io/modelterms7/reference/seg_step.md),
[`nl()`](https://statmodels7.github.io/modelterms7/reference/nl.md)

## Examples

``` r
set.seed(1)
dd <- data.frame(x = sort(runif(200, 0, 10)))
dd$y <- 1 + 0.5 * dd$x + 2 * pmax(dd$x - 6, 0) + rnorm(200, sd = 0.3)

built <- term_build(seg(x), dd)
term_coef_names(built)
#> [1] "seg.beta"   "seg.gamma1" "seg.psi1"  
seg_psi(built)
#> [1] 5.054907


# The term is fitted next. The data are simulated from a known truth, so the
# estimates below can be read against it.
if (requireNamespace("statmodels7", quietly = TRUE)) {
  set.seed(1)
  fd <- data.frame(x = sort(runif(200, 0, 10)))
  fd$y <- 1 + 0.5 * fd$x + 2 * pmax(fd$x - 6, 0) + rnorm(200, sd = 0.3)
  ft <- statmodels7::statmod(y ~ seg(x),
                             distributions7::gaussian1_distrib(), fd)
  # The truth is intercept 1, slope 0.5, change of slope 2 and break-point 6.
  round(coef(ft)$mu, 2)
}
#> (Intercept)    seg.beta  seg.gamma1    seg.psi1 
#>        0.96        0.52        2.04        6.08 
```
