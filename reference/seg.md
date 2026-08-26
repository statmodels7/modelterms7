# Segmented Term: a Broken Line with Estimated Break-Points

A covariate whose slope changes at \\K\\ break-points estimated with
everything else, the relationship staying continuous (muggeo2003).
Written in the equation of any parameter of any distribution, the term
contributes

\$\$\beta x_i + \sum\_{k=1}^{K} \gamma_k\\(x_i - \psi_k)\_{+}, \qquad
(u)\_{+} = \max(0, u),\$\$

to that parameter's linear predictor, so that in a model \\g(\theta_i) =
\eta_i\\ with the rest of the equation supplying \\z_i'\alpha\\,

\$\$\eta_i = z_i'\alpha + \beta x_i + \sum\_{k=1}^{K} \gamma_k\\(x_i -
\psi_k)\_{+}.\$\$

\\\beta\\ is the slope before the first break-point and \\\gamma_k\\ the
change of slope at \\\psi_k\\, so the slope over the \\k\\-th segment is
\\\beta + \gamma_1 + \cdots + \gamma_k\\. The gaussian additive model is
the case \\g = \mathrm{id}\\ on the mean of a `gaussian1_distrib`;
nothing here is particular to it.

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
  the covariate. Where a break-point carries a development they seed it,
  each starting vector solving \\Wp_k \approx \psi_k^{0}\\.

- by:

  A one-sided formula giving every coefficient of the term the same
  development, e.g. `by = ~0 + g` for an independent set per level of a
  factor. A bare variable is rejected; write the formula.

- linear:

  Whether the block carries the linear effect \\\beta x\\. Defaults to
  `TRUE`.

- smoothed:

  `NULL` (the default: the construction exactly as documented above) or
  a penalties7
  [`penalties7::abs_smoother()`](https://statmodels7.github.io/penalties7/reference/abs_smoother.html),
  e.g.
  [`penalties7::smooth_probit()`](https://statmodels7.github.io/penalties7/reference/smooth_probit.html).
  The smoother replaces the step and the hinge by their smooth versions,
  \\(1 + s'(u))/2\\ and \\(u + s(u))/2\\, so every break-point becomes
  an ordinary parameter of a \\C^\infty\\ model: there is no working
  parametrization, no auxiliary coefficient and no scaling schedule
  (`c0` is ignored, with a message), the block is the true Jacobian and
  the term is fitted by Gauss-Newton like
  [`nl()`](https://statmodels7.github.io/modelterms7/reference/nl.md). A
  development of a break-point, `psi ~ random(~1 | id)` and penalized
  ones included, is then legal for every kind, the read-off that
  constrained the discontinuous constructions having gone. The
  smoother's width is resolved at build from the covariate's spacing
  (the median gap between distinct values, within groups where a
  break-point development supplies a partition) unless the object
  carries one, and is reported: it is the width of the transition, the
  bent-cable reading, and the smoothing bias it buys is confined to a
  window of that width (probit, quintic) or decays as \\c/(4\|u\|)\\
  (hyperbolic). The objective is still multimodal in the positions –
  smoothing rounds the local optima and does not remove them, so the
  profile start and the `n_boot` restarts stay necessary; a smoothed fit
  from a bad start has been measured converging to an absurd local
  optimum while reporting success.

- marginal:

  Whether the break-point is a latent variable per group, integrated out
  of the likelihood exactly instead of being an estimated position.
  `FALSE`, the default, is the construction documented above. `TRUE`
  requires the subformula `psi ~ random(~1 | g)` and returns a
  structural term of the likelihood shape
  ([`MarginalBreakTerm()`](https://statmodels7.github.io/modelterms7/reference/MarginalBreakTerm.md));
  see the section of
  [`jump()`](https://statmodels7.github.io/modelterms7/reference/jump.md),
  whose step model is where the marginal buys the most: for a `seg` term
  the native random-changepoint fit (`psi ~ random(~1 | g)` without
  `marginal`) is measured equivalent and remains the recommended route,
  the marginal being the exact-likelihood alternative. One break-point
  here: the conditional is smooth in the position, so the integral runs
  on a Gauss-Kronrod panel per interval between a group's ordered
  observations, and several latents would need a product quadrature
  whose component count no fitting layer can carry.

- n_boot:

  How many bootstrap restarts the fitting layer runs after the iteration
  first converges (Wood 2001, the device `segmented` runs by default):
  each restart re-estimates on a bootstrap resample from the current
  break-points and then on the data again from where the resample ended,
  keeping the better fit. The objective has local optima in the
  break-points, and with several of them a grid start does not reach the
  right basin from every placement. Defaults to 10, `segmented`'s own
  default; 0 disables. The term itself only declares the value; running
  the restarts is the fitting layer's, as with a penalty's
  hyperparameters.

- label:

  A single non-empty string prefixed to the coefficient names.

## Value

An object of class
[`SegTerm()`](https://statmodels7.github.io/modelterms7/reference/SegTerm.md)
(a specification; see
[`term_build()`](https://statmodels7.github.io/modelterms7/reference/term_build.md)).

## Details

### How the break-points are estimated

The contribution is differentiable in \\\psi_k\\ away from the
break-point, with \\\partial/\partial\psi_k =
-\gamma_k\\\mathbb{1}(x\>\psi_k)\\, so the design block is the ordinary
Jacobian, the break-point is an ordinary coefficient, and the increment
a Gauss-Newton step solves for is the update of muggeo2003: his working
variables \\U\\ and \\V\\ are the block's columns and his \\\psi
\leftarrow \psi + \gamma/\delta\\ is that step. Refreshing the term at
the current coefficients
([`term_refresh()`](https://statmodels7.github.io/modelterms7/reference/term_refresh.md))
and fitting is one iteration of it, which is the same contract
[`nl()`](https://statmodels7.github.io/modelterms7/reference/nl.md)
uses.

The objective has local optima in the break-points and the iteration
converges from within a basin around where it starts, so
[`seg_start()`](https://statmodels7.github.io/modelterms7/reference/seg_start.md)
chooses the starting positions on a grid rather than at a conventional
point. A break-point is confined to the interval between the 5th and the
95th percentile of the covariate: outside it the indicator is constant,
the truncated line and that constant are linearly dependent, and the
block is exactly singular, not merely ill-conditioned. A run ending
against the limit has not located a break-point, and
[`seg_psi()`](https://statmodels7.github.io/modelterms7/reference/seg_psi.md)
then reports the limit itself.

The iteration can settle into a cycle of period two in the break-point,
in which case the rule of
[`seg_converged()`](https://statmodels7.github.io/modelterms7/reference/seg_step.md)
is never met while the objective has long since stopped moving; a caller
that can evaluate the objective should stop on its relative change, as
`segmented` does.

## The coefficients

`beta` (present when `linear = TRUE`), `gamma1` ... `gammaK` and `psi1`
... `psiK`, in that order. The break-points are coefficients of the
working block here, which they are not in the discontinuous
constructions, so
[`seg_psi()`](https://statmodels7.github.io/modelterms7/reference/seg_psi.md)
is the function that reports them either way.

## Developing a coefficient with covariates

Every coefficient of the term is a parameter like any other and may be
developed on covariates, through a two-sided formula in `...` whose left
side names it:

    seg(x, psi ~ id)                   a break-point per subject
    seg(x, psi ~ random(~1 | id))      the random-changepoint model
    seg(x, gamma1 ~ z)                 a first slope change varying with z
    seg(x, by = ~0 + g)                every coefficient, per level of g

The right side goes through
[`interpret_formula()`](https://statmodels7.github.io/modelterms7/reference/interpret_formula.md),
so it takes *any* term of this package, penalized ones included, and the
hyperparameters a sub-term declares are reported by
[`term_penalties()`](https://statmodels7.github.io/modelterms7/reference/term_penalties.md)
and estimated by the fitting layer. The left side names either one
coefficient (`gamma2`) or a whole kind (`gamma`, meaning every
\\\gamma_k\\, each with coefficients of its own over the shared design).
`by = ~f` is the shorthand giving every coefficient the same
development, and does not combine with per-coefficient formulas.

A developed coefficient carries one value per observation:
\\\gamma\_{k,i} = w_i'a_k\\ and \\\psi\_{k,i} = w_i'p_k\\, with \\w_i\\
the row of the sub-design. Both enter the block linearly, the truncated
line becoming \\(x_i-\psi\_{k,i})\_{+}w_i'\\ and the Jacobian column
\\-\gamma\_{k,i}\\\mathbb{1}(x_i\>\psi\_{k,i})\\w_i'\\, so the
continuous construction accepts a development of any coefficient and any
combination of them.

A penalty on the changes themselves, the lasso that selects how many
break-points are really there, is the development on an intercept alone,
`seg(x, npsi = 4, gamma ~ 0 + lasso(~1))`: the \\K\\ changes are then
one penalized block under one hyperparameter, the entries of a
subformula shared by every coefficient of a kind being pooled into one.
The `0 +` is not decoration. A subformula follows the intercept
convention of a formula, so `gamma ~ lasso(~1)` carries the subformula's
own unpenalized intercept beside the penalized column, which for an
intercept-only development is the same column twice.

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


# Fitted. The data are simulated from a known truth, so the
# estimates below can be read against it.
if (requireNamespace("statmodels7", quietly = TRUE)) {
  set.seed(1)
  fd <- data.frame(x = sort(runif(200, 0, 10)))
  fd$y <- 1 + 0.5 * fd$x + 2 * pmax(fd$x - 6, 0) + rnorm(200, sd = 0.3)
  ft <- statmodels7::statmod(y ~ seg(x),
                             distributions7::gaussian1_distrib(), fd)
  # truth: intercept 1, slope 0.5, change of slope 2, break-point 6
  round(coef(ft)$mu, 2)
}
#> (Intercept)    seg.beta  seg.gamma1    seg.psi1 
#>        0.96        0.52        2.04        6.08 
```
