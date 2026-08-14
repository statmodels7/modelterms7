# Segmented-with-Jump Term: Slope and Level Both Changing

A covariate whose slope *and* level both change at \\K\\ break-points
estimated with everything else, the union of
[`seg`](https://statmodels7.github.io/modelterms7/reference/seg.md) and
[`jump`](https://statmodels7.github.io/modelterms7/reference/jump.md) at
the same points. Written in the equation of any parameter of any
distribution, the term contributes

\$\$\beta x_i + \sum\_{k=1}^{K} \bigl\[\delta_k\\\mathbb{1}(x_i \geq
\psi_k) + \gamma_k\\(x_i - \psi_k)\_{+}\bigr\]\$\$

to that parameter's linear predictor, so that in a model \\g(\theta_i) =
\eta_i\\ with the rest of the equation supplying \\z_i'\alpha\\,

\$\$\eta_i = z_i'\alpha + \beta x_i + \sum\_{k=1}^{K}
\bigl\[\delta_k\\\mathbb{1}(x_i \geq \psi_k) + \gamma_k\\(x_i -
\psi_k)\_{+}\bigr\].\$\$

\\\beta\\ is the slope before the first break-point, \\\gamma_k\\ the
change of slope at \\\psi_k\\ and \\\delta_k\\ the change of level
there.

## Usage

``` r
jseg(
  x,
  ...,
  npsi = 1,
  psi = NULL,
  by = NULL,
  linear = TRUE,
  c0 = 0.05,
  label = "jseg"
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

- c0:

  The starting value of the scaling factor that separates the
  observations from the break-point, as a fraction of the distance to
  the ends of the range. Defaults to `0.05`, the value fasola2018
  recommend; see Details.

- label:

  A single non-empty string prefixed to the coefficient names.

## Value

An object of class
[`SegTerm`](https://statmodels7.github.io/modelterms7/reference/SegTerm.md)
(a specification; see
[`term_build`](https://statmodels7.github.io/modelterms7/reference/term_build.md)).

## Details

Both devices of the two constructions are used at once: the truncated
line is differentiable in the break-point and the step is linearized by
the identity of fasola2018, on the rescaled covariate its scaling
schedule provides (see
[`jump`](https://statmodels7.github.io/modelterms7/reference/jump.md)).
What is particular to the joint term is the read-off, since the
truncated line depends on the break-point as well and reading only the
step's pair would discard what the change of slope says about it.
Linearizing both parts about the previous position, and writing \\h\\
for the increment, \\U\\ for the truncated line and \\I = Z -
\psi^{0}W\\ for the indicator,

\$\$\gamma\\(x-\psi)\_{+} \approx \gamma U - \gamma h I, \qquad
\delta\\\mathbb{1}(x\>\psi) = \delta Z - \delta\psi W,\$\$

so the fitted coefficients of \\Z\\ and \\W\\ are \\a = \delta - \gamma
h\\ and \\b = \gamma h\psi^{0} - \delta\psi^{0} - \delta h\\, and the
increment solves

\$\$\gamma h^{2} + a h + (b + a\psi^{0}) = 0,\$\$

the root of smaller modulus. At \\\gamma = 0\\ the quadratic degenerates
to \\h = -(b + a\psi^{0})/a\\, that is \\\psi = -b/a\\, so the pure step
is the case this contains rather than an exception to it.

Everything
[`jump`](https://statmodels7.github.io/modelterms7/reference/jump.md)
documents about the scaling schedule, the confinement of a break-point
and the choice of starting positions applies here unchanged.

## The coefficients

`beta` (present when `linear = TRUE`), `gamma1` ... `gammaK`, the
changes of slope, then `delta1` ... `deltaK`, the changes of level, then
`g1` ... `gK`, the auxiliary coefficients from which the break-points
are read. Ask
[`seg_psi`](https://statmodels7.github.io/modelterms7/reference/seg_psi.md)
for the break-points.

## Developing a coefficient with covariates

`beta` and `gammaK` take a development on any design, both entering the
block linearly. `deltaK` and `psiK` enter the quadratic above, which is
a scalar identity: it splits observation by observation only where the
design has one column per group and each observation belongs to one, so
`jseg(x, by = ~0 + g)` is an independent set of everything per level and
a development of `delta` or `psi` on any other design is rejected. The
componentwise reading that would remain diverges whenever a change of
level passes near zero mid-iteration, which on a joint model it
routinely does; it is rejected rather than shipped diverging.

## References

Muggeo, V. M. R. (2003). Estimating regression models with unknown
break-points. *Statistics in Medicine*, 22(19), 3055–3071.

Fasola, S., Muggeo, V. M. R. and Kuchenhoff, H. (2018). A heuristic,
iterative algorithm for change-point detection in abrupt change models.
*Computational Statistics*, 33, 997–1015.

## See also

[`seg`](https://statmodels7.github.io/modelterms7/reference/seg.md),
[`jump`](https://statmodels7.github.io/modelterms7/reference/jump.md),
[`seg_psi`](https://statmodels7.github.io/modelterms7/reference/seg_psi.md),
[`seg_start`](https://statmodels7.github.io/modelterms7/reference/seg_start.md),
[`seg_step`](https://statmodels7.github.io/modelterms7/reference/seg_step.md)

## Examples

``` r
set.seed(1)
dd <- data.frame(x = sort(runif(200, 0, 10)))
dd$y <- 1 + 0.3 * dd$x + 2 * (dd$x > 6) +
  1.5 * pmax(dd$x - 6, 0) + rnorm(200, sd = 0.3)

built <- term_build(jseg(x), dd)
term_coef_names(built)
#> [1] "jseg.beta"   "jseg.gamma1" "jseg.delta1" "jseg.g1"    
```
