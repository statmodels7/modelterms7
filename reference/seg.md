# Segmented, Stepmented and Segmented-with-Jump Terms

A covariate whose effect changes at break-points estimated with
everything else. `seg` changes slope at each break-point and stays
continuous (muggeo2003); `jump` steps to a new level and is
discontinuous (fasola2018); `jseg` does both at the same points.

## Usage

``` r
seg(
  x,
  npsi = 1,
  psi = NULL,
  by = NULL,
  linear = TRUE,
  penalty = c("none", "lasso", "ridge"),
  band = 0.02,
  label = "seg"
)

jump(
  x,
  npsi = 1,
  psi = NULL,
  by = NULL,
  linear = TRUE,
  penalty = c("none", "lasso", "ridge"),
  band = 0.02,
  label = "jump"
)

jseg(
  x,
  npsi = 1,
  psi = NULL,
  by = NULL,
  linear = TRUE,
  penalty = c("none", "lasso", "ridge"),
  band = 0.02,
  label = "jseg"
)
```

## Arguments

- x:

  The covariate, an expression evaluated in the data.

- npsi:

  The number of break-points. Defaults to 1.

- psi:

  Optional starting positions; defaults to evenly spaced quantiles of
  the covariate.

- by:

  An optional factor giving an independent set of break-points per
  level.

- linear:

  Whether the block carries the linear effect. Defaults to `TRUE`.

- penalty:

  One of `"none"` (default), `"lasso"` or `"ridge"`, applied to the
  changes.

- band:

  For a discontinuous term, the half-width of the band around a
  break-point over which the step is replaced by a ramp, as a fraction
  of the covariate's range. Defaults to `0.02`; see Details.

- label:

  A single non-empty string prefixed to the coefficient names.

## Value

An object of class
[`SegTerm`](https://statmodels7.github.io/modelterms7/reference/SegTerm.md)
(a specification; see
[`term_build`](https://statmodels7.github.io/modelterms7/reference/term_build.md)).

## Details

Both constructions rest on the same device: a quantity that depends on
the break-point non-linearly is written as a linear function of it once
something is frozen at the previous iterate, so that a linear fit
returns the new break-point. Iterating that is the estimation algorithm,
and here it is the
[`term_refresh`](https://statmodels7.github.io/modelterms7/reference/term_refresh.md)
contract of
[`nl`](https://statmodels7.github.io/modelterms7/reference/nl.md) with a
different block: refreshing the term at the current coefficients and
fitting is one step of it.

### The continuous case

With \\f(x) = \delta (x-\psi)\_+\\ the contribution is differentiable in
\\\psi\\ away from the break-point, and \\\partial f/\partial\psi =
-\delta\\\mathbb{1}(x\>\psi)\\, so the design block is the ordinary
Jacobian and the break-point is an ordinary coefficient, updated by the
increment a Gauss-Newton step solves for. That is the algorithm of
muggeo2003 written in the coordinates the rest of this package uses: his
working variables \\U\\ and \\V\\ are its columns, and his update \\\psi
\leftarrow \psi + \gamma/\delta\\ is that step.

### The discontinuous case

A jump is not differentiable in \\\psi\\ at all, and yet needs no search
over candidate positions. The identity \$\$\mathbb{1}(x\>\psi) =
\frac{1}{2} + \frac{x-\psi}{2\lvert x-\psi \rvert}\$\$ holds exactly for
\\x \neq \psi\\, and is linear in \\\psi\\ once the weight \\1/(2\lvert
x - \psi\rvert)\\ is held at the previous iterate. Writing \\W =
1/(2\lvert x-\psi^{0}\rvert)\\ and \\Z = xW + 1/2\\, a jump of size
\\\kappa\\ at \\\psi\\ is \\\kappa Z + gW\\ with \\g = -\kappa\psi\\, so
a linear fit on \\(Z, W)\\ returns the break-point as \\\psi =
-g/\kappa\\ (fasola2018). The break-point is therefore not a coefficient
here but a quantity read off two of them, which is why refreshing the
term recovers it before rebuilding the weights.

The weight is unbounded as \\x\\ approaches \\\psi\\, and that matters
more than it looks. Since \\Z - \psi W\\ is the step itself, \\Z\\ is
\\\psi W\\ plus a quantity of order one: let \\W\\ grow without bound
and the two columns become numerically collinear, drowning the very
signal the fit is meant to read. The denominator is therefore held at or
above `band` times the covariate's range, which caps \\W\\.

That is a bandwidth and not a guard. Within the band the step is
replaced by a ramp, so the fixed point of the iteration is that of a
slightly smoothed problem; a narrower band is more faithful and worse
conditioned. The segmented literature makes the same trade by displacing
the observations nearest a break-point instead of capping the weight.

### What the term carries

With `linear = TRUE` the block carries the linear effect too, so the
term is the whole relationship rather than the change in it. `by` gives
an independent set of break-points and changes per level of a factor.
`penalty` puts a penalty on the changes themselves – the slope changes
for `seg`, the jump sizes for `jump` – through a map that selects those
coefficients and leaves the linear effect and the break-points alone;
with `"lasso"` that is a selection of how many break-points are really
there.

A break-point is confined to the interval between the 5th and the 95th
percentile of the covariate. Outside it the block is singular rather
than merely ill-conditioned: with \\\psi\\ below the smallest
observation the indicator is constant, so the truncated line and that
constant are linearly dependent, and with the linear effect present so
are all three columns. A run that ends against the limit has not located
a break-point, and
[`seg_psi`](https://statmodels7.github.io/modelterms7/reference/seg_psi.md)
then returns the limit itself.

The objective has local optima in the break-points, and the iteration
converges from within a basin around the true position rather than from
anywhere. On a joint jump and change of slope at \\x = 5\\ in 500
observations, swept over eight samples and damping factors from 0.05 to
1, every start at 4 or above recovers the break-point at every damping
below 1, every start at 2 or below fails at all of them, and a start at
3 succeeds for some samples and not others. The step also has to be
damped for its own sake: taken whole it overshoots even from a good
start. A run should therefore be started from several positions, which
is what
[`multistart`](https://statmodels7.github.io/optimizers7/reference/multistart.html)
does and what the bootstrap restarting of the segmented literature is
for.

## References

Muggeo, V. M. R. (2003). Estimating regression models with unknown
break-points. *Statistics in Medicine*, 22(19), 3055–3071.

Fasola, S., Muggeo, V. M. R. and Kuchenhoff, H. (2018). A heuristic,
iterative algorithm for change-point detection in abrupt change models.
*Computational Statistics*, 33, 997–1015.

## Examples

``` r
set.seed(1)
dd <- data.frame(x = sort(runif(200, 0, 10)))
dd$y <- 1 + 0.5 * dd$x + 2 * pmax(dd$x - 6, 0) + rnorm(200, sd = 0.3)

built <- term_build(seg(x), dd)
term_coef_names(built)
#> [1] "seg.lin"    "seg.delta1" "seg.psi1"  
seg_psi(built)
#> [1] 5.054907
```
