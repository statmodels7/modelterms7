# Recompute a Break-Point Term at New Coefficients

Moves the break-points to where the coefficients put them, rebuilds the
working block there, and, for a sharp `jump` or `jseg`, advances the
scaling schedule. A smoothed term has no schedule, and a term held by
[`seg_hold()`](https://statmodels7.github.io/modelterms7/reference/seg_hold.md)
keeps its positions and only rebuilds the exact block at the new
coefficients. It is one step of the iteration that estimates the
break-points, and the block that it returns is what the next linear fit
is taken on.

## Arguments

- term:

  A built
  [`SegTerm()`](https://statmodels7.github.io/modelterms7/reference/SegTerm.md).

- coef:

  The coefficients to move to, of length `ncol(term@X)`. Any other
  length signals an error that names the required length.

- ...:

  Unused.

## Value

The term with its block, its break-points, its contribution and its
scaling schedule all recomputed at `coef`.

## How the new positions are found

For `seg`, and for any smoothed term, the break-point is a coefficient,
so the new position is read straight out of `coef`. For a sharp `jump`
and `jseg` it is **read off** two of them: the identity \\1(x \> \psi) =
1/2 + (x-\psi)/(2\|x-\psi\|)\\ makes the step linear in \\\psi\\ once
the weight is frozen, and \\\psi = -g_k/\delta_k\\ follows. `jseg` reads
a quadratic instead, its truncated line depending on the position as
well.

Every new position is clamped into the confinement interval that
[`term_build()`](https://statmodels7.github.io/modelterms7/reference/term_build.md)
set. When the term has more than one break-point and none of its
per-break-point coefficients carries a development, crossed break-points
are relabeled so that each keeps traveling with its own scaling factor
and direction.

## The scaling schedule

The rescaling of Fasola, Muggeo and Kuchenhoff (2018) opens a gap of
relative width \\c\\ around each break-point, and the factor is **halved
whenever that break-point reverses direction**, which is the signal that
the iteration has begun to circle an optimum instead of traveling toward
one. The factor is both the conditioning device and the step control, so
it must advance once per committed step: refreshing inside a line search
would anneal at every trial point, and refreshing from the specification
would freeze it at `c0`.
[`seg_reheat()`](https://statmodels7.github.io/modelterms7/reference/seg_reheat.md)
resets it.

## See also

[`term_refresh()`](https://statmodels7.github.io/modelterms7/reference/term_refresh.md)
for the generic,
[`seg_psi()`](https://statmodels7.github.io/modelterms7/reference/seg_psi.md)
for the positions after a refresh,
[`seg_step()`](https://statmodels7.github.io/modelterms7/reference/seg_step.md)
and
[`term_converged()`](https://statmodels7.github.io/modelterms7/reference/term_converged.md)
for the verdict,
[`seg_reheat()`](https://statmodels7.github.io/modelterms7/reference/seg_reheat.md)
to start the schedule again.

## Examples

``` r
set.seed(1)
d <- data.frame(x = sort(runif(120, 0, 10)))
d$y <- 1 + 0.5 * d$x + 2 * pmax(d$x - 6, 0) + rnorm(120, sd = 0.4)
b <- term_build(seg(x, npsi = 1), d)

# Moving the break-point coefficient moves the block.
cf <- b@blueprint$coef
r <- term_refresh(b, c(cf[1], cf[2], 6))
c(before = seg_psi(b), after = seg_psi(r))
#>   before    after 
#> 4.803127 6.000000 
max(abs(term_matrix(r) - term_matrix(b))) > 0
#> [1] TRUE

# A coefficient vector of the wrong length is an error.
try(term_refresh(b, c(1, 2)))
#> Error : 'coef' must have length 3.
```
