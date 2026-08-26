# Whether a Break-Point Term's Block Is a Jacobian

`TRUE` for
[`seg()`](https://statmodels7.github.io/modelterms7/reference/seg.md)
and for any smoothed construction, `FALSE` for the sharp
[`jump()`](https://statmodels7.github.io/modelterms7/reference/jump.md)
and
[`jseg()`](https://statmodels7.github.io/modelterms7/reference/jseg.md).
It is the predicate a fitting layer routes on: a Jacobian licenses a
Gauss-Newton step and a line search on the model's own objective, and a
working linearization does not.

## Arguments

- term:

  A
  [`SegTerm()`](https://statmodels7.github.io/modelterms7/reference/SegTerm.md),
  built or not: the answer is a property of the construction.

- ...:

  Unused.

## Value

A single logical.

## Details

[`seg()`](https://statmodels7.github.io/modelterms7/reference/seg.md) is
differentiable in its break-points away from them, so its block is the
exact derivative.
[`jump()`](https://statmodels7.github.io/modelterms7/reference/jump.md)
and
[`jseg()`](https://statmodels7.github.io/modelterms7/reference/jseg.md)
freeze the weight \\W = 1/(2\lvert\tilde{x}-\psi\rvert)\\ at the
previous position and read the break-point off two coefficients, so
their block is a working linearization: the fixed-point iteration on it
is not a descent method on the model's objective, its early steps going
uphill on purpose under a large scaling factor, and forcing a sufficient
decrease stalls it.

`smoothed` changes the answer for all three. Replacing the step and the
hinge by an
[`penalties7::abs_smoother()`](https://statmodels7.github.io/penalties7/reference/abs_smoother.html)'s
smooth versions makes every break-point an ordinary parameter with a
true derivative, so a smoothed `jump` answers `TRUE` and is routed like
an [`nl()`](https://statmodels7.github.io/modelterms7/reference/nl.md)
term.

## See also

[`term_jacobian_block()`](https://statmodels7.github.io/modelterms7/reference/term_jacobian_block.md)
for the generic and what the answer decides,
[`term_converged()`](https://statmodels7.github.io/modelterms7/reference/term_converged.md)
for the verdict a working linearization gets instead of a score.

## Examples

``` r
# The continuous construction differentiates; the two discontinuous
# ones report a working linearization.
vapply(list(seg = seg(x), jump = jump(x), jseg = jseg(x)),
       term_jacobian_block, logical(1))
#>   seg  jump  jseg 
#>  TRUE FALSE FALSE 

# Smoothing the step makes every break-point an ordinary parameter.
vapply(list(seg = seg(x, smoothed = penalties7::smooth_probit()),
            jump = jump(x, smoothed = penalties7::smooth_probit()),
            jseg = jseg(x, smoothed = penalties7::smooth_probit())),
       term_jacobian_block, logical(1))
#>  seg jump jseg 
#> TRUE TRUE TRUE 
```
