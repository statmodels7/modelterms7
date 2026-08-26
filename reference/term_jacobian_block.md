# Is a Term's Block the Jacobian of Its Contribution?

`TRUE` when the block a term reports is the exact derivative of its
contribution in its own coefficients, `FALSE` when it is a working
linearization with quantities frozen at the previous iterate. The base
method returns `TRUE`, which a fixed design satisfies trivially and
which
[`nl()`](https://statmodels7.github.io/modelterms7/reference/nl.md) and
[`seg()`](https://statmodels7.github.io/modelterms7/reference/seg.md)
satisfy by construction;
[`jump()`](https://statmodels7.github.io/modelterms7/reference/jump.md)
and
[`jseg()`](https://statmodels7.github.io/modelterms7/reference/jseg.md)
answer `FALSE`.

## Usage

``` r
term_jacobian_block(term, ...)
```

## Arguments

- term:

  A term, built or not: the answer is a property of the construction.

- ...:

  Passed to methods. No shipped method reads anything here.

## Value

A single logical.

## What the answer decides

Where the block is a Jacobian, a scoring step on it is a Gauss-Newton
step, a line search on the model's own objective is licensed, and the
term is fitted inside the same system as everything else.

Where it is a working linearization the fixed-point iteration of Fasola,
Muggeo and Kuchenhoff (2018) is not a descent method on the model's
objective. Its early steps go uphill on purpose, under a scaling factor
that anneals, so forcing a sufficient decrease on it stalls the
iteration. Such a term is fitted by alternating exact working fits at
the frozen block with
[`term_refresh()`](https://statmodels7.github.io/modelterms7/reference/term_refresh.md),
and its convergence is what
[`term_converged()`](https://statmodels7.github.io/modelterms7/reference/term_converged.md)
answers instead of a score.

The two discontinuous constructions are the ones that answer `FALSE`.
Their weight \\W = 1/(2\lvert\tilde x - \psi\rvert)\\ is held at the
previous break-point, and the position is read off two coefficients
rather than being one.

## Smoothing changes the answer

`jump(x, smoothed = ...)` replaces the step by a smooth surrogate, and
the break-point becomes an ordinary parameter with a true Jacobian, so a
smoothed break-point term answers `TRUE`. That is how a fitting layer
routes it without a special case.

## References

Fasola, S., Muggeo, V. M. R. and Kuchenhoff, H. (2018). A heuristic,
iterative algorithm for change-point detection in abrupt change models.
*Computational Statistics*, 33, 997–1015.

## See also

[`term_refresh()`](https://statmodels7.github.io/modelterms7/reference/term_refresh.md)
for the block being recomputed,
[`term_converged()`](https://statmodels7.github.io/modelterms7/reference/term_converged.md)
for the verdict on such a term,
[`seg()`](https://statmodels7.github.io/modelterms7/reference/seg.md),
[`jump()`](https://statmodels7.github.io/modelterms7/reference/jump.md)
and
[`jseg()`](https://statmodels7.github.io/modelterms7/reference/jseg.md).

## Examples

``` r
# The continuous construction differentiates; the two discontinuous
# ones report a working linearization.
vapply(list(seg = seg(x), jump = jump(x), jseg = jseg(x)),
       term_jacobian_block, logical(1))
#>   seg  jump  jseg 
#>  TRUE FALSE FALSE 

# A fixed design is a Jacobian trivially, and so is nl().
term_jacobian_block(linpar(~ x))
#> [1] TRUE
term_jacobian_block(nl(~ a * x, start = list(a = 1)))
#> [1] TRUE

# Smoothing the step makes the break-point an ordinary parameter.
term_jacobian_block(jump(x, smoothed = penalties7::smooth_probit()))
#> [1] TRUE
```
