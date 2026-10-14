# Hold a Break-Point Term at Its Positions

Puts a sharp
[`jump()`](https://statmodels7.github.io/modelterms7/reference/jump.md)
or
[`jseg()`](https://statmodels7.github.io/modelterms7/reference/jseg.md)
term in its held state, or releases it from that state. A held term
keeps its break-points at their current positions, and its block is the
exact design at those positions: the covariate for the linear part,
where the term has one, the truncated line \\(x - \psi_k)\_+\\ for each
change of slope, the indicator \\1(x \> \psi_k)\\ for each change of
level, and a column of zeros in the slot that carried the position.

## Usage

``` r
seg_hold(term, hold = TRUE)
```

## Arguments

- term:

  A built
  [`jump()`](https://statmodels7.github.io/modelterms7/reference/jump.md)
  or
  [`jseg()`](https://statmodels7.github.io/modelterms7/reference/jseg.md)
  term.

- hold:

  `TRUE` to hold the term, `FALSE` to release it. Releasing a term that
  is not held returns it unchanged.

## Value

The term in the requested state, its block rebuilt.

## Details

The working construction of Fasola, Muggeo and Kuchenhoff (2018) reaches
its break-points through a block whose weight is frozen at the previous
iterate. At the fixed point of that iteration the coefficients are
therefore those of the working model, and not the least-squares
coefficients at the positions reached. With the positions held, the term
contributes a linear function of its remaining coefficients, and a fit
of those coefficients is an ordinary fit, which is how `segmented` and
`stepmented` report their coefficients. For a held term,
[`term_held()`](https://statmodels7.github.io/modelterms7/reference/term_held.md)
returns the slots that are held,
[`term_jacobian_block()`](https://statmodels7.github.io/modelterms7/reference/term_jacobian_block.md)
returns `TRUE`,
[`term_refresh()`](https://statmodels7.github.io/modelterms7/reference/term_refresh.md)
does not move the positions,
[`term_converged()`](https://statmodels7.github.io/modelterms7/reference/term_converged.md)
returns `TRUE` and
[`term_stalled()`](https://statmodels7.github.io/modelterms7/reference/term_stalled.md)
returns `FALSE`.

A term whose break-points carry a development is held in the same way.
The positions are then the values \\w_i'p_k\\ of the development at the
current coefficients \\p_k\\.

Releasing the term writes the position back into the slot as \\g_k =
-\delta_k\psi_k\\, with the coefficient vector \\p_k\\ of the
development in place of \\\psi_k\\ for a developed position. The
read-off of the working construction then returns the held positions,
and the iteration can continue from them.

A continuous
[`seg()`](https://statmodels7.github.io/modelterms7/reference/seg.md)
and a smoothed term are rejected, since their positions are ordinary
coefficients, each with a column of its own.

## References

Fasola, S., Muggeo, V. M. R. and Kuchenhoff, H. (2018). A heuristic,
iterative algorithm for change-point detection in abrupt change models.
*Computational Statistics*, 33, 997–1015.

## See also

[`term_held()`](https://statmodels7.github.io/modelterms7/reference/term_held.md),
[`seg_polish_exact()`](https://statmodels7.github.io/modelterms7/reference/seg_polish_exact.md),
[`seg_relocate()`](https://statmodels7.github.io/modelterms7/reference/seg_relocate.md).

## Examples

``` r
set.seed(1)
d <- data.frame(x = sort(runif(100, 0, 10)))
d$y <- 1 + 2 * (d$x > 6) + rnorm(100, sd = 0.3)
b <- seg_hold(term_build(jump(x, psi = 6), d))
term_matrix(b)[c(1, 100), ]
#>      jump.delta1 jump.g1
#> [1,]           0       0
#> [2,]           1       0
seg_psi(b)
#> [1] 6
```
