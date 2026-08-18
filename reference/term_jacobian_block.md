# Is a Term's Block the Jacobian of Its Contribution?

`TRUE` when the design block a term reports is the exact derivative of
its contribution in its own coefficients, `FALSE` when it is a working
linearization with quantities frozen at the previous iterate. The base
method returns `TRUE`, which is what a fixed design satisfies trivially
and what
[`nl`](https://statmodels7.github.io/modelterms7/reference/nl.md) and
[`seg`](https://statmodels7.github.io/modelterms7/reference/seg.md)
satisfy by construction.

## Usage

``` r
term_jacobian_block(term, ...)
```

## Arguments

- term:

  A term (built or not; the answer is a property of the construction).

- ...:

  Passed to methods.

## Value

A single logical.

## Details

The distinction decides how a fitting layer may treat the block. Where
the block is a Jacobian, a scoring step on it is a Gauss–Newton step and
a line search on the model's own objective is licensed, so the term can
be fitted inside the same system as everything else. Where it is a
working linearization –
[`jump`](https://statmodels7.github.io/modelterms7/reference/jump.md)
and
[`jseg`](https://statmodels7.github.io/modelterms7/reference/jseg.md),
whose weight \\W = 1/(2\lvert \tilde x - \psi\rvert)\\ is held at the
previous break-point and whose position is read off two coefficients –
the fixed-point iteration of fasola2018 is not a descent method on the
model's objective, and forcing a sufficient decrease on it stalls the
iteration. Such a term is fitted by alternating exact working fits at
the frozen block with
[`term_refresh`](https://statmodels7.github.io/modelterms7/reference/term_refresh.md),
and its convergence is what
[`term_converged`](https://statmodels7.github.io/modelterms7/reference/term_converged.md)
answers rather than a score.

## References

Fasola, S., Muggeo, V. M. R. and Kuchenhoff, H. (2018). A heuristic,
iterative algorithm for change-point detection in abrupt change models.
*Computational Statistics*, 33, 997–1015.

## See also

[`term_refresh`](https://statmodels7.github.io/modelterms7/reference/term_refresh.md),
[`term_converged`](https://statmodels7.github.io/modelterms7/reference/term_converged.md)

## Examples

``` r
term_jacobian_block(seg(x))
#> [1] TRUE
term_jacobian_block(jump(x))
#> [1] FALSE
```
