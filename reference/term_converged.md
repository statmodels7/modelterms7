# Has a Term's Own Iteration Settled?

`TRUE` when a term whose block is refreshed has nothing further to say
about where its own parameters are. A term whose block does not move
answers `TRUE`, having no iteration of its own.

## Usage

``` r
term_converged(term, ...)
```

## Arguments

- term:

  A built term.

- ...:

  Passed to methods.

## Value

A single logical.

## Details

It exists because a score cannot always answer the question. Where the
block is the Jacobian of the contribution –
[`nl`](https://statmodels7.github.io/modelterms7/reference/nl.md),
[`seg`](https://statmodels7.github.io/modelterms7/reference/seg.md) –
the gradient of the model's objective is the block times the derivative
of the log-likelihood in the predictor, and its vanishing is the test.
Where the block is a working LINEARIZATION with a frozen weight –
[`jump`](https://statmodels7.github.io/modelterms7/reference/jump.md),
[`jseg`](https://statmodels7.github.io/modelterms7/reference/jseg.md) –
it is not: the profile objective of a discontinuous term is a step
function in the break-point, so it has no gradient to vanish, and the
quantity the iteration actually drives to zero is the movement of the
break-point. That movement is what
[`seg_converged`](https://statmodels7.github.io/modelterms7/reference/seg_step.md)
reads, and a fitting layer asks for it here without knowing which
construction it holds.

## See also

[`term_refresh`](https://statmodels7.github.io/modelterms7/reference/term_refresh.md),
[`seg_converged`](https://statmodels7.github.io/modelterms7/reference/seg_step.md)

## Examples

``` r
dd <- data.frame(x = seq(0, 2, length.out = 20))
term_converged(term_build(linpar(~x), dd))
#> [1] TRUE
```
