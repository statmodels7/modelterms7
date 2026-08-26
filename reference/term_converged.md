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
block is the Jacobian of the contribution, as it is for
[`nl()`](https://statmodels7.github.io/modelterms7/reference/nl.md) and
[`seg()`](https://statmodels7.github.io/modelterms7/reference/seg.md),
the gradient of the model's objective is the block times the derivative
of the log-likelihood in the predictor, and its vanishing is the test.

Where the block is a working linearization with a frozen weight, as it
is for
[`jump()`](https://statmodels7.github.io/modelterms7/reference/jump.md)
and
[`jseg()`](https://statmodels7.github.io/modelterms7/reference/jseg.md),
there is no such gradient: the profile objective of a discontinuous term
is a step function in the break-point. What the iteration drives to zero
there is the movement of the break-point itself. That movement is what
[`seg_converged()`](https://statmodels7.github.io/modelterms7/reference/seg_step.md)
reads, and a fitting layer asks for it here without knowing which
construction it holds.

## See also

[`term_refresh()`](https://statmodels7.github.io/modelterms7/reference/term_refresh.md),
[`seg_converged()`](https://statmodels7.github.io/modelterms7/reference/seg_step.md)

## Examples

``` r
dd <- data.frame(x = seq(0, 2, length.out = 20))
term_converged(term_build(linpar(~x), dd))
#> [1] TRUE
```
