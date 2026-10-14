# Has a Term's Own Iteration Settled?

`TRUE` when the iteration that a term carries on its own account has
settled, that is, when its own parameters no longer move. A term with no
such iteration, including every term whose block does not move and
[`nl()`](https://statmodels7.github.io/modelterms7/reference/nl.md),
answers `TRUE` through the base method.

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
the gradient of the model's objective is the transpose of the block
applied to the derivative of the log-likelihood in the predictor, and
the fitting layer tests its vanishing; the term itself has no verdict to
add.

Where the block is a working linearization with a frozen weight, as it
is for a sharp
[`jump()`](https://statmodels7.github.io/modelterms7/reference/jump.md)
or
[`jseg()`](https://statmodels7.github.io/modelterms7/reference/jseg.md),
there is no such gradient: the profile objective of a discontinuous term
is a step function in the break-point. What the iteration drives to zero
there is the movement of the break-point itself. A break-point term of
any construction answers through
[`seg_converged()`](https://statmodels7.github.io/modelterms7/reference/seg_step.md),
which reads that movement, and a held term answers `TRUE`. A fitting
layer calls this generic for each refreshed term, whatever its
construction.

## See also

[`term_refresh()`](https://statmodels7.github.io/modelterms7/reference/term_refresh.md),
[`seg_converged()`](https://statmodels7.github.io/modelterms7/reference/seg_step.md)

## Examples

``` r
dd <- data.frame(x = seq(0, 2, length.out = 20))
term_converged(term_build(linpar(~x), dd))
#> [1] TRUE
```
