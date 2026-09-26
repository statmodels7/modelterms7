# Has a Term's Own Iteration Run Out of Step Control?

`TRUE` when a term whose block is refreshed has not settled and its
iteration has no means left of bringing it to rest, so further refreshes
cannot be expected to settle it. A term with no iteration of its own
answers `FALSE`.

## Usage

``` r
term_stalled(term, ...)
```

## Arguments

- term:

  A built term.

- ...:

  Passed to methods.

## Value

A single logical.

## Details

The discontinuous break-point terms,
[`jump()`](https://statmodels7.github.io/modelterms7/reference/jump.md)
and
[`jseg()`](https://statmodels7.github.io/modelterms7/reference/jseg.md),
control their step with the scaling factor of Fasola, Muggeo and
Kuchenhoff, halved whenever a break-point reverses direction and held at
a floor derived from the conditioning of the working block (see
[`term_refresh()`](https://statmodels7.github.io/modelterms7/reference/term_refresh.md)).
Once every break-point that has not settled sits at that floor the
factor can shrink no further, and the iteration continues without step
control. The profile objective of a discontinuous term is constant
between two consecutive observations, so a break-point in that state can
pass from one observation to the next indefinitely. A fitting layer
reads this to end the working phase and report the fit as not converged.

A term whose block is the Jacobian of its contribution,
[`seg()`](https://statmodels7.github.io/modelterms7/reference/seg.md)
and every smoothed construction among them, has no scaling schedule and
answers `FALSE`.

## See also

[`term_converged()`](https://statmodels7.github.io/modelterms7/reference/term_converged.md),
[`term_refresh()`](https://statmodels7.github.io/modelterms7/reference/term_refresh.md)

## Examples

``` r
dd <- data.frame(x = seq(0, 2, length.out = 20))
term_stalled(term_build(linpar(~x), dd))
#> [1] FALSE
```
