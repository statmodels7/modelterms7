# The Starting Hyperparameters That a Term Asks the Outer Search For

Returns the values at which a fitting layer starts the outer search for
the hyperparameters of `term`, or `NULL` for the layer's own start. A
held value is not a start and always takes precedence.

## Usage

``` r
term_hyper_start(term, ...)
```

## Arguments

- term:

  A term.

- ...:

  Passed to methods.

## Value

`NULL`, or a named numeric vector of starting values on the scale of the
hyperparameters.

## Details

The base method returns `NULL`. For
[`scs()`](https://statmodels7.github.io/modelterms7/reference/scs.md) it
returns the smoothing parameter that a version of
[`term_variants()`](https://statmodels7.github.io/modelterms7/reference/term_variants.md)
carries: a unimodal or sigmoid shape with its location given is fitted
from two starts of the coefficients and two starting smoothing
parameters, 1 and 0.01, because far from the mode of the data the fit
from a smoothing parameter of 1 stays on the flat curve.

## See also

[`term_variants()`](https://statmodels7.github.io/modelterms7/reference/term_variants.md).

## Examples

``` r
term_hyper_start(s(x))
#> NULL
```
