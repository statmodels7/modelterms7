# The Term That Stands in for a Term During a Pilot Fit

Returns the unbuilt term with which a fitting layer replaces `term` in a
pilot fit, or `NULL` when the term needs no pilot. A term returns one
when it carries a discrete choice that the data must settle (the index
of the mode of
[`shape_unimodal()`](https://statmodels7.github.io/modelterms7/reference/shape_unimodal.md)),
which the layer makes with
[`term_variants()`](https://statmodels7.github.io/modelterms7/reference/term_variants.md)
from the contribution that the stand-in fits.

## Usage

``` r
term_pilot(term, ...)
```

## Arguments

- term:

  A term, built or not.

- ...:

  Passed to methods.

## Value

`NULL`, or an unbuilt term carrying the label of `term`.

## See also

[`term_variants()`](https://statmodels7.github.io/modelterms7/reference/term_variants.md),
which reads the pilot.

## Examples

``` r
term_pilot(s(x))
#> NULL
class(term_pilot(scs(x, shape_unimodal())))[1]
#> [1] "modelterms7::SmoothTerm"

# A mode given in the shape needs no pilot.
term_pilot(scs(x, shape_unimodal(mode = 0.4)))
#> NULL
```
