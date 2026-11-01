# The Number of Stages in Which a Term's Versions Are Searched

Returns the number of stages in which a fitting layer searches the
versions of `term`
([`term_variants()`](https://statmodels7.github.io/modelterms7/reference/term_variants.md)),
one part of the choice at each stage with the others held at the version
chosen before.

## Usage

``` r
term_variant_stages(term, ...)
```

## Arguments

- term:

  A term.

- ...:

  Passed to methods.

## Value

A single positive whole number.

## Details

The base method returns 1, a single set of versions. For
[`scs()`](https://statmodels7.github.io/modelterms7/reference/scs.md)
with a located shape without a location and a factor `by`, the method
returns the number of levels, one stage per level.

## See also

[`term_variants()`](https://statmodels7.github.io/modelterms7/reference/term_variants.md).

## Examples

``` r
term_variant_stages(s(x))
#> [1] 1

set.seed(1)
d <- data.frame(x = runif(120), g = gl(3, 40))
tm <- term_build(scs(x, shape_unimodal(), basis7::bspline_smooth(k = 8),
                     by = g), d)
term_variant_stages(tm)
#> [1] 3
```
