# The Versions of a Term That a Fitting Layer Tries

Returns the unbuilt versions of `term` among which a fitting layer
chooses by its criterion, given the contribution of the term's pilot
(see
[`term_pilot()`](https://statmodels7.github.io/modelterms7/reference/term_pilot.md))
at the rows on which `term` was built. The layer fits the model with
each and keeps the best.

## Usage

``` r
term_variants(term, pilot = NULL, ...)
```

## Arguments

- term:

  A built term.

- pilot:

  A numeric vector, the contribution of the pilot at the rows on which
  `term` was built, or `NULL`.

- ...:

  Passed to methods; `stage` and `current` for a term searched in
  stages.

## Value

A named list of unbuilt terms carrying the label of `term`, non-empty at
the first stage.

## Details

The base method returns `term` alone. For
[`scs()`](https://statmodels7.github.io/modelterms7/reference/scs.md)
with
[`shape_unimodal()`](https://statmodels7.github.io/modelterms7/reference/shape_located.md)
or
[`shape_sigmoid()`](https://statmodels7.github.io/modelterms7/reference/shape_located.md)
and a location given, the versions are the term from its two starts,
named `"start: sides"` and `"start: spline"` (see
[`shape_start()`](https://statmodels7.github.io/modelterms7/reference/shape_start.md)).
With no location, the pilot's extremum (unimodal) or steepest point
(sigmoid) gives the Greville abscissa \\\xi\_{M_0}\\ nearest to it, and
the versions are the shape located at the indices \\M_0 - 2, \dots,
M_0 + 2\\ that exist, named `"M = j"`.

A term whose choice has several parts is searched in stages, one part at
a time
([`term_variant_stages()`](https://statmodels7.github.io/modelterms7/reference/term_variant_stages.md)).
The layer then passes `stage = k` and `current`, the version chosen at
the stage before, through `...`. For
[`scs()`](https://statmodels7.github.io/modelterms7/reference/scs.md)
with a factor `by` the parts are the levels: stage \\k\\ moves the index
of level \\k\\ around the index nearest to the pilot's location on the
rows of that level, with the other levels at their indices in `current`
(at the first stage, at their pilot indices), and the versions are named
`"<level>: M = j"`. The index that `current` already holds for level
\\k\\ is left out, because the layer fitted it at the stage before.

## See also

[`term_pilot()`](https://statmodels7.github.io/modelterms7/reference/term_pilot.md),
[`term_variant_stages()`](https://statmodels7.github.io/modelterms7/reference/term_variant_stages.md).

## Examples

``` r
set.seed(1)
d <- data.frame(x = runif(100))
tm <- term_build(scs(x, shape_unimodal(), basis7::bspline_smooth(k = 10)), d)
names(term_variants(tm, pilot = -(d$x - 0.3)^2))
#> [1] "M = 2" "M = 3" "M = 4" "M = 5" "M = 6"
```
