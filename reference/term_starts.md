# Where a Fit of a Term's Own Parameters Starts, Several Times

The starting points a fitting layer tries for a structural term's own
parameters, as a list of vectors of the kind
[`term_start()`](https://statmodels7.github.io/modelterms7/reference/term_start.md)
returns. The layer fits from each and keeps the best.

## Usage

``` r
term_starts(term, ..., target = NULL)
```

## Arguments

- term:

  A built structural term.

- ...:

  Passed to
  [`term_start()`](https://statmodels7.github.io/modelterms7/reference/term_start.md).

- target:

  The response on the scale of the predictor, or `NULL`, as
  [`term_start()`](https://statmodels7.github.io/modelterms7/reference/term_start.md)
  reads it.

## Value

A non-empty list of named numeric vectors on the unconstrained scale,
each of length
[`term_npar()`](https://statmodels7.github.io/modelterms7/reference/term_npar.md)
and named as
[`term_params()`](https://statmodels7.github.io/modelterms7/reference/term_params.md).

## Details

The base method returns one start,
[`term_start()`](https://statmodels7.github.io/modelterms7/reference/term_start.md)'s,
so a term that says nothing is fitted once, as before.
[`regime()`](https://statmodels7.github.io/modelterms7/reference/regime.md)
returns `n_start` of them: its likelihood has several maxima and the
term can say which region of its parameters is worth covering.

The first element is always
[`term_start()`](https://statmodels7.github.io/modelterms7/reference/term_start.md)'s,
so a fit with one start is the fit without this generic.

## See also

[`term_start()`](https://statmodels7.github.io/modelterms7/reference/term_start.md)
for the first start,
[`regime()`](https://statmodels7.github.io/modelterms7/reference/regime.md)
for a term that returns several.

## Examples

``` r
length(term_starts(gas(p = 1, q = 1)))
#> [1] 1
st <- term_starts(regime(k = 2, n_start = 3))
length(st)
#> [1] 3
identical(st[[1]], term_start(regime(k = 2)))
#> [1] TRUE
```
