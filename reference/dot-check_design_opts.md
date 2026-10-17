# Check a Term's Storage and Contrasts

Validates the `sparse` and `contrasts` arguments that the
design-building constructors share, and normalizes `contrasts = NULL` to
an empty list so that the class property is always a list. The
constructor calls it, so an invalid value is rejected at the call that
wrote it.

## Usage

``` r
.check_design_opts(sparse, contrasts, what = "this term")
```

## Arguments

- sparse:

  What the constructor was given: `TRUE`, `FALSE`, or `NULL` for the
  storage to be settled at build. Any other value signals an error whose
  message names `what`.

- contrasts:

  What the constructor was given: a named list, or `NULL`. Any other
  value signals an error. The names are not checked against the
  formula's factors; at build,
  [`stats::model.matrix()`](https://rdrr.io/r/stats/model.matrix.html)
  warns about a name that the data do not carry and ignores that entry.

- what:

  The constructor's name, used in the two messages.

## Value

A list of two: `sparse` unchanged, and `contrasts` as given or an empty
list.

## See also

[`linpar()`](https://statmodels7.github.io/modelterms7/reference/linpar.md)
and
[`penalized_terms()`](https://statmodels7.github.io/modelterms7/reference/penalized_terms.md),
the constructors that call it.
