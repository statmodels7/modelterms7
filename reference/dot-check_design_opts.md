# Check a Term's Storage and Contrasts

Validates the `sparse` and `contrasts` arguments the design-building
constructors share, and normalizes `contrasts = NULL` to an empty list
so that the class property is always a list. Called from the
constructor, so a mistake is reported where it was written.

## Usage

``` r
.check_design_opts(sparse, contrasts, what = "this term")
```

## Arguments

- sparse:

  What the constructor was given: `TRUE`, `FALSE`, or `NULL` for the
  storage to be settled at build. Anything else throws, the message
  naming `what`.

- contrasts:

  What the constructor was given: a named list, or `NULL`. Anything else
  throws. The names are not checked against the formula's factors;
  [`stats::model.matrix()`](https://rdrr.io/r/stats/model.matrix.html)
  does that at build.

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
