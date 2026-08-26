# Print a Censored Response

Prints one line giving the number of observations and how many carry
each status, in the fixed order `observed`, `left`, `right`, `interval`.
A status no observation has is left out, so a response with no censoring
at all prints as `<censored_response> 20 observations: 20 observed`. The
values and the bounds are not shown; read them from the `y`, `lwr` and
`upr` properties.

## Arguments

- x:

  A
  [`censored_response()`](https://statmodels7.github.io/modelterms7/reference/censored_response.md)
  object.

- ...:

  Unused, and accepted so that the signature matches
  [`print()`](https://rdrr.io/r/base/print.html)'s.

## Value

`x`, invisibly. Called for the line it writes to the console.

## See also

[`cens()`](https://statmodels7.github.io/modelterms7/reference/cens.md),
which builds the object and assigns the statuses.

## Examples

``` r
cens(c(0, 1, 5, NA), lwr = 0, upr = 5)
#> <censored_response> 4 observations: 1 observed, 1 left, 1 right, 1 interval

# Only the statuses present are listed.
cens(c(1, 2, 3))
#> <censored_response> 3 observations: 3 observed
```
