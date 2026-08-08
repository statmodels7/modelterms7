# S7 Class for a Censored Response

The response object
[`cens`](https://statmodels7.github.io/modelterms7/reference/cens.md)
constructs: the observed values, the per-observation censoring bounds,
and the status each observation carries (`"observed"`, `"left"`,
`"right"` or `"interval"`). The likelihood assembler of the model layer
consumes it, contributing a density where the observation is exact and a
difference of distribution functions where it is censored.

## Usage

``` r
censored_response(
  y = integer(0),
  lwr = integer(0),
  upr = integer(0),
  status = character(0)
)
```

## Arguments

- y:

  The numeric response values (`NA` for an interval-censored
  observation).

- lwr, upr:

  The numeric censoring bounds, one value per observation.

- status:

  The character vector of per-observation statuses.

## Value

An object of class `censored_response`.

## See also

[`cens`](https://statmodels7.github.io/modelterms7/reference/cens.md)

## Examples

``` r
S7::S7_inherits(cens(c(0, 1.2), lwr = 0), censored_response)
#> [1] TRUE
```
