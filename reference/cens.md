# Censored Response Constructor

Marks a response as censored, for the left-hand side of a model formula:
`cens(y, lwr = 0)` in a formula declares that values at or below the
bound are left-censored there.

## Usage

``` r
cens(y, lwr = -Inf, upr = Inf)
```

## Arguments

- y:

  A numeric vector; `NA` for interval-censored observations.

- lwr:

  A numeric vector of lower bounds, length 1 or `length(y)`. Defaults to
  `-Inf` (no left censoring).

- upr:

  A numeric vector of upper bounds, length 1 or `length(y)`. Defaults to
  `Inf` (no right censoring).

## Value

An object of class
[`censored_response`](https://statmodels7.github.io/modelterms7/reference/censored_response.md).

## Details

The bounds are recycled to the length of `y`, so a scalar bound applies
to every observation and a vector gives per-observation bounds. The
status of each observation follows from the values: an observation with
`y <= lwr` is left-censored (all that is known is \\Y \le lwr\\), one
with `y >= upr` is right-censored, one with `y` strictly inside the
bounds is observed exactly, and one with `y = NA` and both bounds finite
is interval-censored (\\Y \in \[lwr, upr\]\\). An `NA` value without two
finite bounds carries no information and is rejected.

## See also

[`interpret_formula`](https://statmodels7.github.io/modelterms7/reference/interpret_formula.md),
[`check_term`](https://statmodels7.github.io/modelterms7/reference/check_term.md)

## Examples

``` r
r <- cens(c(0, 0.7, 2.4), lwr = 0)
r@status
#> [1] "left"     "observed" "observed"
```
