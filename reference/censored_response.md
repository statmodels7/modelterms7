# S7 Class for a Censored Response

Holds a response whose observations are not all exact: the values, a
lower and an upper bound for each of them, and a status saying which of
the four cases the observation is.
[`cens()`](https://statmodels7.github.io/modelterms7/reference/cens.md)
is the constructor to call, and it derives the statuses from the values
and the bounds. The raw S7 constructor documented here takes all four
vectors and checks only that they agree, so it is the form to use when
the statuses are already known.

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

  A numeric vector of response values, one per observation, `NA` where
  the observation is interval-censored.

- lwr, upr:

  Numeric vectors of censoring bounds, of the same length as `y`.
  Infinite entries are allowed and mean no censoring on that side. The
  validator rejects any pair with `lwr >= upr`.

- status:

  A character vector of the same length as `y`, each entry one of
  `"observed"`, `"left"`, `"right"` or `"interval"`. Any other string
  throws, naming the first one it meets.

## Value

An S7 object of class `censored_response` with properties `y`, `lwr`,
`upr` (numeric, all of one length) and `status` (character, the same
length). It carries no methods beyond
[`print()`](https://rdrr.io/r/base/print.html).

## The four statuses and what each asserts

Writing \\Y\\ for the unobserved response, \\L\\ for `lwr` and \\U\\ for
`upr`, the four values of `status` assert

|              |                   |                         |
|--------------|-------------------|-------------------------|
| status       | what is known     | likelihood contribution |
| `"observed"` | \\Y = y\\         | \\f(y)\\                |
| `"left"`     | \\Y \le L\\       | \\F(L)\\                |
| `"right"`    | \\Y \ge U\\       | \\1 - F(U)\\            |
| `"interval"` | \\L \le Y \le U\\ | \\F(U) - F(L)\\         |

with \\f\\ and \\F\\ the density and the distribution function of the
fitted family. `y` is `NA` for an interval-censored observation, its
value being unknown; every other status carries a number.

## What reads it today

[`interpret_formula()`](https://statmodels7.github.io/modelterms7/reference/interpret_formula.md)
accepts `cens(...)` on the left of a formula and returns the object as
the `response` element of its result. Nothing then assembles the four
contributions above: `statmodels7::statmod()` stops with a message
naming the gap, and no other function in the toolkit reads the class.
The pieces exist:
[`distributions7::distrib_grad_cdf()`](https://statmodels7.github.io/distributions7/reference/distrib_grad_cdf.html)
and `distrib_hess_cdf()` carry the derivatives of \\F\\ in the
parameters. The assembler that would use them is not written, so the
class records what is known about each observation and no more.

## What the validator enforces

`y`, `lwr`, `upr` and `status` must be the same length; every lower
bound must be **strictly** below its upper bound, so `lwr = upr` is
rejected; and every status must be one of the four names. A failure
throws with the offending rule quoted.

## See also

[`cens()`](https://statmodels7.github.io/modelterms7/reference/cens.md),
which derives the statuses instead of taking them;
[`interpret_formula()`](https://statmodels7.github.io/modelterms7/reference/interpret_formula.md),
which accepts it on the left of a formula.

## Examples

``` r
# Built through cens(), which is the ordinary route.
r <- cens(c(0, 1, 5, NA), lwr = 0, upr = 5)
r@status
#> [1] "left"     "observed" "right"    "interval"
r@y                       # NA survives on the interval-censored row
#> [1]  0  1  5 NA
cbind(lwr = r@lwr, upr = r@upr)
#>      lwr upr
#> [1,]   0   5
#> [2,]   0   5
#> [3,]   0   5
#> [4,]   0   5

# The raw constructor takes statuses already known, and checks them.
censored_response(y = c(2, 9), lwr = c(-Inf, -Inf), upr = c(Inf, 9),
                  status = c("observed", "right"))
#> <censored_response> 2 observations: 1 observed, 1 right

# Three ways to fail the validator.
try(censored_response(y = 1, lwr = 2, upr = 1, status = "observed"))
#> Error : <modelterms7::censored_response> object is invalid:
#> - every lower bound must be strictly below its upper bound
try(censored_response(y = 1, lwr = 0, upr = 2, status = "cut"))
#> Error : <modelterms7::censored_response> object is invalid:
#> - unknown status 'cut'
try(censored_response(y = c(1, 2), lwr = 0, upr = 2, status = "observed"))
#> Error : <modelterms7::censored_response> object is invalid:
#> - y, lwr, upr and status must have the same length
```
