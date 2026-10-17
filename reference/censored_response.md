# S7 Class for a Censored Response

Holds a response whose observations are not all exact: the values, a
lower and an upper bound for each of them, and a status that indicates
which of the four cases applies to the observation.
[`cens()`](https://statmodels7.github.io/modelterms7/reference/cens.md)
is the usual constructor and derives the statuses from the values and
the bounds. The raw S7 constructor documented here takes all four
vectors and validates them (equal lengths, bounds that are not missing,
each lower bound strictly below its upper bound, and status names), so
it is the form to use when the statuses are already known.

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
  validator rejects a missing entry and any pair with `lwr >= upr`.

- status:

  A character vector of the same length as `y`, each entry one of
  `"observed"`, `"left"`, `"right"` or `"interval"`. Any other string
  signals an error that names the first one found.

## Value

An S7 object of class `censored_response` with properties `y`, `lwr`,
`upr` (numeric, all of one length) and `status` (character, the same
length). The class has a [`print()`](https://rdrr.io/r/base/print.html)
method and no other methods.

## The four statuses and what each indicates

Writing \\Y\\ for the unobserved response, \\L\\ for `lwr` and \\U\\ for
`upr`, the four values of `status` indicate the following.

|              |                   |                         |
|--------------|-------------------|-------------------------|
| status       | what is known     | likelihood contribution |
| `"observed"` | \\Y = y\\         | \\f(y)\\                |
| `"left"`     | \\Y \le L\\       | \\F(L)\\                |
| `"right"`    | \\Y \ge U\\       | \\1 - F(U)\\            |
| `"interval"` | \\L \le Y \le U\\ | \\F(U) - F(L)\\         |

with \\f\\ and \\F\\ the density and the distribution function of the
fitted family. `y` is `NA` for an interval-censored observation because
its value is unknown; every other status carries a number.

## Use in fitting

[`interpret_formula()`](https://statmodels7.github.io/modelterms7/reference/interpret_formula.md)
accepts `cens(...)` on the left of a formula and returns the object as
the `response` element of its result. The toolkit does not assemble the
four contributions above into a likelihood:
[`statmodels7::statmod()`](https://statmodels7.github.io/statmodels7/reference/statmod.html)
and
[`statmodels7::rstatmod()`](https://statmodels7.github.io/statmodels7/reference/rstatmod.html)
reject a censored response with an error message, and the other
functions of the toolkit do not read the class. The derivatives of \\F\\
in the parameters are available from
[`distributions7::distrib_grad_cdf()`](https://statmodels7.github.io/distributions7/reference/distrib_grad_cdf.html)
and `distrib_hess_cdf()`. Without an assembler, the class stores the
values, the bounds and the statuses, and it holds no other information.

## What the validator enforces

`y`, `lwr`, `upr` and `status` must be the same length; `lwr` and `upr`
must not be missing; every lower bound must be strictly below its upper
bound, so `lwr = upr` is rejected; and every status must be one of the
four names. A failure signals an error that quotes the violated rule.
The validator does not check that a status agrees with the values: a
status of `"interval"` with a non-missing `y` is accepted, and so is a
missing `y` with another status.

## See also

[`cens()`](https://statmodels7.github.io/modelterms7/reference/cens.md),
which derives the statuses instead of taking them;
[`interpret_formula()`](https://statmodels7.github.io/modelterms7/reference/interpret_formula.md),
which accepts it on the left of a formula.

## Examples

``` r
# The object is built through cens(), which is the ordinary route.
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

# The next four calls fail the validator in four different ways.
try(censored_response(y = 1, lwr = NA_real_, upr = 1, status = "observed"))
#> Error : <modelterms7::censored_response> object is invalid:
#> - lwr and upr must not be missing
try(censored_response(y = 1, lwr = 2, upr = 1, status = "observed"))
#> Error : <modelterms7::censored_response> object is invalid:
#> - every lower bound must be strictly below its upper bound
try(censored_response(y = 1, lwr = 0, upr = 2, status = "cut"))
#> Error : <modelterms7::censored_response> object is invalid:
#> - the status 'cut' is not one of the four names
try(censored_response(y = c(1, 2), lwr = 0, upr = 2, status = "observed"))
#> Error : <modelterms7::censored_response> object is invalid:
#> - y, lwr, upr and status must have the same length
```
