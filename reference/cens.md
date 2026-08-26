# Mark a Response as Censored

Builds a
[`censored_response()`](https://statmodels7.github.io/modelterms7/reference/censored_response.md)
from a vector of values and one or two censoring bounds, deriving each
observation's status from where its value falls. `cens(y, lwr = 0)` on
the left of a model formula declares that any value at or below zero is
left-censored there; `cens(y, upr = 8)` does the same at the top; giving
both bounds allows either, and an `NA` value between two finite bounds
is interval-censored.

## Usage

``` r
cens(y, lwr = -Inf, upr = Inf)
```

## Arguments

- y:

  A numeric vector of responses, `NA` at an interval-censored
  observation. Coerced with
  [`as.numeric()`](https://rdrr.io/r/base/numeric.html), so an integer
  vector is accepted and a factor is not.

- lwr:

  A numeric lower bound, of length 1 or `length(y)`, recycled to
  `length(y)`. `-Inf` by default, which is no left censoring. Any other
  length throws `"'lwr' must have length 1 or length(y)."`, and an `NA`
  entry throws `"'lwr' must not contain NA."`.

- upr:

  A numeric upper bound, on the same terms, `Inf` by default. Every
  `upr` must be strictly above its `lwr`; equal bounds throw from the
  class validator.

## Value

A
[`censored_response()`](https://statmodels7.github.io/modelterms7/reference/censored_response.md)
object of length `length(y)`, carrying `y` unchanged (`NA` included),
the recycled `lwr` and `upr`, and the derived `status`.

## The rule that assigns a status

The bounds are recycled to `length(y)` first, so a scalar applies to
every observation and a vector gives one bound per observation. Each
entry is then classified by where it sits, with the bounds inclusive:

|                                |              |                                |
|--------------------------------|--------------|--------------------------------|
| condition                      | status       | what is known                  |
| `y <= lwr`                     | `"left"`     | \\Y \le\\ `lwr`                |
| `y >= upr`                     | `"right"`    | \\Y \ge\\ `upr`                |
| `lwr < y < upr`                | `"observed"` | \\Y = y\\                      |
| `is.na(y)`, both bounds finite | `"interval"` | \\Y \in \[\\`lwr`, `upr`\\\]\\ |

The tests are applied in that order, so a value at or below `lwr` is
left- censored even when it is also at or above `upr`; the validator
forbids `lwr >= upr`, so the two cannot both bind. An `NA` value without
two finite bounds says nothing about \\Y\\ at all and throws.

At the defaults `lwr = -Inf` and `upr = Inf` no value can reach a bound,
every status comes back `"observed"`, and the object carries the same
information the bare vector does.

## Where it can be used

On the left of a formula passed to
[`interpret_formula()`](https://statmodels7.github.io/modelterms7/reference/interpret_formula.md),
which returns the object as its `response` element.
[`statmodels7::statmod()`](https://statmodels7.github.io/statmodels7/reference/statmod.html)
refuses that response with a message naming the gap: the toolkit marks
censoring and does not yet assemble a censored likelihood from it. See
[`censored_response()`](https://statmodels7.github.io/modelterms7/reference/censored_response.md)
for the four contributions such an assembler would need.

## See also

[`censored_response()`](https://statmodels7.github.io/modelterms7/reference/censored_response.md)
for the class and the four likelihood contributions;
[`interpret_formula()`](https://statmodels7.github.io/modelterms7/reference/interpret_formula.md)
for the formula it goes into;
[`distributions7::distrib_grad_cdf()`](https://statmodels7.github.io/distributions7/reference/distrib_grad_cdf.html)
for the distribution-function derivatives a censored likelihood is built
from.

## Examples

``` r
# Left censoring at zero: the first value is at the bound, so it is
# censored there and the other two are exact.
r <- cens(c(0, 0.7, 2.4), lwr = 0)
r@status
#> [1] "left"     "observed" "observed"

# One of each status. The NA is interval-censored between the bounds.
r2 <- cens(c(0, 1, 5, NA), lwr = 0, upr = 5)
r2
#> <censored_response> 4 observations: 1 observed, 1 left, 1 right, 1 interval
data.frame(y = r2@y, lwr = r2@lwr, upr = r2@upr, status = r2@status)
#>    y lwr upr   status
#> 1  0   0   5     left
#> 2  1   0   5 observed
#> 3  5   0   5    right
#> 4 NA   0   5 interval

# Per-observation bounds: only the rows whose own bound binds are censored.
cens(c(1, 2, 3), lwr = c(-Inf, 2, -Inf), upr = c(Inf, Inf, 3))@status
#> [1] "observed" "left"     "right"   

# With no bounds given, every observation is exact.
all(cens(rnorm(20))@status == "observed")
#> [1] TRUE

# An NA with no finite pair of bounds carries no information.
try(cens(c(1, NA)))
#> Error : an NA response is interval-censored and needs finite 'lwr' and 'upr'.

# It goes on the left of a formula.
d <- data.frame(t = c(1, 5, 9, 2), x = c(1, 2, 3, 4))
names(interpret_formula(cens(t, upr = 8) ~ x, d))
#> [1] "response"  "terms"     "intercept" "formula"  
```
