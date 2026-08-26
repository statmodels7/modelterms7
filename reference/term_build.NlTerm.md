# Build a Nonlinear Term

Resolves the parameters, their links and their subformulas against the
data, evaluates the function at the starting values, and takes the
Jacobian there as the design block. It also chooses how each order of
derivative will be obtained and checks that any function the caller
supplied returns the right components.

## Arguments

- term:

  An unbuilt or built
  [`NlTerm()`](https://statmodels7.github.io/modelterms7/reference/NlTerm.md).

- data:

  A data frame carrying the covariates the function names and whatever
  the subformulas name.

- ...:

  Unused.

## Value

The term with `X`, `coef_names`, `blueprint` and, where a subformula
brought one, `penalty` filled.

## What is settled here

On the formula route
[`stats::deriv()`](https://rdrr.io/r/stats/deriv.html) is tried on the
expression, and `deriv_mode` records whether it succeeded. Every
function given to `gradient`, `hessian`, `deriv3` or `deriv4` is called
once and its names checked against the components of that order: a name
that is not one of them, a missing component or a repeated one throws
here. Falling back to the numerical route would leave an exact
derivative silently unused, which is worse than not having it.

Each parameter's subformula goes through
[`interpret_formula()`](https://statmodels7.github.io/modelterms7/reference/interpret_formula.md)
and its terms are built, so their blueprints are recorded and reapplied
at prediction. A structural sub-term, and one whose own block moves with
its coefficients, are rejected: a parameter's submodel must be a fixed
design.

## The starting point

The block is the Jacobian **at the starting coefficients**, so it is
only as good as they are. Where `start` names nothing,
[`term_coef_start()`](https://statmodels7.github.io/modelterms7/reference/term_coef_start.md)
estimates the parameters from the data over a deterministic grid on each
free parameter's chart; zero is a degenerate point for a nonlinear
function, not a neutral one.

## See also

[`nl()`](https://statmodels7.github.io/modelterms7/reference/nl.md),
[`term_refresh()`](https://statmodels7.github.io/modelterms7/reference/term_refresh.md)
for the block at other coefficients,
[`term_predict.NlTerm()`](https://statmodels7.github.io/modelterms7/reference/term_predict.NlTerm.md)
for it at other rows.

## Examples

``` r
set.seed(1)
dd <- data.frame(x = seq(0, 3, length.out = 60),
                 g = factor(rep(c("u", "v"), 30)))
dd$y <- 2 * exp(-1.3 * dd$x) + rnorm(60, sd = 0.05)

b <- term_build(nl(~ a * exp(-r * x), start = list(a = 1, r = 1)), dd)
term_coef_names(b)
#> [1] "nl.a" "nl.r"
dim(term_matrix(b))
#> [1] 60  2

# A developed parameter contributes one coefficient per column of its
# own design.
bs <- term_build(nl(~ a * exp(-r * x), a ~ 0 + g,
                    start = list(a = 1, r = 1)), dd)
term_coef_names(bs)
#> [1] "nl.a.gu" "nl.a.gv" "nl.r"   

# A supplied derivative is checked here, by name.
bad <- function(theta, data) list(a = data$x, wrong = data$x)
try(term_build(nl(~ a * exp(-r * x), gradient = bad,
                  start = list(a = 1, r = 1)), dd))
#> Error : 'gradient' returned the component 'wrong', which is not one of this term's.
#>   Expected: a, r.
```
