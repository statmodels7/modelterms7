# Print a Nonlinear Term

Prints the label and, for a built term, how many coefficients it
carries, which route its derivatives take, and the names of its
parameters. A specification says only that it is one.

## Arguments

- x:

  An
  [`NlTerm()`](https://statmodels7.github.io/modelterms7/reference/NlTerm.md),
  built or not.

- ...:

  Unused, and accepted so that the signature matches
  [`print()`](https://rdrr.io/r/base/print.html)'s.

## Value

`x`, invisibly. Called for the lines it writes.

## Details

The built form is

    <NlTerm> 'nl' built: 2 coefficients; symbolic derivatives
      parameters: a, r

`symbolic` means [`stats::deriv()`](https://rdrr.io/r/stats/deriv.html)
could read the expression, and `numeric` that the derivatives are
differenced. That is worth seeing at a glance: the two routes differ by
six orders of magnitude at the higher orders, which
[`nl()`](https://statmodels7.github.io/modelterms7/reference/nl.md)
tabulates.

The coefficient count is the block's width, so a developed parameter
makes it larger than the number of parameters printed on the second
line.

## See also

[`nl()`](https://statmodels7.github.io/modelterms7/reference/nl.md),
[`nl_fderiv()`](https://statmodels7.github.io/modelterms7/reference/nl_fderiv.md).

## Examples

``` r
set.seed(1)
dd <- data.frame(x = seq(0, 3, length.out = 60),
                 g = factor(rep(c("u", "v"), 30)))
dd$y <- 2 * exp(-1.3 * dd$x) + rnorm(60, sd = 0.05)

# A specification, and the same term built.
nl(~ a * exp(-r * x), start = list(a = 1, r = 1))
#> <NlTerm> 'nl' (specification; call term_build() with data)
term_build(nl(~ a * exp(-r * x), start = list(a = 1, r = 1)), dd)
#> <NlTerm> 'nl' built: 2 coefficients; symbolic derivatives
#>   parameters: a, r

# An opaque function differences its derivatives, and says so.
term_build(nl(function(x, theta) theta$a * exp(-theta$r * x),
              params = c("a", "r"), x = x, start = list(a = 1, r = 1)), dd)
#> <NlTerm> 'nl' built: 2 coefficients; numeric derivatives
#>   parameters: a, r

# Two parameters, three coefficients: `a` is developed over a factor.
term_build(nl(~ a * exp(-r * x), a ~ 0 + g, start = list(a = 1, r = 1)), dd)
#> <NlTerm> 'nl' built: 3 coefficients; symbolic derivatives
#>   parameters: a, r
```
