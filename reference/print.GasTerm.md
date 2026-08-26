# Print a Score-Driven Term

Prints the label, the two orders and, for a built term, over how many
groups the filter runs. A second line names any developed parameters,
and a third lists the parameters themselves.

## Arguments

- x:

  A
  [`GasTerm()`](https://statmodels7.github.io/modelterms7/reference/GasTerm.md),
  built or not.

- ...:

  Unused, and accepted so that the signature matches
  [`print()`](https://rdrr.io/r/base/print.html)'s.

## Value

`x`, invisibly. Called for the lines it writes.

## Details

The form is

    <GasTerm> 'gas': score-driven, p = 1, q = 1; 3 group(s)
      developed: omega
      parameters: omega.(Intercept), omega.g2, alpha1, pacf1

The `developed` line appears only where a subformula was given. The
parameter list is
[`term_params()`](https://statmodels7.github.io/modelterms7/reference/term_params.md),
so a developed parameter shows there as its coefficients.

A built term is described by its group count rather than by the word
"built", the count being the more useful line.

## See also

[`gas()`](https://statmodels7.github.io/modelterms7/reference/gas.md),
[`term_params()`](https://statmodels7.github.io/modelterms7/reference/term_params.md).

## Examples

``` r
set.seed(1)
dd <- data.frame(t = 1:60, y = rnorm(60), id = rep(1:3, each = 20),
                 g = factor(rep(c("u", "v"), 30)))

# A specification, and the same term built over three groups.
gas(p = 1, q = 2)
#> <GasTerm> 'gas': score-driven, p = 1, q = 2 (specification)
#>   parameters: omega, alpha1, pacf1, pacf2
term_build(gas(p = 1, q = 1, by = id, time = t), dd)
#> <GasTerm> 'gas': score-driven, p = 1, q = 1; 3 group(s)
#>   parameters: omega, alpha1, pacf1

# A developed parameter is named on its own line.
term_build(gas(p = 1, q = 1, omega ~ g, by = id, time = t), dd)
#> <GasTerm> 'gas': score-driven, p = 1, q = 1; 3 group(s)
#>   developed: omega
#>   parameters: omega.(Intercept), omega.gv, alpha1, pacf1
```
