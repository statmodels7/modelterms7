# Print a Score-Driven Term

Prints the label, the two orders and, for a built term, over how many
groups the filter runs; an unbuilt term is marked "(specification)".
Where parameters are developed, a line names them; a non-zero scaling
has a line of its own with the power of the information; the parameters
themselves are listed on the last line.

## Arguments

- x:

  A
  [`GasTerm()`](https://statmodels7.github.io/modelterms7/reference/GasTerm.md),
  built or not.

- ...:

  Unused, and accepted so that the signature matches
  [`print()`](https://rdrr.io/r/base/print.html)'s.

## Value

`x`, invisibly. The function is called for the lines that it writes.

## Details

The form is

    <GasTerm> 'gas': score-driven, p = 1, q = 1; 3 group(s)
      developed: omega
      parameters: omega.(Intercept), omega.gv, kappa1, pacf1

The `developed` line appears only where a subformula was given. The
parameter list is
[`term_params()`](https://statmodels7.github.io/modelterms7/reference/term_params.md),
so a developed parameter shows there as its coefficients.

## See also

[`gas()`](https://statmodels7.github.io/modelterms7/reference/gas.md),
[`term_params()`](https://statmodels7.github.io/modelterms7/reference/term_params.md).

## Examples

``` r
set.seed(1)
dd <- data.frame(t = 1:60, y = rnorm(60), id = rep(1:3, each = 20),
                 g = factor(rep(c("u", "v"), 30)))

# The first call shows a specification, and the second the same term built
# over three groups.
gas(p = 1, q = 2)
#> <GasTerm> 'gas': score-driven, p = 1, q = 2 (specification)
#>   parameters: omega, kappa1, pacf1, pacf2
term_build(gas(p = 1, q = 1, by = id, time = t), dd)
#> <GasTerm> 'gas': score-driven, p = 1, q = 1; 3 group(s)
#>   parameters: omega, kappa1, pacf1

# A developed parameter is named on its own line.
term_build(gas(p = 1, q = 1, omega ~ g, by = id, time = t), dd)
#> <GasTerm> 'gas': score-driven, p = 1, q = 1; 3 group(s)
#>   developed: omega
#>   parameters: omega.(Intercept), omega.gv, kappa1, pacf1
```
