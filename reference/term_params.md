# Parameters of a Structural Term

Returns the names of a structural term's own parameters, in the order
that its filter, its likelihood and their derivatives use. A structural
term contributes no design block, so these are not coefficients: they
are estimated beside the distribution's, on the unconstrained scale that
[`term_links()`](https://statmodels7.github.io/modelterms7/reference/term_links.md)
defines, and
[`term_npar()`](https://statmodels7.github.io/modelterms7/reference/term_npar.md)
counts them.

## Usage

``` r
term_params(term, ...)
```

## Arguments

- term:

  An object inheriting from
  [`structural_term()`](https://statmodels7.github.io/modelterms7/reference/structural_term.md).
  A class that does not implement the generic signals the error
  `"the term class 'X' does not implement term_params()."`.

- ...:

  Passed to methods. The shipped methods do not use these arguments.

## Value

A character vector, one name per parameter, of length
[`term_npar()`](https://statmodels7.github.io/modelterms7/reference/term_npar.md).

## Details

The names are the term's own vocabulary.
[`gas()`](https://statmodels7.github.io/modelterms7/reference/gas.md)
gives `omega` for the level, `kappa1` ... `kappap` for the score
loadings and `pacf1` ... `pacfq` for the persistence.
[`regime()`](https://statmodels7.github.io/modelterms7/reference/regime.md)
gives `level1`, the gaps `gap2` ... `gapk` between consecutive levels,
and the free entries `alr*.*` of its transition matrix. A marginal
break-point term gives the names listed in
[`MarginalBreakTerm()`](https://statmodels7.github.io/modelterms7/reference/MarginalBreakTerm.md).
The names index everything else about the term:
[`term_start()`](https://statmodels7.github.io/modelterms7/reference/term_start.md)
returns one value per name,
[`term_readable()`](https://statmodels7.github.io/modelterms7/reference/term_readable.md)
maps them onto the quantities that a reader reads, and the `index` of a
[`term_penalties()`](https://statmodels7.github.io/modelterms7/reference/term_penalties.md)
entry gives positions in this vector.

A subformula expands in place the parameter that it develops, so
`gas(p = 1, q = 1)` has three parameters and
`gas(p = 1, q = 1, omega ~ z)` has four: `omega.(Intercept)` and
`omega.z` where the level was.

The method on
[`structural_term()`](https://statmodels7.github.io/modelterms7/reference/structural_term.md)
signals an error that names the class. A structural class supplies its
own method.

## See also

[`term_links()`](https://statmodels7.github.io/modelterms7/reference/term_links.md)
for the chart on which each parameter is carried,
[`term_start()`](https://statmodels7.github.io/modelterms7/reference/term_start.md)
for where the parameters begin,
[`term_readable()`](https://statmodels7.github.io/modelterms7/reference/term_readable.md)
for the quantities to which they map,
[`term_coef_names()`](https://statmodels7.github.io/modelterms7/reference/term_coef_names.md)
for the additive branch's equivalent.

## Examples

``` r
# The score-driven vocabulary consists of a level, a loading and a
# persistence.
term_params(gas(p = 1, q = 1))
#> [1] "omega"  "kappa1" "pacf1" 
term_params(gas(p = 2, q = 2))
#> [1] "omega"  "kappa1" "kappa2" "pacf1"  "pacf2" 

# A subformula expands the parameter that it develops, in place.
set.seed(1)
d <- data.frame(y = rnorm(30), z = rnorm(30), t = 1:30)
term_params(term_build(gas(p = 1, q = 1, omega ~ z, time = t), d))
#> [1] "omega.(Intercept)" "omega.z"           "kappa1"           
#> [4] "pacf1"            

# On this branch term_npar() counts the names that term_params() returns.
g <- gas(p = 1, q = 2)
c(npar = term_npar(g), names = length(term_params(g)))
#>  npar names 
#>     4     4 
```
