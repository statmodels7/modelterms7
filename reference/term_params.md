# Parameters of a Structural Term

The names of a structural term's own parameters, in the order its filter
and its derivative recursions expect them. A structural term contributes
no design block, so these are not coefficients: they are estimated
beside the distribution's, on the unconstrained scale
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
  A class that does not implement the generic throws
  `"the term class 'X' does not implement term_params()."`.

- ...:

  Passed to methods. No shipped method reads anything here.

## Value

A character vector, one name per parameter, of length
[`term_npar()`](https://statmodels7.github.io/modelterms7/reference/term_npar.md).

## Details

The names are the term's own vocabulary.
[`gas()`](https://statmodels7.github.io/modelterms7/reference/gas.md)
answers `omega` for the level, `alpha1` ... `alphap` for the score
loadings and `pacf1` ... `pacfq` for the persistence;
[`regime()`](https://statmodels7.github.io/modelterms7/reference/regime.md)
answers its levels and the free entries of its transition matrix. They
are what indexes everything else about the term:
[`term_start()`](https://statmodels7.github.io/modelterms7/reference/term_start.md)
returns one value per name,
[`term_readable()`](https://statmodels7.github.io/modelterms7/reference/term_readable.md)
carries them onto the quantities a reader reads, and a
[`term_penalties()`](https://statmodels7.github.io/modelterms7/reference/term_penalties.md)
entry's `index` gives positions in this vector.

A **subformula** expands the parameter it develops in place, so
`gas(p = 1, q = 1)` has three parameters and
`gas(p = 1, q = 1, omega ~ z)` has four: `omega.(Intercept)` and
`omega.z` where the level was.

The method on
[`structural_term()`](https://statmodels7.github.io/modelterms7/reference/structural_term.md)
throws, naming the class: a structural class supplies this itself.

## See also

[`term_links()`](https://statmodels7.github.io/modelterms7/reference/term_links.md)
for the chart each rides,
[`term_start()`](https://statmodels7.github.io/modelterms7/reference/term_start.md)
for where they begin,
[`term_readable()`](https://statmodels7.github.io/modelterms7/reference/term_readable.md)
for the quantities they map to,
[`term_coef_names()`](https://statmodels7.github.io/modelterms7/reference/term_coef_names.md)
for the additive branch's equivalent.

## Examples

``` r
# The score-driven vocabulary: a level, a loading, a persistence.
term_params(gas(p = 1, q = 1))
#> [1] "omega"  "alpha1" "pacf1" 
term_params(gas(p = 2, q = 2))
#> [1] "omega"  "alpha1" "alpha2" "pacf1"  "pacf2" 

# A subformula expands the parameter it develops, in place.
set.seed(1)
d <- data.frame(y = rnorm(30), z = rnorm(30), t = 1:30)
term_params(term_build(gas(p = 1, q = 1, omega ~ z, time = t), d))
#> [1] "omega.(Intercept)" "omega.z"           "alpha1"           
#> [4] "pacf1"            

# It is what term_npar() counts on this branch.
g <- gas(p = 1, q = 2)
c(npar = term_npar(g), names = length(term_params(g)))
#>  npar names 
#>     4     4 
```
