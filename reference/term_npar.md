# Number of Parameters of a Built Term

How many parameters of its own a built term carries: the number of
columns of the design block for an additive term, and the number of
entries of
[`term_params()`](https://statmodels7.github.io/modelterms7/reference/term_params.md)
for a structural one, which has no block. It is the length of the vector
a
[`term_penalties()`](https://statmodels7.github.io/modelterms7/reference/term_penalties.md)
entry indexes into, the length
[`edf()`](https://statmodels7.github.io/modelterms7/reference/edf.md)
measures against, and the length a fit reserves for the term.

## Usage

``` r
term_npar(term, ...)
```

## Arguments

- term:

  A built term (see
  [`term_build()`](https://statmodels7.github.io/modelterms7/reference/term_build.md)).
  An unbuilt additive term throws
  `"the term has not been built; call term_build(term, data) first."`.

- ...:

  Passed to methods. No shipped method reads anything here.

## Value

A single whole number.

## Details

The two methods are the two branches. On
[`additive_term()`](https://statmodels7.github.io/modelterms7/reference/additive_term.md)
it is `ncol(term@X)`, so it equals `length(term_coef_names(term))` and a
specification throws. On
[`structural_term()`](https://statmodels7.github.io/modelterms7/reference/structural_term.md)
it is `length(term_params(term))`, which is the count of the term's own
parameters after any subformula has expanded: `gas(p = 1, q = 1)` has
three, and `gas(p = 1, q = 1, omega ~ z)` has four, the level's
intercept and slope in place of the level.

## See also

[`term_coef_names()`](https://statmodels7.github.io/modelterms7/reference/term_coef_names.md)
and
[`term_params()`](https://statmodels7.github.io/modelterms7/reference/term_params.md)
for the names behind the count,
[`term_matrix()`](https://statmodels7.github.io/modelterms7/reference/term_matrix.md)
for the block,
[`edf()`](https://statmodels7.github.io/modelterms7/reference/edf.md)
for what the term spends of it.

## Examples

``` r
d <- data.frame(x = rnorm(20), z = rnorm(20),
                g = factor(rep(letters[1:4], 5)))

# An additive term counts columns.
b <- term_build(linpar(~ x + g), d)
c(npar = term_npar(b), names = length(term_coef_names(b)),
  cols = ncol(term_matrix(b)))
#>  npar names  cols 
#>     5     5     5 

# A structural term counts its own parameters, and a subformula
# replaces one of them by the coefficients developing it.
term_npar(term_build(gas(p = 1, q = 1), d))
#> [1] 3
gz <- term_build(gas(p = 1, q = 1, omega ~ z), d)
c(npar = term_npar(gz), params = length(term_params(gz)))
#>   npar params 
#>      4      4 
term_params(gz)
#> [1] "omega.(Intercept)" "omega.z"           "alpha1"           
#> [4] "pacf1"            
```
