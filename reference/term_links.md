# Links of a Structural Term's Parameters

One linkfunctions7 link per parameter of
[`term_params()`](https://statmodels7.github.io/modelterms7/reference/term_params.md),
carrying that parameter from its own admissible set onto the whole real
line. A fitting layer optimizes on that unconstrained scale, so an
optimizer never has to be told about the constraint.

## Usage

``` r
term_links(term, ...)
```

## Arguments

- term:

  An object inheriting from
  [`structural_term()`](https://statmodels7.github.io/modelterms7/reference/structural_term.md).
  A class that does not implement the generic throws
  `"the term class 'X' does not implement term_links()."`.

- ...:

  Passed to methods. No shipped method reads anything here.

## Value

A named list of linkfunctions7 link objects, one per entry of
[`term_params()`](https://statmodels7.github.io/modelterms7/reference/term_params.md)
and named by it.

## Details

The charts are chosen so that a proposal from anywhere lands somewhere
admissible.
[`gas()`](https://statmodels7.github.io/modelterms7/reference/gas.md)
puts its level on the identity, its score loadings on the **log** so a
loading stays positive, and its persistence coordinates on the
**rhobit** so each partial autocorrelation stays inside \\(-1, 1)\\ and
the filter stays stationary. The `links` argument of a constructor
overrides them per parameter.

The persistence is the case worth understanding. The stationary region
in the autoregressive coefficients is not a box, so no collection of
scalar links covers it; the partial autocorrelations are each in \\(-1,
1)\\ independently, and Levinson-Durbin carries them onto the
coefficients. That is why
[`term_readable()`](https://statmodels7.github.io/modelterms7/reference/term_readable.md)
exists: the coordinate and the quantity a reader reads are different
things.

Where a parameter carries a subformula the link is applied **inside**
the development, \\\psi\_{j,t} = g_j^{-1}(z_t'\gamma_j)\\, so the
parameter stays admissible at every observation and the coefficients
\\\gamma_j\\ are unconstrained.

The method on
[`structural_term()`](https://statmodels7.github.io/modelterms7/reference/structural_term.md)
throws, naming the class.

## See also

[`term_params()`](https://statmodels7.github.io/modelterms7/reference/term_params.md)
for the names,
[`term_start()`](https://statmodels7.github.io/modelterms7/reference/term_start.md)
for the point on this scale a fit begins at,
[`term_readable()`](https://statmodels7.github.io/modelterms7/reference/term_readable.md)
for the quantities the coordinates map to.

## Examples

``` r
# Identity for the level, log for the loading, rhobit for the persistence.
vapply(term_links(gas(p = 1, q = 2)), function(l) l@link_name, character(1))
#>      omega     alpha1      pacf1      pacf2 
#> "identity"      "log"   "rhobit"   "rhobit" 

# Each carries its parameter's own set onto the whole line.
lk <- term_links(gas(p = 1, q = 1))
vapply(lk, function(l) paste(l@link_bounds, collapse = ", "), character(1))
#>       omega      alpha1       pacf1 
#> "-Inf, Inf"    "0, Inf"     "-1, 1" 

# So any coordinate at all gives an admissible parameter. At a
# coordinate of 40 the persistence prints as 1 and is not: the
# inverse link is clamped strictly inside its bounds.
rho <- linkfunctions7::linkinv(lk$pacf1, c(-40, 0, 40))
all(rho > -1 & rho < 1)
#> [1] TRUE
1 - rho[3]
#> [1] 2.220446e-16
linkfunctions7::linkinv(lk$alpha1, c(-40, 0, 40))
#> [1] 4.248354e-18 1.000000e+00 2.353853e+17
```
