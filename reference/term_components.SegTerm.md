# How a Break-Point Term's Columns Divide

One entry per coefficient of the term, `beta`, `gamma1` ... , `delta1`
... , `psi1` ... , giving the columns that coefficient owns and the
sub-terms developing it. A coefficient with no development owns one
column; a developed one owns as many as its own design has.

## Arguments

- term:

  A built
  [`SegTerm()`](https://statmodels7.github.io/modelterms7/reference/SegTerm.md).
  An unbuilt one gives an empty list.

- ...:

  Unused.

## Value

A named list, one entry per coefficient and named by it, each with
`name`, `index`, `subs` and `sub_index` as
[`term_components()`](https://statmodels7.github.io/modelterms7/reference/term_components.md)
describes.

## Details

It is what lets a consumer report a fitted break-point term coefficient
by coefficient, and what
[`term_penalties()`](https://statmodels7.github.io/modelterms7/reference/term_penalties.md)
uses to place a sub-term's penalty on exactly the coordinates it covers.
Those coordinates are **named** as a subset of the term's own parameters
and never selected with a map: a separable penalty under a selection map
is the generalized-lasso problem, which has no proximal operator.

## See also

[`term_components()`](https://statmodels7.github.io/modelterms7/reference/term_components.md)
for the contract,
[`term_penalties()`](https://statmodels7.github.io/modelterms7/reference/term_penalties.md)
for the penalties the sub-terms bring.

## Examples

``` r
set.seed(1)
d <- data.frame(x = sort(runif(120, 0, 10)), id = factor(rep(1:4, each = 30)))
d$y <- 1 + 0.5 * d$x + 2 * pmax(d$x - 6, 0) + rnorm(120, sd = 0.4)

# A break-point per subject: psi1 owns four columns, the others one each.
b <- term_build(seg(x, psi ~ id), d)
term_coef_names(b)
#> [1] "seg.beta"             "seg.gamma1"           "seg.psi1.(Intercept)"
#> [4] "seg.psi1.id2"         "seg.psi1.id3"         "seg.psi1.id4"        
lapply(term_components(b), function(z) z$index)
#> $beta
#> [1] 1
#> 
#> $gamma1
#> [1] 2
#> 
#> $psi1
#> [1] 3 4 5 6
#> 
```
