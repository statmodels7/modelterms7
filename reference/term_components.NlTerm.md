# The Columns Each Parameter of a Nonlinear Term Owns

One entry per parameter of \\f\\, giving the columns of the block that
parameter owns and the sub-terms developing it. A parameter with no
subformula owns one column; a developed one owns as many as its own
design has.

## Arguments

- term:

  A built
  [`NlTerm()`](https://statmodels7.github.io/modelterms7/reference/NlTerm.md).
  An unbuilt one gives an empty list.

- ...:

  Unused.

## Value

A named list, one entry per parameter of \\f\\ and named by it, each
with `name`, `index`, `subs` and `sub_index` as
[`term_components()`](https://statmodels7.github.io/modelterms7/reference/term_components.md)
describes.

## Details

A consumer reads it to report a fitted nonlinear term parameter by
parameter, and
[`edf()`](https://statmodels7.github.io/modelterms7/reference/edf.md)
counts over it. A coefficient name is built for a reader, so the
division cannot be recovered by parsing one back; the term has to say
it.

## See also

[`term_components()`](https://statmodels7.github.io/modelterms7/reference/term_components.md)
for the contract,
[`term_penalties()`](https://statmodels7.github.io/modelterms7/reference/term_penalties.md)
for the penalties those sub-terms bring.

## Examples

``` r
set.seed(1)
dd <- data.frame(x = seq(0.2, 3, length.out = 20),
                 g = factor(rep(c("a", "b"), 10)))
dd$y <- 2 * exp(-1.3 * dd$x)

# `a` developed over a two-level factor, `r` scalar: two columns and one.
b <- term_build(nl(~ a * exp(-r * x), a ~ 0 + g,
                   start = list(a = 1, r = 1.3)), dd)
term_coef_names(b)
#> [1] "nl.a.ga" "nl.a.gb" "nl.r"   
lapply(term_components(b), function(z) z$index)
#> $a
#> [1] 1 2
#> 
#> $r
#> [1] 3
#> 
```
