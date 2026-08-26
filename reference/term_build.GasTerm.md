# Build a Score-Driven Term

Resolves the grouping, the ordering and every subformula against the
data, builds the sub-terms of the developments and records all of it in
the blueprint. The recursion itself runs later, in
[`term_filter()`](https://statmodels7.github.io/modelterms7/reference/term_filter.md).

## Arguments

- term:

  A
  [`GasTerm()`](https://statmodels7.github.io/modelterms7/reference/GasTerm.md).

- data:

  A data frame carrying whatever `by`, `time` and the subformulas name.

- ...:

  Unused.

## Value

The term with `blueprint` filled.
[`term_is_built()`](https://statmodels7.github.io/modelterms7/reference/term_is_built.md)
reads that property on this branch, so it is `TRUE` for the result.

## Details

`by` and `time` are evaluated in the data. Each must give one
non-missing value per row, and the blueprint records the row indices of
each group in time order, and the filter iterates over exactly that.

Each subformula's right-hand side goes through
[`interpret_formula()`](https://statmodels7.github.io/modelterms7/reference/interpret_formula.md)
and its terms are built, so their blueprints are recorded and a
prediction reapplies them. A structural sub-term, and one whose own
block moves with its coefficients, are rejected: a parameter's
development must be a fixed design. The penalties those sub-terms carry
become the term's own
[`term_penalties()`](https://statmodels7.github.io/modelterms7/reference/term_penalties.md)
entries, keyed `parameter::subterm`.

Without `time` the rows are taken as they come. That is a real choice:
the recursion is about order, so a data frame that is not already sorted
gives a different model.

## See also

[`gas()`](https://statmodels7.github.io/modelterms7/reference/gas.md),
[`term_filter()`](https://statmodels7.github.io/modelterms7/reference/term_filter.md),
[`term_penalties()`](https://statmodels7.github.io/modelterms7/reference/term_penalties.md).

## Examples

``` r
set.seed(1)
dd <- data.frame(t = 1:60, y = rnorm(60), id = rep(1:3, each = 20),
                 g = factor(rep(c("u", "v"), 30)))

# One series, then a panel of three.
lengths(term_build(gas(p = 1, q = 1, time = t), dd)@blueprint$order)
#>  1 
#> 60 
lengths(term_build(gas(p = 1, q = 1, by = id, time = t), dd)@blueprint$order)
#>  1  2  3 
#> 20 20 20 

# A development expands the parameter and brings its penalty.
gb <- term_build(gas(p = 1, q = 1, omega ~ ridge(~ g), time = t), dd)
term_params(gb)
#> [1] "omega.(Intercept)" "omega.ridge.gu"    "omega.ridge.gv"   
#> [4] "alpha1"            "pacf1"            
vapply(term_penalties(gb), function(e) e$name, character(1))
#> [1] "omega::ridge(~g)"
```
