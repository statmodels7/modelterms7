# The Hyperparameters a Term Holds

The values the caller fixed in the constructor, by hyperparameter. Those
a term does not name are estimated by whatever criterion the fit runs.

## Usage

``` r
term_hyper(term, ...)
```

## Arguments

- term:

  A term, built or not.

- ...:

  Passed to methods.

## Value

A named list, one entry per penalty of the term, each a named list of
held values. Empty where the term holds nothing.

## Details

A term carrying several penalties answers per penalty, under the same
names
[`term_penalties`](https://statmodels7.github.io/modelterms7/reference/term_penalties.md)
gives its entries, and every entry of that enumeration carries its own
held values in the field `fixed` – so a term that copies the entries of
its sub-terms, which is what a structural term with subformulas does,
propagates them without knowing they exist.

## See also

[`term_penalties`](https://statmodels7.github.io/modelterms7/reference/term_penalties.md),
[`ridge`](https://statmodels7.github.io/modelterms7/reference/ridge.md)

## Examples

``` r
term_hyper(lasso(~x, lambda = 3))
#> [[1]]
#> [[1]]$lambda
#> [1] 3
#> 
#> 
term_hyper(lasso(~x))
#> list()
```
