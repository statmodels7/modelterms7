# The Grid a Term Asks For

Reports how many values a path visits for each of the term's
hyperparameters, one entry per penalty. A hyperparameter the term names
nothing for is swept at the fitting layer's own default.

## Usage

``` r
term_grid(term, ...)
```

## Arguments

- term:

  A term, built or not.

- ...:

  Passed to methods. No shipped method reads anything here.

## Value

A named list, one entry per penalty of the term, each a named list of
single whole numbers keyed by hyperparameter. Empty where the term names
no grid at all.

## Details

How finely a hyperparameter is swept belongs to the term because the
term is where the penalty is named: a block of four columns and one of
four hundred want different grids, and an outer criterion, which is put
to every hyperparameter of the model, does not know which it is looking
at.

Only a hyperparameter with a **path** has a grid.
[`ridge()`](https://statmodels7.github.io/modelterms7/reference/ridge.md),
[`s()`](https://statmodels7.github.io/modelterms7/reference/s.md),
[`te()`](https://statmodels7.github.io/modelterms7/reference/te.md) and
[`random()`](https://statmodels7.github.io/modelterms7/reference/random.md)
report nothing, their hyperparameters being estimated at the mode by a
marginal criterion;
[`lasso()`](https://statmodels7.github.io/modelterms7/reference/lasso.md)
reports its `n_lambda`, and
[`enet()`](https://statmodels7.github.io/modelterms7/reference/enet.md),
[`scad()`](https://statmodels7.github.io/modelterms7/reference/scad.md)
and
[`mcp()`](https://statmodels7.github.io/modelterms7/reference/mcp.md)
report both of theirs. Those constructors write their defaults into the
term, so `lasso(~ x)` reports `lambda = 25` where an unset argument
would give an empty list.

The keys are
[`term_penalties()`](https://statmodels7.github.io/modelterms7/reference/term_penalties.md)'s
entry names, `""` for a penalty over the whole term.

## See also

[`term_path_min()`](https://statmodels7.github.io/modelterms7/reference/term_path_min.md)
for how far down the path reaches,
[`term_search()`](https://statmodels7.github.io/modelterms7/reference/term_search.md)
for how several of them are combined,
[`term_values()`](https://statmodels7.github.io/modelterms7/reference/term_values.md)
for a grid written out instead of counted.

## Examples

``` r
# A path term reports its grid, at the default and when set.
term_grid(lasso(~ x))
#> [[1]]
#> [[1]]$lambda
#> [1] 25
#> 
#> 
term_grid(lasso(~ x, n_lambda = 50))
#> [[1]]
#> [[1]]$lambda
#> [1] 50
#> 
#> 

# Two axes, and their two sizes: the kink axis is swept far more finely.
term_grid(enet(~ x))
#> [[1]]
#> [[1]]$lambda
#> [1] 25
#> 
#> [[1]]$alpha
#> [1] 5
#> 
#> 

# A penalty with no kink has no path, so nothing to count.
term_grid(ridge(~ x))
#> list()
term_grid(s(x, k = 5))
#> list()
```
