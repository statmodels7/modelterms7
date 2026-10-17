# The Grid Sizes of a Term

Reports the number of values that a path visits for each hyperparameter
of the term that is swept, one entry per penalty. A hyperparameter that
the term holds has no grid and is not reported. A swept hyperparameter
for which the term records no size is covered at the default size of the
fitting layer.

## Usage

``` r
term_grid(term, ...)
```

## Arguments

- term:

  A term, built or not.

- ...:

  Passed to methods. The shipped methods do not use it.

## Value

A named list, one entry per penalty of the term, each a named list of
single whole numbers keyed by hyperparameter. The list is empty where
the term records no grid size for any hyperparameter.

## Details

How finely a hyperparameter is swept belongs to the term because the
term is where the penalty is named: a block of four columns and one of
four hundred need different grids, and an outer criterion applies to
every hyperparameter of the model alike, so a size set there could not
distinguish between them.

Only a hyperparameter with a path has a grid, and a held hyperparameter
has none.
[`ridge()`](https://statmodels7.github.io/modelterms7/reference/ridge.md),
[`s()`](https://statmodels7.github.io/modelterms7/reference/s.md) and
[`te()`](https://statmodels7.github.io/modelterms7/reference/te.md) with
their default penalty, and
[`random()`](https://statmodels7.github.io/modelterms7/reference/random.md)
with its default Gaussian prior, report nothing, because a marginal
criterion estimates their hyperparameters at the mode. A smoother whose
penalty factory has a kink is swept by a path but reports nothing
either, and the path uses its default size.
[`lasso()`](https://statmodels7.github.io/modelterms7/reference/lasso.md)
reports its `n_lambda`,
[`enet()`](https://statmodels7.github.io/modelterms7/reference/enet.md)
reports `n_lambda` and `n_alpha`, and
[`scad()`](https://statmodels7.github.io/modelterms7/reference/scad.md)
and
[`mcp()`](https://statmodels7.github.io/modelterms7/reference/mcp.md)
report `n_lambda` together with `n_a` or `n_gamma` where the shape is
estimated, since the shape is held by default. These constructors write
their sizes into the term, so `lasso(~ x)` reports `lambda = 25` where
an unset argument would give an empty list. A hyperparameter written out
as a vector has no grid to build and is not reported here; the path
visits the written-out values, which
[`term_values()`](https://statmodels7.github.io/modelterms7/reference/term_values.md)
reports.

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

# The two axes have two sizes, and the kink axis is swept more finely.
term_grid(enet(~ x))
#> [[1]]
#> [[1]]$lambda
#> [1] 25
#> 
#> [[1]]$alpha
#> [1] 5
#> 
#> 

# A held hyperparameter has no grid, and the held shape of SCAD has none.
term_grid(enet(~ x, alpha = 0.5))
#> [[1]]
#> [[1]]$lambda
#> [1] 25
#> 
#> 
term_grid(scad(~ x))
#> [[1]]
#> [[1]]$lambda
#> [1] 25
#> 
#> 
term_grid(scad(~ x, a = NULL))
#> [[1]]
#> [[1]]$lambda
#> [1] 25
#> 
#> [[1]]$a
#> [1] 5
#> 
#> 

# A penalty with no kink has no path, so there is nothing to count.
term_grid(ridge(~ x))
#> list()
term_grid(s(x, basis7::bspline_smooth(k = 5)))
#> list()
```
