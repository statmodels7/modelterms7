# How Far Down Its Path a Term Reaches

Reports the fraction of the emptying value the path descends to. The
path runs from the kink that leaves every coefficient of the block at
zero down to that fraction of it, so a smaller number reaches a denser
fit and a larger one stops sooner.

## Usage

``` r
term_path_min(term, ...)
```

## Arguments

- term:

  A term, built or not.

- ...:

  Passed to methods. No shipped method reads anything here.

## Value

A named list, one entry per penalty of the term, each a single number in
\\(0, 1)\\. Empty where the term names none.

## Details

It belongs to the term for the same reason the grid size does: how far
the useful range of a hyperparameter extends is a property of the block,
and a criterion applies to every term of the model at once.

It is **one number per penalty**, one per hyperparameter being
unnecessary: only the path over the size of the kink uses it. A bounded
hyperparameter is swept over its own interval, and a shape that does not
move the kink over a geometric grid above its lower bound; a fraction of
an emptying value means nothing in either.

[`lasso()`](https://statmodels7.github.io/modelterms7/reference/lasso.md),
[`enet()`](https://statmodels7.github.io/modelterms7/reference/enet.md),
[`scad()`](https://statmodels7.github.io/modelterms7/reference/scad.md)
and
[`mcp()`](https://statmodels7.github.io/modelterms7/reference/mcp.md)
write their `min_ratio` default of `1e-4` into the term, so they report
it whether or not the caller set one.
[`ridge()`](https://statmodels7.github.io/modelterms7/reference/ridge.md),
[`s()`](https://statmodels7.github.io/modelterms7/reference/s.md),
[`te()`](https://statmodels7.github.io/modelterms7/reference/te.md) and
[`random()`](https://statmodels7.github.io/modelterms7/reference/random.md)
have no path and report nothing.

The keys are
[`term_penalties()`](https://statmodels7.github.io/modelterms7/reference/term_penalties.md)'s
entry names, `""` for a penalty over the whole term.

## See also

[`term_grid()`](https://statmodels7.github.io/modelterms7/reference/term_grid.md)
for how many values the path visits,
[`term_search()`](https://statmodels7.github.io/modelterms7/reference/term_search.md)
for how several hyperparameters are combined.

## Examples

``` r
# The default, and a deeper path.
term_path_min(lasso(~ x))
#> [[1]]
#> [1] 1e-04
#> 
term_path_min(lasso(~ x, min_ratio = 1e-6))
#> [[1]]
#> [1] 1e-06
#> 

# No path, nothing to report.
term_path_min(ridge(~ x))
#> list()
```
