# Which of a Term's Parameters Acts as an Intercept

Returns the name of the parameter, if any, that shifts the term's
contribution by a constant and is therefore the same direction as an
intercept in the equation in which the term sits.

## Usage

``` r
term_level_param(term, ...)
```

## Arguments

- term:

  A term.

- ...:

  Passed to methods.

## Value

A character vector: `character(0)` for a term that shifts nothing by a
constant, `"level1"` for a regime, `"omega"` for an undeveloped
[`gas()`](https://statmodels7.github.io/modelterms7/reference/gas.md)
level, and otherwise the names of the level's coordinates whose design
column is constant.

## Details

It exists so that a fitting layer can resolve the confounding instead of
rejecting the model. A score-driven level \\\omega\\ and a regime's
first level both add a constant to their equation's predictor: with an
intercept there too, adding \\c\\ to one and subtracting the matching
amount from the other leaves every predictor unchanged, and the
likelihood is flat along that direction. Which of the two is dropped is
determined in the layer, because only the layer has the full list of
terms in the equation.

The base method returns `character(0)`, which means that the term does
not shift its contribution by a constant and the question does not
arise. The method for
[`regime()`](https://statmodels7.github.io/modelterms7/reference/regime.md)
returns `"level1"`. The method for
[`gas()`](https://statmodels7.github.io/modelterms7/reference/gas.md)
returns `"omega"` when the level is not developed and otherwise every
coordinate of the development whose design column is constant, which can
be more than one.

## See also

[`term_level_design()`](https://statmodels7.github.io/modelterms7/reference/term_level_design.md)
for the same question with the level developed,
[`term_params()`](https://statmodels7.github.io/modelterms7/reference/term_params.md)
for the names from which the answer is taken.

## Examples

``` r
term_level_param(gas(p = 1, q = 1))
#> [1] "omega"
term_level_param(regime(k = 2))
#> [1] "level1"
term_level_param(linpar(~x))
#> character(0)
```
