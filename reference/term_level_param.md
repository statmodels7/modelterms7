# Which of a Term's Parameters Acts as an Intercept

The name of the parameter, if any, that shifts the term's contribution
by a constant, and is therefore the same direction as an intercept in
the equation the term sits in.

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

A single string, or `character(0)`.

## Details

It exists so that a fitting layer can resolve the confounding rather
than refuse the model. A score-driven level \\\omega\\ and a regime's
first level both add a constant to their equation's predictor: with an
intercept there too, adding \\c\\ to one and subtracting the matching
amount from the other leaves every predictor unchanged, and the
likelihood is flat along that direction. Which of the two is dropped is
the layer's decision and not the term's, since only the layer knows what
else the equation carries.

The base method returns `character(0)`, meaning the term shifts nothing
by a constant and no such question arises.

## See also

[`term_params`](https://statmodels7.github.io/modelterms7/reference/term_params.md)

## Examples

``` r
term_level_param(gas(p = 1, q = 1))
#> [1] "omega"
term_level_param(linpar(~x))
#> character(0)
```
