# How a Term's Columns Divide Among Its Own Parameters

One entry per parameter the term is written in, saying which columns of
its block belong to that parameter and which sub-terms develop it, or an
empty list for a term whose columns answer to nothing above them.

## Usage

``` r
term_components(term, ...)
```

## Arguments

- term:

  A built term.

- ...:

  Passed to methods.

## Value

A list, one entry per own parameter, each with `name` (the parameter),
`index` (its columns in the term's block), `subs` (the sub-terms
developing it, empty where there are none) and `sub_index` (the columns
belonging to each of those sub-terms). Empty for a term whose columns do
not divide.

## Details

A term written in parameters of its own –
[`nl`](https://statmodels7.github.io/modelterms7/reference/nl.md) in the
parameters of \\f\\,
[`seg`](https://statmodels7.github.io/modelterms7/reference/seg.md) in a
slope, a change and a break-point – may develop any of them over
covariates, and then its block carries several groups of columns that
mean different things. What divides them is the TERM's answer and cannot
be recovered from the coefficient names: a name is built for a reader
and parsing one back is the shape of mistake this package avoids
everywhere else.

A parameter may be developed by SEVERAL sub-terms at once, and they need
not be of one kind: `seg(x, psi ~ random(~1 | id))` develops the
break-point with an unpenalized intercept AND a random block, so the
component's `subs` has two entries and only the second carries a
penalty. A consumer that reports a component therefore reports a
sequence and not a single kind.

A STRUCTURAL term contributes no design columns, and there `index` gives
positions in
[`term_params`](https://statmodels7.github.io/modelterms7/reference/term_params.md)
instead: the vector its state, its readable quantities and its variance
matrix are all indexed by. In both cases the field names the term's own
coefficients.

The base method returns an empty list, which says that the term's
columns are its own and divide no further. That is the honest answer for
[`linpar`](https://statmodels7.github.io/modelterms7/reference/linpar.md),
[`s`](https://statmodels7.github.io/modelterms7/reference/s.md),
[`random`](https://statmodels7.github.io/modelterms7/reference/random.md)
and the penalized constructors, whose columns are one block with one
meaning.

## See also

[`term_penalties`](https://statmodels7.github.io/modelterms7/reference/term_penalties.md),
[`term_coef_names`](https://statmodels7.github.io/modelterms7/reference/term_coef_names.md)

## Examples

``` r
dd <- data.frame(x = seq(0.2, 3, length.out = 20),
                 g = factor(rep(c("a", "b"), 10)))
dd$y <- 2 * exp(-1.3 * dd$x)
b <- term_build(nl(~ a * exp(-r * x), a ~ 0 + g, start = list(r = 1.3)), dd)
lapply(term_components(b), function(z) z$index)
#> $a
#> [1] 1 2
#> 
#> $r
#> [1] 3
#> 

# a term whose columns are one block answers with nothing
term_components(term_build(linpar(~ x), dd))
#> list()
```
