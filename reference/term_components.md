# How a Term's Columns Divide Among Its Own Parameters

Returns one entry per parameter in which a term is written, giving the
columns that belong to that parameter and the sub-terms that develop it.
The method returns an empty list for a term whose columns are one block
with one meaning.

## Usage

``` r
term_components(term, ...)
```

## Arguments

- term:

  A built term. An unbuilt one returns an empty list.

- ...:

  Passed to methods. The shipped methods do not use it.

## Value

A named list with one element per own parameter, named by that
parameter. Each element is a list with

- `name`:

  the parameter's name, the same as the element's.

- `index`:

  its columns in the term's block, or its positions in
  [`term_params()`](https://statmodels7.github.io/modelterms7/reference/term_params.md)
  for a structural term.

- `subs`:

  the sub-terms developing it, an empty list for a parameter that has no
  sub-terms.

- `sub_index`:

  one integer vector per sub-term, splitting `index` among them; empty
  where `subs` is.

The result is an empty list for a term whose columns do not divide.

## Why the term supplies the division

A term may be written in parameters of its own:
[`nl()`](https://statmodels7.github.io/modelterms7/reference/nl.md) in
the parameters of \\f\\,
[`seg()`](https://statmodels7.github.io/modelterms7/reference/seg.md) in
a slope, a change and a break-point. Any of them may be developed over
covariates, and the block then carries several groups of columns meaning
different things. Only the term records the group to which each column
belongs. A coefficient name is built for a reader, so a consumer reads
the division from this function and does not recover it by parsing
names.

## A parameter may have several sub-terms, of different kinds

`seg(x, psi ~ random(~ 1 | id))` develops the break-point with an
unpenalized intercept and a random block, so that component's `subs` has
two entries and only the second carries a penalty. A consumer reporting
a component therefore reports a sequence of sub-terms.

`sub_index` splits `index` among those sub-terms, in the order in which
the block binds them, which is
[`component_sub_index()`](https://statmodels7.github.io/modelterms7/reference/component_sub_index.md)
applied to their coefficient counts.

## Structural terms

A structural term contributes no design columns, and `index` gives
positions in
[`term_params()`](https://statmodels7.github.io/modelterms7/reference/term_params.md),
which is the vector that indexes its state, its readable quantities and
its variance matrix. In both branches `index` refers to the term's own
coefficients.

The base method returns an empty list, which is the answer for
[`linpar()`](https://statmodels7.github.io/modelterms7/reference/linpar.md),
[`s()`](https://statmodels7.github.io/modelterms7/reference/s.md),
[`random()`](https://statmodels7.github.io/modelterms7/reference/random.md)
and the five penalized constructors, whose columns are one block with
one meaning.

## See also

[`component_sub_index()`](https://statmodels7.github.io/modelterms7/reference/component_sub_index.md)
for the split,
[`term_penalties()`](https://statmodels7.github.io/modelterms7/reference/term_penalties.md)
for the penalties that those sub-terms bring,
[`term_coef_names()`](https://statmodels7.github.io/modelterms7/reference/term_coef_names.md)
for the names of the columns that are divided.

## Examples

``` r
set.seed(1)
dd <- data.frame(x = seq(0.2, 3, length.out = 20),
                 g = factor(rep(c("a", "b"), 10)))
dd$y <- 2 * exp(-1.3 * dd$x)

# The nl() term has two parameters, and the first is developed over a
# factor, which gives two columns for `a` and one for `r`.
b <- term_build(nl(~ a * exp(-r * x), a ~ 0 + g, start = list(r = 1.3)), dd)
lapply(term_components(b), function(z) z$index)
#> $a
#> [1] 1 2
#> 
#> $r
#> [1] 3
#> 
term_coef_names(b)
#> [1] "nl.a.ga" "nl.a.gb" "nl.r"   

# This lists the sub-terms of the developed parameter and splits its
# columns among them.
term_components(b)$a$sub_index
#> [[1]]
#> [1] 1 2
#> 

# A term whose columns are one block returns an empty list.
term_components(term_build(linpar(~ x), dd))
#> list()
```
