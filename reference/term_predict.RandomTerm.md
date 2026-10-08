# A Random-Effect Block at New Rows

Rebuilds the within-group design at `newdata` and interacts it with the
group indicators **of the levels recorded at build time**, so the block
has the same columns in the same order however few levels the new rows
happen to use. A level the term never saw is refused.

## Arguments

- term:

  A built
  [`RandomTerm()`](https://statmodels7.github.io/modelterms7/reference/RandomTerm.md).
  An unbuilt one throws
  `"the term has not been built; call term_build(term, data) first."`.

- newdata:

  A data frame carrying the grouping variable and the within-group
  covariates. Its grouping factor need carry only the levels its own
  rows use.

- unseen:

  What a row of a level the term never saw gets: `"error"` (the default)
  signals the error described above, `"zero"` gives it a row of zeros,
  which is the prediction with that group's effect at zero. It is what
  `statmodels7`'s `predict(random = "zero")` and
  `predict(random = "marginal")` ask for.

- ...:

  Unused.

## Value

A `dgCMatrix` of `nrow(newdata)` rows and
[`term_npar()`](https://statmodels7.github.io/modelterms7/reference/term_npar.md)
columns, with the term's coefficient names as column names.

## Details

The refusal is the right answer rather than a limitation: a coefficient
was never fitted for an unseen group, so there is nothing to predict
with. The message names the level. Predicting a new group's response
means predicting at the population value, which is the model without
this term's contribution.

The block comes back sparse, as the fitted one is.

## See also

[`term_predict()`](https://statmodels7.github.io/modelterms7/reference/term_predict.md)
for the generic,
[`term_build.RandomTerm()`](https://statmodels7.github.io/modelterms7/reference/term_build.RandomTerm.md)
for what recorded the levels.

## Examples

``` r
set.seed(1)
dd <- data.frame(x = rnorm(9), g = factor(rep(c("a", "b", "c"), 3)))
b <- term_build(random(~ 1 | g), dd)

# A subset using two levels still gets all three columns.
nd <- droplevels(dd[dd$g != "c", ])
c(levels_here = nlevels(nd$g), cols = ncol(term_predict(b, nd)))
#> levels_here        cols 
#>           2           3 

# On the fitting data it returns the block itself.
all.equal(term_predict(b, dd), term_matrix(b))
#> [1] TRUE

# A level the fit never saw has no coefficient, so it is refused.
bad <- dd
levels(bad$g) <- c("a", "b", "zz")
try(term_predict(b, bad))
#> Error : grouping level 'zz' was not present at build time.

# Unless it is asked to give such a row no effect of its own.
as.matrix(term_predict(b, bad, unseen = "zero"))[bad$g == "zz", ]
#>      random.a random.b random.c
#> [1,]        0        0        0
#> [2,]        0        0        0
#> [3,]        0        0        0
```
