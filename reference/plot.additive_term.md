# Plot an Additive Term's Coefficients

[`autoplot()`](https://statmodels7.github.io/modelterms7/reference/autoplot.md)
returns a graph of the coefficients of a built additive term as a stem
plot, and [`plot()`](https://rdrr.io/r/graphics/plot.default.html)
prints that graph: one point per coefficient against its position in the
block, a stem down to zero, a grey line at zero, and the coefficient
names along the horizontal axis. It shows which coefficients a penalty
has driven to zero and how large the survivors are, so the effect of a
lasso or a smooth is visible at a glance.

The coefficients are supplied by the caller. A term holds a design block
and no fit, so there is nothing to display without them, and
`coef = NULL` is rejected.

## Usage

``` r
# S3 method for class 'additive_term'
plot(x, coef = NULL, ...)

autoplot_additive_term(object, coef = NULL, ...)
```

## Arguments

- x, object:

  A built additive term. An unbuilt term signals the error
  `"the term has not been built; call term_build(term, data) first."`.

- coef:

  The coefficients to draw, a numeric vector of length `term_npar(x)`.
  Required: `NULL` signals the error
  `"'coef' is required: a term is displayed at fitted coefficients."`,
  and a vector of any other length signals an error that gives the
  required length.

- ...:

  Unused; an argument passed here signals an error.

## Value

[`autoplot()`](https://statmodels7.github.io/modelterms7/reference/autoplot.md)
returns a ggplot2 object;
[`plot()`](https://rdrr.io/r/graphics/plot.default.html) prints the same
object and returns it invisibly.

## Details

The title of the graph is the term's `label` when it has one, so
`"ridge"` or `"s(x)"`, and the class name otherwise.

The horizontal axis is the coefficient's position in the block, so the
graph shows the block and not the covariate. For
[`s()`](https://statmodels7.github.io/modelterms7/reference/s.md), whose
block is a Demmler-Reinsch reparametrization, the columns are ordered
from the smoothest to the most wiggly, and a penalized fit shows a
decaying profile.

The graph is drawn with ggplot2, which is suggested; the methods signal
an error that names it when it is not installed. The graph has one
element, the coefficients, which the vertical axis names, and so it has
no legend.

## See also

[`term_coef_names()`](https://statmodels7.github.io/modelterms7/reference/term_coef_names.md)
for the axis labels,
[`edf()`](https://statmodels7.github.io/modelterms7/reference/edf.md)
for what the same coefficients cost in degrees of freedom.

## Examples

``` r
set.seed(5)
d <- data.frame(x1 = rnorm(30), x2 = rnorm(30), x3 = rnorm(30))
b <- term_build(lasso(~ x1 + x2 + x3), d)

# The graph shows two survivors and one coefficient at exactly zero.
plot(b, coef = c(0.8, 0, -0.35))


# A smooth's block runs from least to most wiggly.
d2 <- data.frame(x = seq(0, 1, length.out = 60))
bs <- term_build(s(x, basis7::bspline_smooth(k = 8)), d2)
plot(bs, coef = c(1.2, 0.9, -0.4, 0.2, -0.1, 0.05, -0.02))


# Coefficients are required, and must fit the block.
try(plot(b))
#> Error : 'coef' is required: a term is displayed at fitted coefficients.
try(plot(b, coef = 1))
#> Error : 'coef' must have length 3.
```
