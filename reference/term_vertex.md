# Whether a Shape-Constrained Term Is at the Vertex of Its Shape

Returns the names of the coefficients of the levels of a term whose
constrained coordinates are all at their bound, and an empty vector when
there are none.

## Usage

``` r
term_vertex(term, coef = NULL, ...)
```

## Arguments

- term:

  A built term.

- coef:

  The coefficients, one per coefficient of the term, or `NULL` for the
  coefficients that the term carries.

- ...:

  Passed to methods.

## Value

A character vector of coefficient names, with the attribute `complete`,
`TRUE` when every level that carries a shape is at the vertex, and for
[`scte()`](https://statmodels7.github.io/modelterms7/reference/scte.md)
the attribute `part`.

## Details

Each coordinate that a shape constrains is \\b_j = l_j +
g^{-1}(\gamma_j)\\, with \\l_j\\ its bound (zero, or a minimum slope).
When the data go against the shape every \\b_j\\ goes to its bound and
every \\\gamma_j\\ to \\-\infty\\: the curve is the one at the bound
(constant, linear for a plain curvature, or of slope `min_slope`) and
the columns of the block on those coordinates vanish. A level is at the
vertex when the derivative \\(g^{-1})'(\gamma_j)\\ of every constrained
coordinate is below \\10^{-10}\\ times its value at the default start:
the inner fit stops somewhere on the way to \\-\infty\\ (measured
between -67 and -320), and a curve whose every constrained slope is that
small is the curve at the bound to that accuracy. The base method
returns an empty vector, which is right for a term without a shape.

Under two monotone shapes,
[`scte()`](https://statmodels7.github.io/modelterms7/reference/scte.md)
has a part that can reach its bound on its own: the elements constrained
in both directions, which carry the interaction. When the term is not at
its vertex and all of them are at their bound, the method returns their
names with the attribute `part` set to `"interaction"`; the surface is
then additive on the scale of the link. A vertex of the whole term has
`part` set to `"term"`.

## Examples

``` r
set.seed(1)
dd <- data.frame(x = runif(50))
tm <- term_build(scs(x, shape_increasing()), dd)
term_vertex(tm, coef = rep(-300, length(term_coef_names(tm))))
#> [1] "scs(x).g1" "scs(x).g2" "scs(x).g3" "scs(x).g4" "scs(x).g5" "scs(x).g6"
#> [7] "scs(x).g7" "scs(x).g8" "scs(x).g9"
#> attr(,"complete")
#> [1] TRUE
term_vertex(tm, coef = rep(0, length(term_coef_names(tm))))
#> character(0)
#> attr(,"complete")
#> [1] FALSE
```
