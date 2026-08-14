# Starting Positions for a Break-Point Term

Chooses the starting positions of a
[`seg`](https://statmodels7.github.io/modelterms7/reference/seg.md),
[`jump`](https://statmodels7.github.io/modelterms7/reference/jump.md) or
[`jseg`](https://statmodels7.github.io/modelterms7/reference/jseg.md)
term by scoring an equally spaced grid on the least-squares profile of
the term's own columns, and returns the specification with `psi` set to
the best combination found.

## Usage

``` r
seg_start(spec, data, y, k = 10)
```

## Arguments

- spec:

  An unbuilt
  [`SegTerm`](https://statmodels7.github.io/modelterms7/reference/SegTerm.md).

- data:

  A data frame in which the covariate is evaluated.

- y:

  The response, one value per row of `data`.

- k:

  The number of grid points. Defaults to 10.

## Value

The specification, with `psi` set.

## Details

fasola2018 recommend fixing the starting value by evaluating the
objective on a small grid spanned over the range of the covariate rather
than at a single conventional point, and the recommendation matters more
than it sounds: the objective has local optima in the break-point, and
the iteration converges from within a basin around the position it
starts at. Measured on a joint jump and change of slope in 500
observations, over eight samples, the fraction of runs recovering the
break-point is 0 to 0.5 depending on where a single start is placed and
1 from the grid.

Writing \\X(\psi)\\ for the design the term produces at a candidate
position, the position chosen is

\$\$\hat\psi = \arg\min\_{\psi \in \mathcal{G}} \bigl\lVert y -
X(\psi)\\ \widehat{\beta}(\psi) \bigr\rVert^{2}, \qquad
\widehat{\beta}(\psi) = \arg\min\_{\beta} \lVert y - X(\psi)\beta
\rVert^{2},\$\$

over an equally spaced grid \\\mathcal{G}\\ of `k` points in the range
of the covariate. The inner minimization is a linear fit, so the whole
rule costs `k` of them.

The grid is scored on the residual sum of squares of an intercept, the
term's columns at each candidate position and, where the term carries
one, the linear effect. That is the exact profile for a gaussian
response and an adequate starting rule for any other, the quantity being
used to place a starting value and not to fit. The positions found seed
a development as well, each starting vector projecting the position onto
the sub-design. With several break-points every increasing combination
of grid points is scored, so `k` should be kept small.

## References

Fasola, S., Muggeo, V. M. R. and Kuchenhoff, H. (2018). A heuristic,
iterative algorithm for change-point detection in abrupt change models.
*Computational Statistics*, 33, 997–1015.

## See also

[`seg`](https://statmodels7.github.io/modelterms7/reference/seg.md),
[`jump`](https://statmodels7.github.io/modelterms7/reference/jump.md),
[`jseg`](https://statmodels7.github.io/modelterms7/reference/jseg.md)

## Examples

``` r
set.seed(1)
dd <- data.frame(x = sort(runif(200, 0, 10)))
dd$y <- 0.3 * dd$x + 2 * (dd$x > 6.5) + rnorm(200, sd = 0.3)
seg_start(jump(x), dd, dd$y)@spec$psi
#> [1] 6.462078
```
