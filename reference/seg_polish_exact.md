# Polish a Break-Point Term's Positions Over Every Interval

The exact minimum of the least-squares profile of a sharp
[`seg()`](https://statmodels7.github.io/modelterms7/reference/seg.md),
[`jump()`](https://statmodels7.github.io/modelterms7/reference/jump.md)
or
[`jseg()`](https://statmodels7.github.io/modelterms7/reference/jseg.md)
term, one break-point at a time with the others held. The profile of a
change of level is constant between consecutive observations of the
covariate, and the profile of a change of slope has one stationary point
inside each such interval, so its minimum over one break-point is found
interval by interval, and this function does that for every interval
inside the confinement limits.

## Usage

``` r
seg_polish_exact(term, y, sweeps = 10, weights = NULL)
```

## Arguments

- term:

  A built break-point term (see
  [`term_build()`](https://statmodels7.github.io/modelterms7/reference/term_build.md)).

- y:

  A numeric vector, one value per observation of the build data: the
  response, net of whatever the caller wants held.

- sweeps:

  At most how many passes over the positions. Defaults to 10; the
  descent usually stops moving in two or three.

- weights:

  Optional non-negative weights, one per observation; the profile is
  then weighted least squares.

## Value

The term at the polished positions (see
[`seg_relocate()`](https://statmodels7.github.io/modelterms7/reference/seg_relocate.md)).

## Details

Write \\F\\ for the columns that do not move with the break-point being
polished (an intercept, the covariate where the term carries a linear
part, and the columns of the other break-points), and \\a(q)\\ for the
columns that do: the indicator \\1(x \> q)\\, and for a
[`jseg()`](https://statmodels7.github.io/modelterms7/reference/jseg.md)
also the truncated line \\(x - q)1(x \> q)\\. The residual sum of
squares at \\q\\ is that of the response on \\F\\ less the part \\a(q)\\
explains of what \\F\\ leaves, and every cross product \\a(q)\\ enters
is a sum over the observations above \\q\\ of a fixed quantity, or of
one fixed quantity minus \\q\\ times another. Sorted once, those sums
are cumulative sums, so all the intervals cost \\O(np^2)\\ together
rather than one linear fit each. The sweeps over the break-points repeat
until none moves. Each break-point of a
[`jump()`](https://statmodels7.github.io/modelterms7/reference/jump.md)
or a
[`jseg()`](https://statmodels7.github.io/modelterms7/reference/jseg.md)
is placed at the midpoint of its interval.

For a
[`seg()`](https://statmodels7.github.io/modelterms7/reference/seg.md)
the column that moves is \\(x - q)\_+\\. Inside the interval between the
consecutive values \\u_k \< u\_{k+1}\\ it is \\t - q s\\, with \\s = 1(x
\> u_k)\\ and \\t = x s\\, so the part of the residual sum of squares it
removes is \\(c_t - q c_s)^2 / (a\_{tt} - 2 q a\_{ts} + q^2 a\_{ss})\\,
the \\c\\ and \\a\\ being the cross products of \\t\\, \\s\\ and the
response after the fixed columns are projected out. Its one stationary
point is \\q^\* = -b/\gamma\\ of the least-squares fit on \\t\\ and
\\s\\, so the minimum over the interval is at \\q^\*\\ when it falls
inside, or at an end, which is an observed value.

[`seg_polish()`](https://statmodels7.github.io/modelterms7/reference/seg_polish.md)
sweeps a grid of `k` points instead, which reaches a neighbourhood of
the minimum and not the interval: measured on 16 samples of 200
observations, a fit polished by the grid reached the interval of the
global minimum in 11.

The profile is least squares of `y` on the term's own columns plus an
intercept, exact for a gaussian response; a fitting layer for any other
family accepts the positions only where its own objective improves.

A held term whose break-point is developed, \\\psi_i = w_i'p\\ with
\\w_i\\ the row of the sub-design, is polished in the coefficients
\\p\\. Along a line \\p + tv\\ the position of observation \\i\\ crosses
\\x_i\\ at \\t_i = (x_i - \psi_i)/(w_i'v)\\, so the profile is evaluated
once between each pair of consecutive crossings, and for a
[`jseg()`](https://statmodels7.github.io/modelterms7/reference/jseg.md)
whose change of slope is not developed like its position it is then
minimized inside the best interval, where it varies smoothly. Where the
sub-design partitions the observations into groups (`psi ~ g`,
`by = ~ 0 + g`) the lines move one group's position at a time, which is
a search over every interval of that group, and the sweep also starts
from each group's own minimum, the minimum of the profile on that
group's observations alone; the better of the two results is kept. With
a continuous covariate in the sub-design the lines are the coordinate
directions and eight fixed directions in each plane of two coordinates,
so the result is the best point found along those lines and not a global
minimum. A developed term is polished only when held (see
[`seg_hold()`](https://statmodels7.github.io/modelterms7/reference/seg_hold.md)).

## See also

[`seg_polish()`](https://statmodels7.github.io/modelterms7/reference/seg_polish.md),
[`seg_hold()`](https://statmodels7.github.io/modelterms7/reference/seg_hold.md),
[`seg_profile_rss()`](https://statmodels7.github.io/modelterms7/reference/seg_profile_rss.md).

## Examples

``` r
set.seed(1)
dd <- data.frame(x = sort(runif(300, 0, 10)))
dd$y <- 2 * (dd$x > 3) - 1.5 * (dd$x > 7) + rnorm(300, sd = 0.3)
b <- term_build(jump(x, npsi = 2, psi = c(2, 5)), dd)
seg_psi(seg_polish_exact(b, dd$y))
#> [1] 3.025479 6.975970
```
