# Polish a Break-Point Term's Positions Over Every Interval

Minimizes the least-squares profile of a sharp
[`seg()`](https://statmodels7.github.io/modelterms7/reference/seg.md),
[`jump()`](https://statmodels7.github.io/modelterms7/reference/jump.md)
or
[`jseg()`](https://statmodels7.github.io/modelterms7/reference/jseg.md)
term over the positions of its break-points, one break-point at a time
with the others held. The profile of a
[`jump()`](https://statmodels7.github.io/modelterms7/reference/jump.md)
or a
[`jseg()`](https://statmodels7.github.io/modelterms7/reference/jseg.md)
is constant between consecutive distinct values of the covariate, and
the profile of a
[`seg()`](https://statmodels7.github.io/modelterms7/reference/seg.md)
has one stationary point inside each such interval. The minimum over one
break-point is therefore found interval by interval, and this function
does that for every interval inside the confinement limits.

## Usage

``` r
seg_polish_exact(term, y, sweeps = 10, weights = NULL)
```

## Arguments

- term:

  A built break-point term (see
  [`term_build()`](https://statmodels7.github.io/modelterms7/reference/term_build.md)).

- y:

  A numeric vector of finite values, one per observation of the build
  data: the response, net of whatever the caller wants held.

- sweeps:

  At most how many passes over the positions. Defaults to 10. The
  descent stops once a full pass moves no position.

- weights:

  Optional non-negative weights, one per observation; the profile is
  then weighted least squares.

## Value

The term at the polished positions (see
[`seg_relocate()`](https://statmodels7.github.io/modelterms7/reference/seg_relocate.md)),
or the term itself when no interval lies inside the confinement limits.

## Details

Write \\F\\ for the columns that do not move with the break-point being
polished (an intercept, the covariate where the term carries a linear
part, and the columns of the other break-points), and \\a(q)\\ for the
columns that do: the indicator \\1(x \> q)\\, and for a
[`jseg()`](https://statmodels7.github.io/modelterms7/reference/jseg.md)
also the truncated line \\(x - q)1(x \> q)\\. The residual sum of
squares at \\q\\ is that of the response on \\F\\ less the part of the
remainder that \\a(q)\\ explains, and every cross product that involves
\\a(q)\\ is a sum over the observations above \\q\\ of a fixed quantity,
or of one fixed quantity minus \\q\\ times another. Sorted once, those
sums are cumulative sums, so all the intervals cost \\O(np^2)\\
together, where one linear fit per interval would cost more. The sweeps
over the break-points repeat, at most `sweeps` times, until none moves.
Inside an interval between consecutive distinct values of the covariate,
the span of the columns of a
[`jump()`](https://statmodels7.github.io/modelterms7/reference/jump.md)
or a
[`jseg()`](https://statmodels7.github.io/modelterms7/reference/jseg.md)
does not depend on \\q\\, so each break-point of these terms is placed
at the midpoint of its interval.

For a
[`seg()`](https://statmodels7.github.io/modelterms7/reference/seg.md)
the column that moves is \\(x - q)\_+\\. Inside the interval between the
consecutive values \\u_k \< u\_{k+1}\\ it is \\t - q s\\, with \\s = 1(x
\> u_k)\\ and \\t = x s\\, so the part of the residual sum of squares it
removes is \\(c_t - q c_s)^2 / (a\_{tt} - 2 q a\_{ts} + q^2 a\_{ss})\\,
the \\c\\ and \\a\\ being the cross products of \\t\\, \\s\\ and the
response after the fixed columns are projected out. Its one stationary
point is \\q^\* = -b/\gamma\\ of the least-squares fit on \\t\\ and
\\s\\. The minimum over the interval is at \\q^\*\\, when that falls
inside the interval, or at one of the two ends, whichever gives the
smaller residual sum of squares. The ends are consecutive observed
values of the covariate or confinement limits.

[`seg_polish()`](https://statmodels7.github.io/modelterms7/reference/seg_polish.md)
sweeps instead a grid of `k` quantiles of the covariate, so that the
position it reaches is accurate only to the resolution of the grid.

The profile is the least-squares fit of `y` on an intercept and the
term's own columns, which is exact for a gaussian response. For any
other family, the fitting layer uses the positions that the profile
gives only where its own objective improves. The response must be
finite, and missing values are rejected.

A held term whose break-point is developed, \\\psi_i = w_i'p\\ with
\\w_i\\ the row of the sub-design, is polished in the coefficients
\\p\\. Along a line \\p + tv\\ the position of observation \\i\\ crosses
\\x_i\\ at \\t_i = (x_i - \psi_i)/(w_i'v)\\, so the profile is evaluated
once between each pair of consecutive crossings. For a
[`jseg()`](https://statmodels7.github.io/modelterms7/reference/jseg.md)
the profile is then also minimized inside the best interval, where it
varies smoothly. Where the sub-design partitions the observations into
groups (`psi ~ g`, `by = ~ 0 + g`), the lines move the position of one
group at a time, which is a search over every interval of that group.
The sweep is also started from each group's own minimum, which is the
minimum of the profile on the observations of that group alone, and the
better of the two results is kept. Otherwise the lines are the
coordinate directions and, in each plane of two coordinates, six further
fixed directions, so that the result is the best point found along those
lines and not a global minimum. A developed
[`jump()`](https://statmodels7.github.io/modelterms7/reference/jump.md)
or
[`jseg()`](https://statmodels7.github.io/modelterms7/reference/jseg.md)
is polished only when held (see
[`seg_hold()`](https://statmodels7.github.io/modelterms7/reference/seg_hold.md)),
and a developed
[`seg()`](https://statmodels7.github.io/modelterms7/reference/seg.md) is
rejected.

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
