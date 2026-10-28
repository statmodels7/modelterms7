# Confine the Free Break-Points of a Term Away From the Held Ones

Prepares a continuous
[`seg()`](https://statmodels7.github.io/modelterms7/reference/seg.md)
term for a fit in which some of its coefficients are held at given
values. Each break-point whose position is not held stays in the
interval between the held positions that surround it, with at least
`min_obs` observations between it and each of them, and the coefficients
of a break-point of which a coefficient is held are not exchanged with
those of another when the positions cross.

## Usage

``` r
seg_pin(term, held, start = NULL, min_obs = 2L)
```

## Arguments

- term:

  A built
  [`seg()`](https://statmodels7.github.io/modelterms7/reference/seg.md)
  term whose break-point positions carry no development.

- held:

  A named numeric vector of the coefficients of the term that are held.
  The names are those of `term@coef_names`. The values of the
  coefficients that carry the positions of break-points are positions;
  the other values are not read, and the coefficient of the linear part
  does not affect the break-points.

- start:

  The positions of the break-points from which the fit starts, one for
  each break-point, ordered as the term orders them. The default is the
  positions of the term.

- min_obs:

  The smallest number of observations between a free and a held
  break-point, a positive integer.

## Value

The term, with the intervals of its break-points and the groups within
which they exchange labels recorded.
[`seg_limits()`](https://statmodels7.github.io/modelterms7/reference/seg_limits.md)
returns the intervals.

## Details

Two break-points of a continuous term at one position make two columns
of the block equal. As the free break-point approaches the held one, the
coefficients of the two changes of slope grow without bound in opposite
directions, because the pair then imitates a jump of the level. Without
this function the relabeling of crossed break-points in
[`term_refresh()`](https://statmodels7.github.io/modelterms7/reference/term_refresh.md)
exchanges the coefficients of the held break-point with those of the
free one, so that the value held for the first is written onto the
second and the fit alternates between two states without converging.

A change of slope is estimated from the observations between its
break-point and the next one, so that two observations are the least
that identifies it, and the interval of a free break-point leaves that
many between it and a held one. The interval is the one, among those
that the held positions cut out of the confinement interval of the term,
that contains the position in `start` or, when that position leaves
fewer observations than that, the nearest interval with room. The side
is therefore the one on which the break-point starts, and the
break-point cannot cross a held one during the fit. The held positions
are clamped to the confinement interval of the term. Free break-points
that share an interval exchange their labels when their positions cross,
as in an ordinary fit.

The columns of a free break-point at an end of its interval are zero, as
they are at the confinement limits of the term, so that the coefficient
stays at the end.

## See also

[`seg_limits()`](https://statmodels7.github.io/modelterms7/reference/seg_limits.md),
[`term_refresh()`](https://statmodels7.github.io/modelterms7/reference/term_refresh.md),
[`seg_hold()`](https://statmodels7.github.io/modelterms7/reference/seg_hold.md).

## Examples

``` r
set.seed(1)
d <- data.frame(x = sort(runif(100, 0, 10)))
b <- term_build(seg(x, npsi = 2, psi = c(3, 7)), d)
p <- seg_pin(b, held = c(seg.psi2 = 5), start = c(3, 7))
seg_limits(p)
#>          lower    upper
#> psi1 0.8356852 4.935413
#> psi2 0.8356852 9.084412
```
