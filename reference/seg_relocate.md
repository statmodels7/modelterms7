# Move a Break-Point Term to Given Positions

Returns the term in which the break-points are placed at `psi`, the
changes are kept, the scaling schedule is reset and the block is
rebuilt, so that the term is ready to iterate from there. A restart
proposal needs this function: a bootstrap resample perturbs the
objective by \\1/\sqrt{n}\\ and stops escaping a deep basin as the
sample grows, so the restarting loop also proposes fresh positions
outright, and this function places the term at them.

## Usage

``` r
seg_relocate(term, psi)
```

## Arguments

- term:

  A built break-point term (see
  [`term_build()`](https://statmodels7.github.io/modelterms7/reference/term_build.md)).

- psi:

  A numeric vector, one finite position per break-point.

## Value

The term at the new positions.

## Details

The positions are sorted and confined to the term's own interval. For a
sharp
[`jump()`](https://statmodels7.github.io/modelterms7/reference/jump.md)
or
[`jseg()`](https://statmodels7.github.io/modelterms7/reference/jseg.md)
the auxiliary coefficients are set to \\g_k = -\delta_k\psi_k\\, where a
change of level smaller than \\10^{-8}\\ in absolute value is replaced
by one first, so that the read-off returns exactly the positions given.
A continuous or smoothed term stores the positions themselves. On a term
held by
[`seg_hold()`](https://statmodels7.github.io/modelterms7/reference/seg_hold.md)
it moves the held positions. A term whose per-break-point coefficients
carry a development is rejected: its positions are one per observation
and a single vector does not place them.

## See also

[`seg_reheat()`](https://statmodels7.github.io/modelterms7/reference/seg_reheat.md),
[`seg_start()`](https://statmodels7.github.io/modelterms7/reference/seg_start.md)

## Examples

``` r
set.seed(1)
dd <- data.frame(x = sort(runif(100, 0, 10)))
b <- term_build(jump(x, psi = 4), dd)
seg_psi(seg_relocate(b, 6))
#> [1] 6
```
