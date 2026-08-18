# Move a Break-Point Term to Given Positions

The term with its break-points placed at `psi`, the changes kept, the
scaling schedule fresh and the block rebuilt – ready to iterate from
there. A restart proposal needs it: a bootstrap resample perturbs the
objective by \\1/\sqrt{n}\\ and stops escaping a deep basin as the
sample grows, so the restarting loop also proposes fresh positions
outright, which this is the operation behind.

## Usage

``` r
seg_relocate(term, psi)
```

## Arguments

- term:

  A built break-point term (see
  [`term_build`](https://statmodels7.github.io/modelterms7/reference/term_build.md)).

- psi:

  A numeric vector, one position per break-point.

## Value

The term at the new positions.

## Details

The positions are sorted and confined to the term's own interval. For a
discontinuous construction the auxiliary coefficients are set to \\g_k =
-\delta_k\psi_k\\, a change of level at zero replaced by one, so the
read-off returns exactly the positions given. A term whose
per-break-point coefficients carry a development is rejected: its
positions are one per observation and a single vector does not place
them.

## See also

[`seg_reheat`](https://statmodels7.github.io/modelterms7/reference/seg_reheat.md),
[`seg_start`](https://statmodels7.github.io/modelterms7/reference/seg_start.md)

## Examples

``` r
set.seed(1)
dd <- data.frame(x = sort(runif(100, 0, 10)))
b <- term_build(jump(x, psi = 4), dd)
seg_psi(seg_relocate(b, 6))
#> [1] 6
```
