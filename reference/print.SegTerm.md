# Print a Break-Point Term

Prints the label, the construction and how many break-points the term
carries, followed by where they currently sit. A smoothed term adds a
line naming its smoother and the transition width, and a term with a
developed coefficient says which one.

## Arguments

- x:

  A
  [`SegTerm()`](https://statmodels7.github.io/modelterms7/reference/SegTerm.md),
  built or not.

- ...:

  Unused, and accepted so that the signature matches
  [`print()`](https://rdrr.io/r/base/print.html)'s.

## Value

`x`, invisibly. Called for the lines it writes.

## Details

The form is

    <SegTerm> 'seg': seg, 1 break-point
      at: 6.501

The positions come from
[`seg_psi()`](https://statmodels7.github.io/modelterms7/reference/seg_psi.md),
so they are the current ones and move with every
[`term_refresh()`](https://statmodels7.github.io/modelterms7/reference/term_refresh.md).
Where a break-point carries a development there is no single number, a
position then being one value per observation, and the range is printed
instead.

The smoother line reports the width because it changes what the term
means: the transition is a real feature of the model, not an
implementation detail, and a fit at one width is not a fit at another.

## See also

[`seg_psi()`](https://statmodels7.github.io/modelterms7/reference/seg_psi.md)
for the positions,
[`seg()`](https://statmodels7.github.io/modelterms7/reference/seg.md),
[`jump()`](https://statmodels7.github.io/modelterms7/reference/jump.md)
and
[`jseg()`](https://statmodels7.github.io/modelterms7/reference/jseg.md).

## Examples

``` r
set.seed(1)
d <- data.frame(x = sort(runif(120, 0, 10)), id = factor(rep(1:4, each = 30)))
d$y <- 1 + 0.5 * d$x + 2 * pmax(d$x - 6, 0) + rnorm(120, sd = 0.4)

# A specification, and the same term built at its starting quantile.
seg(x, npsi = 2)
#> <SegTerm> 'seg': seg (specification)
term_build(seg(x, npsi = 2), d)
#> <SegTerm> 'seg': seg, 2 break-points
#>   at: 3.774, 6.642

# A smoothed term names its smoother and its width.
term_build(jump(x, smoothed = penalties7::smooth_probit()), d)
#> <SegTerm> 'jump': jump, 1 break-point
#>   smoothed (probit, h = 0.0607)
#>   at: 4.803

# A developed break-point has a range instead of a number.
term_build(seg(x, psi ~ id), d)
#> <SegTerm> 'seg': seg, 1 break-point; psi1 developed
#>   at (developed): [4.803, 4.803]
```
