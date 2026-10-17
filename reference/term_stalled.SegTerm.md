# Whether a Break-Point Term Has Stalled

Returns `TRUE` for a sharp
[`jump()`](https://statmodels7.github.io/modelterms7/reference/jump.md)
or
[`jseg()`](https://statmodels7.github.io/modelterms7/reference/jseg.md)
term that
[`seg_hold()`](https://statmodels7.github.io/modelterms7/reference/seg_hold.md)
does not hold, when at least one break-point has not settled and every
break-point that has not settled sits at the floor of its scaling
factor, which
[`term_refresh()`](https://statmodels7.github.io/modelterms7/reference/term_refresh.md)
records at each refresh. The method always returns `FALSE` for
[`seg()`](https://statmodels7.github.io/modelterms7/reference/seg.md)
and for a smoothed construction, which have no scaling schedule, for a
term held by
[`seg_hold()`](https://statmodels7.github.io/modelterms7/reference/seg_hold.md),
and for a term not yet refreshed twice, whose step has not been
measured.

## Arguments

- term:

  A built
  [`SegTerm()`](https://statmodels7.github.io/modelterms7/reference/SegTerm.md).

- ...:

  Unused.

## Value

A single logical.
