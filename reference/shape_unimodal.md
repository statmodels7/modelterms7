# A Unimodal Shape

Constructs the shape that makes a smooth increase up to a mode and
decrease after it, for
[`scs()`](https://statmodels7.github.io/modelterms7/reference/scs.md).

## Usage

``` r
shape_unimodal(link = linkfunctions7::log_link(), mode = NULL, index = NULL)
```

## Arguments

- link:

  A linkfunctions7 link onto \\(0, \infty)\\.

- mode:

  `NULL`, or a single number: the value of the covariate at which the
  mode is placed.

- index:

  `NULL`, or the index \\M\\ itself, which a fitting layer sets on the
  versions it tries. It takes precedence over `mode`.

## Value

An S7 object of class `ShapeUnimodal`, inheriting from
[shape_constraint](https://statmodels7.github.io/modelterms7/reference/shape_constraint.md).

## Details

The mode sits at a Greville abscissa \\\xi_M\\ of the basis: the slopes
of the coefficients are \\+g^{-1}(\gamma_j)\\ up to index \\M\\ and
\\-g^{-1}(\gamma_j)\\ after it, so the coefficients, and with them the
curve, rise to a single maximum. \\M = 1\\ is a decreasing curve and \\M
= K\\ an increasing one. The penalty of
[`scs()`](https://statmodels7.github.io/modelterms7/reference/scs.md) is
one random walk on \\\gamma\\ through the mode, whose null space is the
same for every \\M\\, so fits at different \\M\\ are compared by their
criterion.

With `mode` given, \\M\\ is the index of the Greville abscissa nearest
to it and the model is fitted once. With `mode = NULL` the fitting layer
chooses \\M\\: a pilot fit with
[`s()`](https://statmodels7.github.io/modelterms7/reference/s.md) in
place of the term gives the abscissa nearest to its maximum, \\M_0\\,
and the model is fitted at the indices \\M_0 - 2, \dots, M_0 + 2\\; the
best criterion decides (see
[`term_pilot()`](https://statmodels7.github.io/modelterms7/reference/term_pilot.md)
and
[`term_variants()`](https://statmodels7.github.io/modelterms7/reference/term_variants.md)).

## See also

[`scs()`](https://statmodels7.github.io/modelterms7/reference/scs.md),
[`shape_increasing()`](https://statmodels7.github.io/modelterms7/reference/shape_constructors.md).

## Examples

``` r
shape_unimodal()
#> <shape: unimodal, link = log>
shape_unimodal(mode = 0.4)
#> <shape: unimodal, link = log, mode = 0.4>
```
