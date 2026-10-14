# Check a Term's Path Depth

Validates the fraction of the emptying value a path descends to.

## Usage

``` r
check_min_ratio(v, what = "this term")
```

## Arguments

- v:

  What the constructor was given, or `NULL`.

- what:

  The term's label, for the message.

## Value

A numeric vector of length one, or of length zero.

## Details

The fraction is one number per term, because only the path over the size
of the kink uses it. A bounded hyperparameter is swept over its own
interval, and a shape that does not move the kink over a geometric grid
above its lower bound; a fraction of an emptying value has no meaning on
either axis.

## See also

[`check_grid()`](https://statmodels7.github.io/modelterms7/reference/check_grid.md),
[`term_path_min()`](https://statmodels7.github.io/modelterms7/reference/term_path_min.md)
