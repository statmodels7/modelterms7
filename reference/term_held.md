# Coefficients Held by a Term

Returns the positions, among the coefficients of a term, of those that
the term holds fixed. The fitting layer excludes these coefficients from
the solve, from the information matrix that is inverted and from the
reading of the mode. The base method returns an empty vector, since all
the coefficients of an ordinary term are estimated from the derivative
of the objective.

## Usage

``` r
term_held(term, ...)
```

## Arguments

- term:

  A term.

- ...:

  Passed to methods.

## Value

An integer vector of positions in the term's coefficients, possibly
empty.

## Details

The one shipped term with a method of its own is a
[`jump()`](https://statmodels7.github.io/modelterms7/reference/jump.md)
or
[`jseg()`](https://statmodels7.github.io/modelterms7/reference/jseg.md)
term put in its held state by
[`seg_hold()`](https://statmodels7.github.io/modelterms7/reference/seg_hold.md).
Its break-points are fixed at a minimum of the profile objective. That
objective is constant between consecutive observations, so it has no
derivative from which to estimate the positions, and the coefficient
slot that carries each position in the working construction becomes a
column of zeros.

## See also

[`seg_hold()`](https://statmodels7.github.io/modelterms7/reference/seg_hold.md),
[`term_kinks()`](https://statmodels7.github.io/modelterms7/reference/term_kinks.md).

## Examples

``` r
set.seed(1)
d <- data.frame(x = sort(runif(100, 0, 10)))
d$y <- 1 + 2 * (d$x > 6) + rnorm(100, sd = 0.3)
b <- term_build(jump(x, psi = 6), d)
term_held(b)
#> integer(0)
term_held(seg_hold(b))
#> [1] 2
```
