# Print a Model Term

Prints one line naming the term's class, its label when it has one, and
whether it is built. A built term reports how many coefficients its
block carries; a specification says so and names the call that would
build it. This is the default for every term class, and several classes
override it to add what they alone carry.

## Arguments

- x:

  A term, built or not.

- ...:

  Unused, and accepted so that the signature matches
  [`print()`](https://rdrr.io/r/base/print.html)'s.

## Value

`x`, invisibly. Called for the line it writes.

## Details

The two forms are

    <SmoothTerm> 's(x)' built: 4 coefficients
    <SmoothTerm> 's(x)' (specification; call term_build() with data)

The class name is the S7 class's own, the label comes from the `label`
property and is omitted when empty, and the count is `ncol(x@X)`.

A structural term has no `X` to count, so it gets the word `built` and
no count. Every shipped structural term registers a method of its own,
so this line is reached only by a structural class written outside the
package.

The classes that override this add something of their own:
[`print.PenalizedTerm()`](https://statmodels7.github.io/modelterms7/reference/print.PenalizedTerm.md)
names the penalty and the standardization,
[`gas()`](https://statmodels7.github.io/modelterms7/reference/gas.md)'s
and
[`regime()`](https://statmodels7.github.io/modelterms7/reference/regime.md)'s
report their parameters,
[`nl()`](https://statmodels7.github.io/modelterms7/reference/nl.md)'s
its formula and
[`seg()`](https://statmodels7.github.io/modelterms7/reference/seg.md)'s
its break-points.

## See also

[`term_is_built()`](https://statmodels7.github.io/modelterms7/reference/term_is_built.md)
for the predicate it branches on,
[`term_coef_names()`](https://statmodels7.github.io/modelterms7/reference/term_coef_names.md)
for the names behind the count.

## Examples

``` r
d <- data.frame(x = 1:4)

# A specification, and the same term built.
linpar(~ x)
#> <LinparTerm> (specification; call term_build() with data)
term_build(linpar(~ x), d)
#> <LinparTerm> built: 2 coefficients

# The label is shown when there is one.
s(x, k = 5)
#> <SmoothTerm> 's(x)' (specification; call term_build() with data)
linpar(~ x, label = "lin")
#> <LinparTerm> 'lin' (specification; call term_build() with data)
```
