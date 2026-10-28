# Print a Regime Term

Prints two lines: the label and the number of regimes, followed by the
term's parameters. The first line of a built term gives the number of
groups over which the recursion runs, and the first line of a
specification has no group count.

## Arguments

- x:

  A
  [`RegimeTerm()`](https://statmodels7.github.io/modelterms7/reference/RegimeTerm.md),
  built or not.

- ...:

  Unused, and accepted so that the signature matches
  [`print()`](https://rdrr.io/r/base/print.html)'s.

## Value

`x`, invisibly. The function is called for the two lines that it writes.

## Details

The two forms are

    <RegimeTerm> 'regime': 2 regimes (specification)
      parameters: level1, gap2, alr1.1, alr2.1

    <RegimeTerm> 'regime': 2 regimes; 1 group(s)
      parameters: level1, gap2, alr1.1, alr2.1

The group count is `length(build_info$order)`, which is 1 unless the
term was given a `by`. The parameter list is
[`term_params()`](https://statmodels7.github.io/modelterms7/reference/term_params.md),
which is the same before and after a build because only `k` determines
it.

The description of a built term gives its count of groups and omits the
word "built".

## See also

[`regime()`](https://statmodels7.github.io/modelterms7/reference/regime.md),
[`term_params()`](https://statmodels7.github.io/modelterms7/reference/term_params.md).

## Examples

``` r
set.seed(1)
dd <- data.frame(t = 1:40, id = rep(1:2, each = 20),
                 y = c(rnorm(20), rnorm(20, 3)))

# The first call prints a specification, and the second prints the same
# term built over two groups.
regime(2)
#> <RegimeTerm> 'regime': 2 regimes (specification)
#>   parameters: level1, gap2, alr1.1, alr2.1
term_build(regime(2, by = id, time = t), dd)
#> <RegimeTerm> 'regime': 2 regimes; 2 group(s)
#>   parameters: level1, gap2, alr1.1, alr2.1

# Three regimes carry one level, two gaps and six log-ratios.
regime(3)
#> <RegimeTerm> 'regime': 3 regimes (specification)
#>   parameters: level1, gap2, gap3, alr1.1, alr1.2, alr2.1, alr2.2, alr3.1, alr3.2
```
