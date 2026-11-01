# How a Term's Block Moves Along Four Directions

Returns the fourth derivative of a term's block along four directions,
one entry per observation and column of the block, the companion of
[`term_block_deriv3()`](https://statmodels7.github.io/modelterms7/reference/term_block_deriv3.md)
at the next order.

## Usage

``` r
term_block_deriv4(term, coef = NULL, v, u, s, r, ...)
```

## Arguments

- term:

  A built term.

- coef:

  The coefficients, one per coefficient of the term, or `NULL` for the
  coefficients that the term carries.

- v, u, s, r:

  Numeric vectors as long as the term's coefficients.

- ...:

  Passed to methods.

## Value

A numeric matrix, one row per observation and one column per coefficient
of the term.

## Details

The base method returns zeros. For
[`nl()`](https://statmodels7.github.io/modelterms7/reference/nl.md) the
entry is the fifth derivative of the predictor along the column and the
four directions, written out over the 52 partitions of the five: grouped
by the block that holds the column, of size five (one partition), four
(four), three (six with the other two directions together and six with
them apart), two (four choices, each with the other three directions in
one block, in a pair and a single in three ways, or apart) and one (the
fifteen partitions of the four directions). Each block of a partition
carries the derivative of one inverse link to its size and each
partition the derivative of \\f\\ to its number of blocks. A parameter
whose development moves is rejected, as in
[`term_block_deriv3()`](https://statmodels7.github.io/modelterms7/reference/term_block_deriv3.md).

## See also

[`term_block_deriv3()`](https://statmodels7.github.io/modelterms7/reference/term_block_deriv3.md)

## Examples

``` r
dd <- data.frame(x = seq(0.2, 3, length.out = 20))
dd$y <- 2 * exp(-1.3 * dd$x)
b <- term_build(nl(~ a * exp(-r * x), start = list(a = 2, r = 1.3)), dd)
dim(term_block_deriv4(b, v = c(1, 0), u = c(0, 1), s = c(1, 1),
                      r = c(1, -1)))
#> [1] 20  2
```
