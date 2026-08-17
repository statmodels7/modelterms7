# Settle Whether a Block Is Built Sparse

Passes an explicit `TRUE` or `FALSE` through, and where the caller left
`NULL` decides from the size of the block.

## Usage

``` r
.resolve_sparse(sparse, n, ncol_ind)
```

## Arguments

- sparse:

  `TRUE`, `FALSE`, or `NULL` to decide here.

- n:

  The number of rows.

- ncol_ind:

  The columns coming from indicators, from
  [`.indicator_cols`](https://statmodels7.github.io/modelterms7/reference/dot-indicator_cols.md).

## Value

A single logical.

## Details

The dense indicator part holds `n * ncol_ind` cells where the sparse one
holds one non-zero per row, so that product is what the two routes are
separated by, and the threshold is read off it rather than off a count
of levels. Measured end to end on `y ~ 0 + g + s(x)` over eighteen
combinations of sample size and level count, the routes cross at about
\\10^5\\ cells: at \\n = 1000\\ the sparse route loses at every level
count up to sixty (\\6 \times 10^4\\ cells, 0.93 times the dense route),
at \\n = 5000\\ it crosses between fifteen and twenty-five levels, and
at \\n = 20000\\ between six and ten. The same threshold accounts for
the large cases: four hundred levels at \\n = 20000\\ are \\8 \times
10^6\\ cells and run 43.75 times faster, with the log-likelihood
identical.

A design carrying no factor has no indicator part, so the product is
zero and the block is built dense, which is what the measurements ask
for there (0.66 to 0.90 times the dense route on purely continuous
covariates).

## See also

[`.indicator_cols`](https://statmodels7.github.io/modelterms7/reference/dot-indicator_cols.md),
[`linpar`](https://statmodels7.github.io/modelterms7/reference/linpar.md)
