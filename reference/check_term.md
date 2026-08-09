# Numerical Validation of a Model Term

Runs a battery of structural checks on a term specification against a
data frame: that it builds, that the block's dimensions, names and count
agree, that the smoothness flag is a logical scalar, that
[`term_predict`](https://statmodels7.github.io/modelterms7/reference/term_predict.md)
on the same data reproduces the block exactly, and that prediction on a
subset of rows equals the corresponding rows of the block. The last
check is the blueprint's: a term that re-derives factor levels from the
new data instead of reusing the levels recorded at build time fails it
as soon as the subset drops a level, which is why the subset is chosen
to drop one whenever the data carry a factor.

## Usage

``` r
check_term(term, data, verbose = TRUE)
```

## Arguments

- term:

  A term specification (an object inheriting from
  [`model_term`](https://statmodels7.github.io/modelterms7/reference/model_term.md)).

- data:

  A data frame.

- verbose:

  Logical; print one line per check.

## Value

Invisibly, a data frame with columns `check`, `status` (`"OK"` or
`"FAILED"`) and `info`.

## Details

Writing \\X = \\ `term_matrix(term_build(term, data))`, the two
identities checked are

\$\$\texttt{term\\predict}(\text{term}, \text{data}) = X, \qquad
\texttt{term\\predict}(\text{term}, \text{data}\[S, \]) = X\[S, \],\$\$

for a row subset \\S\\. The second does not follow from the first: a
term that rebuilds its encoding from the rows it is given satisfies the
first, since the rows are then the same, and fails the second as soon as
\\S\\ omits a factor level or narrows the range a basis is placed on.
The subset is dropped of its unused levels before it is passed, so a
plain row subset cannot pass by carrying the original levels along with
it.

## Examples

``` r
dd <- data.frame(x = 1:6, g = factor(rep(c("a", "b", "c"), 2)))
res <- check_term(linpar(~ x + g), dd, verbose = FALSE)
all(res$status == "OK")
#> [1] TRUE
```
