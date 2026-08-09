# The Progress of a Break-Point Iteration

`seg_step` returns how far each break-point moved at the last call to
[`term_refresh`](https://statmodels7.github.io/modelterms7/reference/term_refresh.md),
and `seg_converged` compares the largest of those with the tolerance of
fasola2018, a hundredth of the distance between consecutive distinct
observations of the covariate. A term that has been built but not yet
refreshed has taken no step, so `seg_step` returns `NA` and
`seg_converged` returns `FALSE`.

## Usage

``` r
seg_step(term)

seg_converged(term)
```

## Arguments

- term:

  A built
  [`SegTerm`](https://statmodels7.github.io/modelterms7/reference/SegTerm.md).

## Value

`seg_step` returns a numeric vector with one entry per break-point and
per level of `by`; `seg_converged` returns a single logical.

## Details

With \\x\_{(1)} \< \cdots \< x\_{(m)}\\ the distinct covariate values,
the run stops at

\$\$\max_k \lvert \psi_k^{(t)} - \psi_k^{(t-1)} \rvert \< \Delta, \qquad
\Delta = 0.01 \cdot \operatorname{median}\_{i}\\ (x\_{(i+1)} -
x\_{(i)}).\$\$

fasola2018 take the smallest of those gaps rather than their median,
which agrees with this on the evenly spaced covariates of their examples
and is of order \\m^{-2}\\ on a random one, hence unreachable.

The rule is one of resolution: below that distance the objective of a
discontinuous term, a step function of the break-point, cannot change.
It therefore tightens as the sample grows while the precision the fixed
point is reached at does not, so on a large sample the last iterations
move the break-point by a little more than the rule allows and the run
continues past the point where the estimate has settled. A caller that
can evaluate the objective should stop on its relative change instead,
which is what `segmented` does and what costs a continuous term nothing:
its iteration can settle into a cycle of period two in the break-point,
in which case this rule is never met while the objective has long since
stopped moving.

## See also

[`seg`](https://statmodels7.github.io/modelterms7/reference/seg.md),
[`seg_psi`](https://statmodels7.github.io/modelterms7/reference/seg_psi.md),
[`seg_start`](https://statmodels7.github.io/modelterms7/reference/seg_start.md)
[`seg_start`](https://statmodels7.github.io/modelterms7/reference/seg_start.md)

## Examples

``` r
set.seed(1)
dd <- data.frame(x = sort(runif(100, 0, 10)))
dd$y <- 2 * (dd$x > 6) + rnorm(100, sd = 0.3)
b <- term_build(jump(x, psi = 4, linear = FALSE), dd)
cf <- b@blueprint$coef
for (it in 1:30) {
  b <- term_refresh(b, cf)
  X <- term_matrix(b)
  cf <- as.numeric(qr.solve(crossprod(X), crossprod(X, dd$y)))
  if (seg_converged(b)) break
}
c(psi = seg_psi(b, cf), step = seg_step(b))
#>         psi        step 
#> 5.996835956 0.001012559 
```
