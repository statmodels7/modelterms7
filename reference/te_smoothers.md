# One Smoother per Margin of a Tensor Product

Resolves the `smooths` argument of
[`te()`](https://statmodels7.github.io/modelterms7/reference/te.md) into
one smoother per covariate, and rejects a margin asking for what the
product does not read from it.

## Usage

``` r
te_smoothers(smooths, nv)
```

## Arguments

- smooths:

  One smoother, or a list of one per covariate.

- nv:

  How many covariates.

## Value

A list of `nv` smoothers.

## Details

The product reads a margin's basis and its roughness matrix, so `k`,
`degree`, `order`, `measure` and the interval all act. The constraint,
the null space and the reparametrization are the **product's**: the
tensor basis contains the constant whatever its margins do, so the block
carries the sum-to-zero constraint over the observed covariates rather
than a marginal one, and the marginal linear effects are not separated
out as [`s()`](https://statmodels7.github.io/modelterms7/reference/s.md)
separates its one.

Those three are therefore rejected at a non-default value rather than
accepted and ignored, which would report a fit of a model the caller did
not ask for.
