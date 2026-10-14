# Build a Smooth Term

Builds the basis of an
[`s()`](https://statmodels7.github.io/modelterms7/reference/s.md) or
[`te()`](https://statmodels7.github.io/modelterms7/reference/te.md) term
at the observed covariates, applies the reparametrization that
construction calls for, multiplies by a `by` variable where there is
one, and attaches the roughness penalty. Two quantities are computed
**from the data** here and recorded in the blueprint, so that
[`term_predict()`](https://statmodels7.github.io/modelterms7/reference/term_predict.md)
reapplies them instead of deriving them again: the Demmler-Reinsch
transform for
[`s()`](https://statmodels7.github.io/modelterms7/reference/s.md), and
the centering constraint for
[`te()`](https://statmodels7.github.io/modelterms7/reference/te.md).

## Arguments

- term:

  An unbuilt or built
  [`SmoothTerm()`](https://statmodels7.github.io/modelterms7/reference/SmoothTerm.md).

- data:

  A data frame carrying the covariates and the `by` variable.

- ...:

  Unused.

## Value

The term with `X`, `coef_names` and `blueprint` filled, and with
`penalty` filled except where the term declares penalties of its own (a
penalty factory on the smoother, or `by_hyper = "level"`), in which case
`penalty` is `NULL` and
[`term_penalties()`](https://statmodels7.github.io/modelterms7/reference/term_penalties.md)
returns them.

## The default basis

Each margin is built by its smoother,
[`basis7::bspline_smooth()`](https://statmodels7.github.io/basis7/reference/bspline_smooth.html)
by default. Unless the smoother sets `lower` and `upper`, the basis is
placed on the observed range of the covariate, padded by a thousandth of
that range at each end so that the extreme observations are strictly
inside.

## One covariate

The block is built by
[`basis7::smoother_build()`](https://statmodels7.github.io/basis7/reference/smoother_build.html),
which takes the Demmler-Reinsch inner product at the observed values.
With the default
[`basis7::bspline_smooth()`](https://statmodels7.github.io/basis7/reference/bspline_smooth.html)
(`null_space = "keep"`) a column \\(x - \bar{x})/s_x\\ is prepended and
the penalty is \\\mathrm{diag}(0, 1, \dots, 1)\\; the center and the
scale go into the blueprint with the transform. With
`null_space = "drop"` that column is absent and the penalty is the
identity.

## Several covariates

[`basis7::tensor_basis()`](https://statmodels7.github.io/basis7/reference/tensor_basis.html)
combines the margins, and one penalty component per margin is carried
into the product as \\I \otimes \cdots \otimes P_v \otimes \cdots
\otimes I\\, each \\P_v\\ the margin's roughness matrix at the order and
the measure that its smoother carries (the second-derivative Gram at the
defaults), normalized to a maximum entry of one.

The Kronecker product is taken over the **reversed** blocks, because
[`basis7::tensor_basis()`](https://statmodels7.github.io/basis7/reference/tensor_basis.html)
varies the first margin fastest. The block is then centered by
[`basis7::constrain_basis()`](https://statmodels7.github.io/basis7/reference/constrain_basis.html)
over the observed covariates, so it has one column fewer than the
product of the marginal dimensions.

## The `by` variable

A factor `by` interacts the basis with the level indicators, giving `m`
copies of the block and a penalty repeated blockwise; the levels are
recorded so a prediction uses the same set. A numeric `by` multiplies
the basis after an unpenalized constant column, named `const`, has been
put in front of it, and the block is not centered. The storage is
settled here and recorded, an explicit `sparse = TRUE` being rejected
where there is no factor `by`.

## See also

[`s()`](https://statmodels7.github.io/modelterms7/reference/s.md) and
[`te()`](https://statmodels7.github.io/modelterms7/reference/te.md) for
the two constructions,
[`term_predict.SmoothTerm()`](https://statmodels7.github.io/modelterms7/reference/term_predict.SmoothTerm.md)
for the block at new rows.

## Examples

``` r
set.seed(1)
dd <- data.frame(x = sort(runif(80)), z = runif(80))

# What the build records for a one-covariate smooth.
b <- term_build(s(x, basis7::bspline_smooth(k = 8)), dd)
names(b@blueprint)
#> [1] "core"      "marg"      "spec"      "vars"      "by"        "by_levels"
#> [7] "sparse"    "penalties" "nblock"   
b@blueprint$core$kind
#> [1] "smoother"

# And for a tensor product: one column fewer than 4 x 4.
bt <- term_build(te(x, z, smooths = basis7::bspline_smooth(k = 4)), dd)
term_npar(bt)
#> [1] 15
```
