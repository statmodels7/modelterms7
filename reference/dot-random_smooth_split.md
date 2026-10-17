# The Null Space and the Penalized Coordinates of a Smooth

Splits the coefficients of a smooth into the directions on which every
component of its penalty vanishes and the directions that the penalty
covers, and returns the linear map from the smooth's coefficients to
those two groups.

## Usage

``` r
.random_smooth_split(comps, names, tol = 1e-10)
```

## Arguments

- comps:

  The component matrices, as
  [`.random_smooth_components()`](https://statmodels7.github.io/modelterms7/reference/dot-random_smooth_components.md)
  returns them.

- names:

  The smooth's coordinate names, without its label.

- tol:

  The relative eigenvalue tolerance for a zero.

## Value

A list with `T`, the map (\\p\\ rows, one column per coordinate of the
random effect, null space first); `nnull`, the dimension of the null
space; `npen`, the number of penalized coordinates; `comps`, the
components on the penalized coordinates; and `names`.

## Details

The null space of a sum of positive semidefinite matrices is the
intersection of the components' null spaces, so it is read once from the
components normalized to a unit largest entry and does not move with the
smoothing parameters.

Where the smoother has already placed the null space in its leading
columns with a zero border in every component, as every basis7 smoother
does under `null_space = "keep"`, the map keeps those columns and their
names. Otherwise, as for a tensor product, an orthonormal basis of the
null space and one of its complement are taken from the eigenvectors of
the normalized sum, and the columns are named `null1`, ... and `z1`,
....

With a single component the penalized coordinates are then whitened:
with \\P = R^\top R\\ the penalty on the penalized coordinates
\\\gamma\\, the coordinates \\\delta = R\gamma\\ carry the identity
penalty, so the random effect over them is independent and Gaussian with
one standard deviation. A Demmler-Reinsch smoother already has the
identity there, and its coordinates are unchanged. With several
components the penalized coordinates are kept and each component is
carried onto them.
