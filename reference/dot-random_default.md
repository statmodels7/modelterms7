# The Centered Gaussian Defaults

Returns `gaussian1_distrib` at one within-group column, so the
hyperparameter is the standard deviation of the effects, and the
multivariate Gaussian on an unstructured covariance for several
correlated ones.

## Usage

``` r
.random_default(d, correlated)
```

## Arguments

- d:

  The number of within-group columns.

- correlated:

  Whether the effects may correlate.

## Value

A distributions7 object.

## Details

The matrix parametrization is the family's, and the family built here is
`mvgaussian1_distrib()`, so the structure is the effects' covariance and
its hyperparameters are the log-Cholesky coordinates of that covariance:
the logarithms of the diagonal entries of its Cholesky factor, and the
off-diagonal entries of that factor. The precision reading is a
different prior, because the two differ in the sign of the
log-determinant term, and it is written by passing
`mvgaussian2_distrib()` to `distrib` instead.
