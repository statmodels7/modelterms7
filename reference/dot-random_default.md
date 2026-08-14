# The Centered Gaussian Defaults

`gaussian1_distrib` at one within-group column, so the hyperparameter IS
the standard deviation of the effects, and the multivariate Gaussian on
an unstructured covariance for several correlated ones.

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

The structure's role is DECLARED rather than left at "either": a
structure that does not say which matrix of the prior it is cannot be
read as either, the two differing in the sign of the log-determinant
term.
