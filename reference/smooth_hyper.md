# The Held Smoothing Parameters of a Smooth

Carries the `lambda` a smooth's constructor was given onto the names its
penalty uses, and checks it against the positivity every smoothing
parameter obeys.

## Usage

``` r
smooth_hyper(lambda, names, what = "this smooth")
```

## Arguments

- lambda:

  What the constructor was given, or `NULL`.

- names:

  The penalty's own hyperparameter names.

- what:

  The term's label, for the message.

## Value

A named list of held values.

## Details

A one-dimensional smooth and an isotropic tensor product carry one
smoothing parameter, named `lambda`; an anisotropic tensor product
carries one per margin, `lambda1`, `lambda2` and so on, which is what
[`additive_penalty`](https://statmodels7.github.io/penalties7/reference/additive_penalty.html)
names them. An unnamed vector is read in that order and must be as long
as there are margins; a named one may hold some and leave the others to
be estimated.

## See also

[`check_hyper`](https://statmodels7.github.io/modelterms7/reference/check_hyper.md),
[`s`](https://statmodels7.github.io/modelterms7/reference/s.md),
[`te`](https://statmodels7.github.io/modelterms7/reference/te.md)
