# What a Prior on the Effects Has to Satisfy

Checks a caller-supplied effects distribution where the term is built
and names what is wrong with it.

## Usage

``` r
.random_check_prior(d, dim_needed, what)
```

## Arguments

- d:

  A distributions7 object.

- dim_needed:

  The number of within-group columns.

- what:

  How to name the argument in a message.

## Value

Invisibly `TRUE`; called for its errors.

## Details

What is rejected is a free location. It is confounded with the intercept
of the equation the term sits in, leaving a flat direction along which
the fit has no answer.

A location held at a value is identified whatever that value is, and it
is not policed: it shrinks the effects toward that value, which is a
modeling statement. Nor could the value be policed in general. Where the
prior is a transformation of another family, the parameter interpreted
as its mean is the mean on the original scale, and holding the mean of a
gamma at one is what centers its logarithm.
