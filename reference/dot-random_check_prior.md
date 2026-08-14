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

What is rejected is a FREE location: it is confounded with the intercept
of the equation the term sits in, which is a flat direction and not a
model. A location HELD at a value is identified, whatever the value, and
is not policed – it shrinks the effects towards that value, which is a
modelling statement rather than a defect. Nor could the value be policed
in general: where the prior is a transformation of another family, the
parameter interpreted as its mean is the mean on the ORIGINAL scale, and
holding the mean of a gamma at one is what centers its logarithm.
