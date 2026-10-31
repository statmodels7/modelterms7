# The Product of Two Jets, as a Jet

Returns the jet of \\AB\\: its value, its first derivative and its
matrix along every subset of the directions.

## Usage

``` r
.gas_jet_times(A, B)

.gas_jet_add(A, B, w = 1)

.gas_jet_const(v0, d1, m, dirs)
```

## Arguments

- A, B:

  Two jets over the same directions.

- w:

  A number, the weight of the second jet.

- v0, d1:

  A value and a first derivative.

- m:

  The number of coordinates.

- dirs:

  The directions.

## Value

A jet.
