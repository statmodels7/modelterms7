# A Jet over a Set of Directions, and Its Product Rule

`.gas_jetn()` holds the derivatives of one factor of the recursion along
every subset of up to three directions, each with two free slots, and
`.gas_prodn()` returns the derivative of the product of two such factors
along all the directions, with two free slots.

## Usage

``` r
.gas_jetn(v0, d1, mats, dirs)

.gas_subsets(nd)

.gas_jet_get(A, S, free)

.gas_prodn(A, B)

.gas_prodS(A, B, S)

.gas_subsets_of(S)

.gas_jet_key(S)
```

## Arguments

- v0:

  The value.

- d1:

  The first derivative, a vector.

- mats:

  A list of matrices named by subset: `"s"` for the second derivative,
  `"s1"`, `"s2"`, `"s12"` and so on, the digits being the positions of
  the directions in increasing order.

- dirs:

  The directions.

- A, B:

  Two jets over the same directions.

- S:

  The subset of the directions along which to differentiate, by
  position.

## Value

`.gas_jetn()` a jet; `.gas_prodn()` a matrix.

## Details

A jet has the value, the first derivative and one matrix per subset
\\S\\ of the directions, the derivative of order \\\|S\| + 2\\ along
\\S\\ with two free slots. The derivative along \\S\\ with one free slot
is the matrix of \\S\\ without one of its directions applied to that
direction, and with no free slot it is that vector applied to another
direction. The product rule deals the directions and the two free slots
to the two factors: \$\$D_S(AB) = \sum\_{T\subseteq
S}\big(A\_{T,2}B\_{S\setminus T,0} + A\_{T,0}B\_{S\setminus T,2} +
A\_{T,1}B\_{S\setminus T,1}^{\top} + B\_{S\setminus
T,1}A\_{T,1}^{\top}\big),\$\$ which is 16 terms with two directions, as
[`.gas_prod4()`](https://statmodels7.github.io/modelterms7/reference/dot-gas_prod4.md)
writes them, and 32 with three.
