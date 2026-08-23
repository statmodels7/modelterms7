# The Columns of a Component That Belong to Each of Its Sub-Terms

Splits a developed parameter's columns among the sub-terms developing
it, in the order the block binds them.

## Usage

``` r
component_sub_index(index, subs)
```

## Arguments

- index:

  The component's columns in the term's block.

- subs:

  The sub-terms developing the parameter.

## Value

A list of integer vectors, one per sub-term, empty where there are no
sub-terms.

## Details

A developed parameter's block is its sub-terms' blocks bound side by
side in the order they were given, so the division is their coefficient
counts cumulated. It is computed by the term rather than left to a
consumer because it rests on how the block is assembled, which is the
term's business and not something a name can be parsed for.
