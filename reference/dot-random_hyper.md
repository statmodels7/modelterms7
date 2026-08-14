# Where a Held Hyperparameter Belongs

Splits the term's `hyper` over its penalty entries and checks every name
against the penalty that carries it.

## Usage

``` r
.random_hyper(entries, hyper, label)
```

## Arguments

- entries:

  The entries from
  [`.random_entries`](https://statmodels7.github.io/modelterms7/reference/dot-random_entries.md).

- hyper:

  The term's `hyper`, already normalized.

- label:

  The term's label, for the message.

## Value

The entries, with `fixed` filled in and checked.

## Details

A name is qualified by the within-group column where there is one
penalty per column. An unqualified one is an ERROR that lists what there
is, not a value recycled over every column: a caller who wants the same
value everywhere writes it into the distribution, where it stops being a
hyperparameter at all, and silent recycling is the trap this file's
history records for `ifelse`.
