# Where a Held Hyperparameter Belongs

Splits a term's `hyper` over the entries it declares through
[`term_penalties()`](https://statmodels7.github.io/modelterms7/reference/term_penalties.md)
and checks every name against the penalty that carries it. The same
applies to `id`.

## Usage

``` r
.entry_hyper(entries, hyper, ids, label, what = "the effects' distribution")
```

## Arguments

- entries:

  The entries a term declares, each a list with at least `name` and
  `penalty`.

- hyper:

  The term's `hyper`, already normalized.

- ids:

  The term's `ids`, which are checked against the same names.

- label:

  The term's label, for the message.

- what:

  What the hyperparameters belong to, named in the message.

## Value

The entries, with `fixed`, `values` and `ids` filled in and checked.

## Details

Which hyperparameters a term has is not always known when the term is
written. A random effect's are the effects' distribution's, and a
smooth's are its penalty's, which may come from a factory called at a
coefficient count the data settle. The check is therefore made here, at
the build, which is the first point at which the penalties exist.

A name is qualified by the entry it belongs to where there is more than
one: the within-group column for a random effect, the level of the
factor for a smooth with one smoothing parameter per level. An
unqualified name is an error that lists what there is, not a value
recycled over every entry, since silent recycling is the trap this
package's history records for `ifelse`.
