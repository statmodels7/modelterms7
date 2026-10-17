# Where a Held Hyperparameter Belongs

Splits a term's `hyper` over the entries that the term declares through
[`term_penalties()`](https://statmodels7.github.io/modelterms7/reference/term_penalties.md)
and checks every name against the penalty that carries it. The same
applies to `id`.

## Usage

``` r
.entry_hyper(entries, hyper, ids, label, what = "the effects' distribution")
```

## Arguments

- entries:

  The entries that a term declares, where each entry is a list with at
  least `name` and `penalty`.

- hyper:

  The term's `hyper`, already normalized.

- ids:

  The term's `ids`, which are checked against the same names.

- label:

  The term's label, for the message.

- what:

  The object to which the hyperparameters belong, named in the message.

## Value

The entries, with `fixed`, `values` and `ids` filled in and checked.

## Details

The hyperparameters of a term are not always determined when the term is
written. A random effect's are the effects' distribution's, and a
smooth's are its penalty's, which may come from a factory called at a
coefficient count that the data determine. This function therefore makes
the check at the build, which is the first point at which the penalties
exist.

A name is qualified by the entry to which it belongs where there is more
than one: the within-group column for a random effect, the level of the
factor for a smooth with one smoothing parameter per level. An
unqualified name is an error that lists the names available. A value is
never recycled over the entries, because a value applied silently to
every entry would be wrong for the entries to which it does not belong.
