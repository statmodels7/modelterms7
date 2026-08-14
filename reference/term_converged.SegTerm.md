# Whether a Break-Point Term's Break-Points Have Settled

[`seg_converged`](https://statmodels7.github.io/modelterms7/reference/seg_step.md),
so that a fitting layer reads the rule the construction is stopped on
without knowing it is holding a break-point term. A term that has not
been refreshed has not moved and reports `FALSE`, the first refresh
being the one that measures nothing.

## Arguments

- term:

  A built
  [`SegTerm`](https://statmodels7.github.io/modelterms7/reference/SegTerm.md).

- ...:

  Unused.

## Value

A single logical.
