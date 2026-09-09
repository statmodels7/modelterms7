# Report an Argument That s() and te() No Longer Take

Names the replacement for each of the arguments the smooth constructors
carried before the construction became an object. They are reported here
rather than left to R's own "unused argument", which names the argument
and not what to write instead.

## Usage

``` r
retired_smooth_args(dots, fn)
```

## Arguments

- dots:

  The `...` of the constructor.

- fn:

  `"s"` or `"te"`, for the message.

## Value

`NULL`, invisibly. Called for the error it signals.

## Details

The four decisions a smooth is made of – the basis, the penalty, the
null space and the reparametrization – now travel together on a basis7
smoother, so the settings that described the first two belong there. `k`
and `degree` are a B-spline's and were never general; `linear` named the
null space as "linear", which it is only when the penalty is of order 2.
