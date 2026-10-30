# The Packages That a Graph Requires

Signals an error that names the missing packages when ggplot2, or
patchwork with `patchwork = TRUE`, is not installed. The graphs of the
package are drawn with ggplot2, which is suggested and not imported.

## Usage

``` r
need_ggplot(patchwork = FALSE)
```

## Arguments

- patchwork:

  Whether the graph also requires patchwork.

## Value

`NULL`, invisibly.
