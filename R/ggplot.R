#' The autoplot() Generic of ggplot2
#'
#' @description
#' The generic `autoplot()` belongs to \pkg{ggplot2}, which this package
#' suggests and does not import. The methods written on it here are
#' registered by `S7::methods_register()` when \pkg{ggplot2} is loaded.
#'
#' @return An S7 external generic, the description of `ggplot2::autoplot()`
#'   that `S7::method()` registers methods on.
#'
#' @keywords internal
autoplot <- S7::new_external_generic("ggplot2", "autoplot", "object")


#' The Packages That a Graph Requires
#'
#' @description
#' Signals an error that names the missing packages when \pkg{ggplot2}, or
#' \pkg{patchwork} with `patchwork = TRUE`, is not installed. The graphs of
#' the package are drawn with \pkg{ggplot2}, which is suggested and not
#' imported.
#'
#' @param patchwork Whether the graph also requires \pkg{patchwork}.
#'
#' @return `NULL`, invisibly.
#'
#' @keywords internal
need_ggplot <- function(patchwork = FALSE) {
  pkgs <- c("ggplot2", if (patchwork) "patchwork")
  miss <- pkgs[!vapply(pkgs, requireNamespace, logical(1), quietly = TRUE)]
  if (length(miss)) {
    stop(sprintf("The graph is drawn with %s; install %s.",
                 paste(miss, collapse = " and "),
                 if (length(miss) == 1L) "it" else "them"), call. = FALSE)
  }
  invisible(NULL)
}
