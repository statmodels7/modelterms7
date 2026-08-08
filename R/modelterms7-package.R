#' @keywords internal
"_PACKAGE"

# Registers the S7 methods written on base generics (print) with the S3
# dispatch table of the installed namespace; without it they are found
# under pkgload and silently absent from an installed package.
#' @noRd
.onLoad <- function(...) {
  S7::methods_register()
}
