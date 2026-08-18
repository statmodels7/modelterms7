#' @keywords internal
#' @useDynLib modelterms7, .registration = TRUE
#' @importFrom Rcpp sourceCpp
#' @importFrom RcppParallel RcppParallelLibs
#' @importFrom graphics abline axis par plot segments
"_PACKAGE"

# Registers the S7 methods written on base generics (print) with the S3
# dispatch table of the installed namespace; without it they are found
# under pkgload and silently absent from an installed package.
#' @noRd
.onLoad <- function(...) {
  S7::methods_register()
}
