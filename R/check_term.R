#' @include term_classes.R generics.R
NULL

# What counts as a design block. A term's `X` is `class_any` on purpose, and
# a grouping indicator is built SPARSE -- a row belongs to one group, so the
# density is 1/m and the dense form is quadratic in the wrong place. So the
# question a validator may ask is whether the object is two-dimensional and
# numeric, not whether it is a base matrix: `is.matrix()` is FALSE for every
# Matrix class and would fail a term for being efficient.
.is_block <- function(x) {
  (is.matrix(x) && is.numeric(x)) ||
    (isS4(x) && length(dim(x)) == 2L)
}

# Coerce to a block WITHOUT densifying: a Matrix is already one and is left
# alone, and only something that is neither is sent through as.matrix().
.as_block <- function(x) if (.is_block(x)) x else as.matrix(x)

# Rows whose removal is most likely to expose a blueprint defect: when the
# data carry a factor, dropping every row of one level makes a rebuilt (as
# opposed to reapplied) encoding lose a column, so the subset check below
# fails for exactly the mistake it exists to catch.
.check_subset_rows <- function(data) {
  n <- nrow(data)
  for (nm in names(data)) {
    v <- data[[nm]]
    if (is.factor(v) && nlevels(v) > 1L) {
      keep <- which(v != levels(v)[nlevels(v)])
      if (length(keep) > 0L && length(keep) < n) return(keep)
    }
  }
  seq_len(max(1L, n %/% 2L))
}

#' Numerical Validation of a Model Term
#'
#' @description
#' Runs a battery of structural checks on a term specification against a
#' data frame: that it builds, that the block's dimensions, names and
#' count agree, that the smoothness flag is a logical scalar, that
#' \code{\link{term_predict}} on the same data reproduces the block
#' exactly, and that prediction on a subset of rows equals the
#' corresponding rows of the block. The last check is the blueprint's: a
#' term that re-derives factor levels from the new data instead of reusing
#' the levels recorded at build time fails it as soon as the subset drops
#' a level, which is why the subset is chosen to drop one whenever the
#' data carry a factor.
#'
#' @details
#' Writing \eqn{X = } \code{term_matrix(term_build(term, data))}, the two
#' identities checked are
#'
#' \deqn{\texttt{term\_predict}(\text{term}, \text{data}) = X,
#'   \qquad
#'   \texttt{term\_predict}(\text{term}, \text{data}[S, ]) = X[S, ],}
#'
#' for a row subset \eqn{S}. The second does not follow from the first: a
#' term that rebuilds its encoding from the rows it is given satisfies the
#' first, since the rows are then the same, and fails the second as soon as
#' \eqn{S} omits a factor level or narrows the range a basis is placed on.
#' The subset is dropped of its unused levels before it is passed, so a
#' plain row subset cannot pass by carrying the original levels along with
#' it.
#'
#' @param term A term specification (an object inheriting from
#'   \code{\link{model_term}}).
#' @param data A data frame.
#' @param verbose Logical; print one line per check.
#'
#' @return Invisibly, a data frame with columns \code{check},
#'   \code{status} (\code{"OK"} or \code{"FAILED"}) and \code{info}.
#'
#' @examples
#' dd <- data.frame(x = 1:6, g = factor(rep(c("a", "b", "c"), 2)))
#' res <- check_term(linpar(~ x + g), dd, verbose = FALSE)
#' all(res$status == "OK")
#'
#' @seealso \code{\link{interpret_formula}}, \code{\link{cens}}
#' @export
check_term <- function(term, data, verbose = TRUE) {
  if (!S7::S7_inherits(term, model_term)) {
    stop("'term' must inherit from 'model_term'.", call. = FALSE)
  }
  if (!is.data.frame(data)) {
    stop("'data' must be a data frame.", call. = FALSE)
  }

  rows <- list()
  add <- function(check, ok, info = "") {
    rows[[length(rows) + 1L]] <<- data.frame(
      check = check,
      status = if (isTRUE(ok)) "OK" else "FAILED",
      info = info, stringsAsFactors = FALSE)
  }
  finish <- function() {
    out <- do.call(rbind, rows)
    if (verbose) {
      for (i in seq_len(nrow(out))) {
        cat(sprintf("  [%s] %s%s\n", out$status[i], out$check[i],
                    if (nzchar(out$info[i])) paste0(" -- ", out$info[i]) else ""))
      }
    }
    invisible(out)
  }

  built <- tryCatch(term_build(term, data), error = function(e) e)
  if (inherits(built, "error")) {
    add("build", FALSE, conditionMessage(built))
    return(finish())
  }
  X <- term_matrix(built)
  add("build",
      .is_block(X) && nrow(X) == nrow(data),
      sprintf("%d x %d block%s", nrow(X), ncol(X),
              if (isS4(X)) ", sparse" else ""))

  cn <- term_coef_names(built)
  add("names", length(cn) == ncol(X) && anyDuplicated(cn) == 0L)
  add("npar", term_npar(built) == ncol(X))

  sm <- tryCatch(term_smooth(built), error = function(e) NA)
  add("smooth", is.logical(sm) && length(sm) == 1L && !is.na(sm),
      if (isTRUE(sm)) "smooth" else if (identical(sm, FALSE)) "non-smooth" else "")

  bare <- function(m) array(as.numeric(m), dim(as.matrix(m)))
  Xr <- tryCatch(term_predict(built, data), error = function(e) e)
  add("reproduce",
      !inherits(Xr, "error") &&
        isTRUE(all.equal(bare(Xr), bare(X), tolerance = 1e-12)),
      if (inherits(Xr, "error")) conditionMessage(Xr) else "")

  idx <- .check_subset_rows(data)
  # droplevels reproduces how new data actually arrive: a factor there
  # carries only the levels it uses, and only the blueprint knows the rest
  nd <- droplevels(data[idx, , drop = FALSE])
  Xs <- tryCatch(term_predict(built, nd), error = function(e) e)
  add("subset",
      !inherits(Xs, "error") && .is_block(Xs) &&
        nrow(Xs) == length(idx) && ncol(Xs) == ncol(X) &&
        isTRUE(all.equal(bare(Xs), bare(X[idx, , drop = FALSE]),
                         tolerance = 1e-12)),
      if (inherits(Xs, "error")) conditionMessage(Xs) else
        sprintf("%d rows", length(idx)))

  finish()
}
