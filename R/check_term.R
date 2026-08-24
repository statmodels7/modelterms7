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

#' Structural Checks on a Model Term
#'
#' @description
#' Builds `term` against `data` and runs six checks on the result: that the
#' build succeeds and returns a two-dimensional numeric block with one row per
#' observation, that the coefficient names are unique and as numerous as the
#' columns, that [term_npar()] agrees with that count, that [term_smooth()]
#' answers with a single non-missing logical, that [term_predict()] on the same
#' data reproduces the block, and that [term_predict()] on a subset of rows
#' returns the corresponding rows of it. One row of the result per check,
#' printed as it goes and returned invisibly.
#'
#' The last check is the one worth running. A term is supposed to record its
#' encoding at build time and reapply it; a term that re-derives the encoding
#' from whatever rows it is handed passes every other check and fails this one.
#'
#' @details
#' # The two identities
#'
#' Write \eqn{X} for `term_matrix(term_build(term, data))` and \eqn{S} for a
#' subset of the row indices. The last two checks are
#'
#' \deqn{\mathrm{predict}(\mathrm{term}, \mathrm{data}) = X,
#'   \qquad
#'   \mathrm{predict}(\mathrm{term}, \mathrm{data}[S, ]) = X[S, ],}
#'
#' both compared by [all.equal()] at a relative tolerance of `1e-12`. The second
#' does not follow from the first: a term that rebuilds its encoding satisfies
#' the first, the rows being the same ones, and fails the second as soon as
#' \eqn{S} omits a factor level or narrows the range a basis is placed on.
#'
#' # Which rows the subset takes
#'
#' The subset is chosen to make that failure reachable. If any column of `data`
#' is a factor with two or more levels, every row of its last level is dropped,
#' so a rebuilt encoding is one column short. Failing that, the first
#' `nrow(data) %/% 2` rows are taken, which still narrows the range of a numeric
#' covariate.
#'
#' [droplevels()] is applied to the subset before it is passed on, because that
#' is how new data really reach [term_predict()]: a factor column there carries
#' only the levels its own rows use, and the rest are known to the blueprint
#' alone. Without the call a plain row subset carries the original level set
#' along with it and the check cannot fail.
#'
#' # What is not checked
#'
#' The battery covers the design block. It reads neither the penalty a
#' penalized term attaches nor the hyperparameters it declares, so
#' [term_penalties()] and [term_hyper()] are never called.
#'
#' `check_term()` applies to an additive term. A structural term ([gas()],
#' [regime()]) contributes no design columns and registers no [term_matrix()]
#' method, so the call stops with S7's method-not-found error at the first
#' check.
#'
#' @param term A term specification: any object inheriting from [model_term()],
#'   built or unbuilt. Anything else throws
#'   `"'term' must inherit from 'model_term'."`. An already-built term is
#'   accepted and is rebuilt against `data`.
#' @param data A data frame carrying every variable the term names. Anything
#'   else throws `"'data' must be a data frame."`. Any number of rows is
#'   accepted; with one row the subset check compares the block against itself
#'   and is uninformative.
#' @param verbose Print one line per check as it completes, `TRUE` by default,
#'   in the form `  [OK] subset -- 4 rows`. The rows are returned either way, so
#'   `verbose = FALSE` is the form to use inside a test.
#'
#' @return Invisibly, a data frame of `check`, `status` and `info`, all
#'   character, one row per check:
#'   \describe{
#'     \item{`check`}{`"build"`, `"names"`, `"npar"`, `"smooth"`,
#'       `"reproduce"` and `"subset"`, in that order.}
#'     \item{`status`}{`"OK"` or `"FAILED"`.}
#'     \item{`info`}{The block's dimensions for `build`, with `", sparse"`
#'       appended when it is an S4 matrix; `"smooth"` or `"non-smooth"` for
#'       `smooth`; the number of rows kept for `subset`; the condition message
#'       where a check threw; and `""` otherwise.}
#'   }
#'   A build that throws gives a single row, `check = "build"` and
#'   `status = "FAILED"`, with the message in `info`; no later check is
#'   attempted.
#'
#' @seealso [term_build()] and [term_predict()], the two generics the checks
#'   compare; [term_matrix()] for the block itself; [interpret_formula()] for
#'   reading a whole formula into terms.
#'
#' @examples
#' dd <- data.frame(x = 1:6, g = factor(rep(c("a", "b", "c"), 2)))
#'
#' # A parametric block passes all six.
#' res <- check_term(linpar(~ x + g), dd)
#' all(res$status == "OK")
#'
#' # What the subset check is testing. The chosen rows drop level "c", and
#' # the reapplied block keeps the four columns the blueprint recorded.
#' b  <- term_build(linpar(~ x + g), dd)
#' nd <- droplevels(dd[c(1, 2, 4, 5), ])
#' levels(nd$g)
#' dim(term_predict(b, nd))              # 4 x 4, the blueprint's levels
#' dim(model.matrix(~ x + g, nd))        # 4 x 3, what a rebuild gives
#'
#' # The same failure for a basis: rebuilding places the knots on the
#' # narrower range, so the columns are different functions of x.
#' d2  <- data.frame(x = seq(0, 1, length.out = 40))
#' bs  <- term_build(s(x, k = 6), d2)
#' X   <- term_matrix(bs)
#' sub <- 1:20
#' max(abs(term_predict(bs, d2[sub, , drop = FALSE]) - X[sub, ]))
#' max(abs(term_matrix(term_build(s(x, k = 6), d2[sub, , drop = FALSE])) -
#'         X[sub, ]))
#'
#' # A sparse block is accepted, and the info column says it is sparse.
#' check_term(random(~ 1 | g), dd)
#'
#' # A build that throws gives one row and stops there.
#' print(check_term(linpar(~ x + nonexistent), dd, verbose = FALSE))
#'
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
