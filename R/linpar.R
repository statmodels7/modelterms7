#' @include term_classes.R generics.R
NULL

#' Unpenalized Parametric Term
#'
#' @description
#' Specifies an unpenalized parametric block: the model matrix of a one-sided
#' formula, with the usual [stats::model.matrix()] conventions for factors,
#' contrasts, interactions and the intercept. It is the plainest term in the
#' package, and the one [interpret_formula()] collects a formula's bare
#' covariates into.
#'
#' The term's contribution to the predictor is \eqn{X\beta}, with no penalty,
#' so every column costs one degree of freedom and [edf()] returns the column
#' count exactly.
#'
#' @details
#' # The block, and what a build records
#'
#' Building runs [stats::model.frame()] with `na.action = na.pass` and
#' `drop.unused.levels = FALSE`, so a row with a missing covariate keeps its
#' place and the block stays aligned with the response, and a factor level
#' present in the data but used by no row still gets its column.
#'
#' The blueprint records the terms object, the factor levels
#' ([stats::.getXlevels()]), the contrasts actually used and the storage
#' settled on. [term_predict()] reapplies all four, so a factor column at new
#' rows is encoded against the levels seen at build time. A level the blueprint
#' does not know is rejected by [stats::model.frame()] with
#' `"factor g has new levels zz"`.
#'
#' # Several parametric blocks
#'
#' `y ~ x1 + x2` and `y ~ linpar(~ x1 + x2)` produce the same block, so the
#' explicit constructor is for callers wanting more than one parametric block
#' with distinct labels, or wanting to set `sparse` or `contrasts` on it.
#' Arguments for the block [interpret_formula()] builds implicitly go through
#' that function's own `linpar` argument.
#'
#' @section Sparse storage:
#' `sparse = TRUE` builds through [Matrix::sparse.model.matrix()], which builds
#' the block sparse. Building a dense matrix and compressing it would cost the
#' memory the choice exists to avoid. Measured at 20000 rows and a factor of
#' 1000 levels, the two routes give identical numbers at 0.007 s and 1.8 MB
#' against 0.164 s and 161.3 MB; a design of 20000 by 19014 that would be 3.0
#' GB dense builds sparse in 0.3 s and 3.1 MB.
#'
#' It pays where the formula carries a **factor of many levels**, whose
#' indicator columns hold one non-zero per row. On numeric covariates the block
#' is dense whatever is asked for, and the sparse storage then costs more than
#' it saves.
#'
#' `sparse = NULL`, the default, settles it at build from the size of the
#' indicator part: `n * ncol_ind > 1e5`, where `ncol_ind` counts the columns
#' coming from factors ([.indicator_cols()]). Measured over fifteen
#' combinations, building the block and forming its crossproduct, the two
#' routes cross between \eqn{10^5} and \eqn{3 \times 10^5} cells and the sparse
#' route then wins by orders: at 20000 rows it is 1.4 times faster at 15
#' levels, 14 times at 60 and 445 times at 400.
#'
#' The settled storage is part of the blueprint, so [term_predict()] returns
#' the same kind of block at new rows.
#'
#' @param formula A one-sided formula, `~ x1 + x2`. A two-sided formula throws
#'   `"'formula' must be one-sided, e.g. ~ x1 + x2."`, and anything that is not
#'   a formula throws. Its environment is kept and used for symbols the data do
#'   not carry. A formula with no columns at all, `~ 0`, fails in the class
#'   validator at build time.
#' @param label A single character string, `""` by default. When non-empty it
#'   is prefixed to every coefficient name as `label.name`, and it is the title
#'   [plot()] uses. Anything that is not one string throws.
#' @param sparse `TRUE` to build a `dgCMatrix`, `FALSE` for a base matrix, or
#'   `NULL`, the default, to settle it from the design by the rule above.
#'   Anything else throws
#'   `"'sparse' in 'linpar' must be TRUE, FALSE, or NULL to settle it from the design."`.
#' @param contrasts A named list of contrasts for the formula's factors, of the
#'   kind [stats::model.matrix()]'s `contrasts.arg` takes, or `NULL`, the
#'   default, for the session's `options("contrasts")`. Anything that is not a
#'   list throws.
#'
#' @return An unbuilt [LinparTerm()]: a specification, with `X`, `coef_names`
#'   and `blueprint` empty until [term_build()] fills them, and `penalty`
#'   `NULL` permanently.
#'
#' @seealso [ridge()], [lasso()], [scad()], [mcp()] and [enet()] for the
#'   penalized blocks; [s()] and [random()] for the penalized structures;
#'   [interpret_formula()], which builds one of these implicitly;
#'   [term_build()] and [term_predict()].
#'
#' @examples
#' dd <- data.frame(x = 1:8, g = factor(rep(c("a", "b", "c", "d"), 2)))
#'
#' built <- term_build(linpar(~ x + g), dd)
#' term_matrix(built)
#' term_coef_names(built)
#'
#' # Unpenalized: every column costs one degree of freedom.
#' c(npar = term_npar(built), edf = edf(built))
#'
#' # A label prefixes the names, which is how two blocks stay apart.
#' term_coef_names(term_build(linpar(~ x, label = "lin"), dd))
#'
#' # Contrasts are recorded at build and reapplied at prediction.
#' bc <- term_build(linpar(~ g, contrasts = list(g = "contr.sum")), dd)
#' term_coef_names(bc)
#'
#' # A missing covariate keeps its row, so the block stays aligned.
#' term_matrix(term_build(linpar(~ x),
#'                        data.frame(x = c(1, NA, 3))))
#'
#' # Storage is settled from the design: a small factor stays dense.
#' small <- data.frame(g = factor(rep(1:3, 10)))
#' class(term_matrix(term_build(linpar(~ g), small)))
#'
#' # Two hundred levels over two thousand rows is 4e5 cells, so sparse.
#' set.seed(1)
#' big <- data.frame(g = factor(sample(1:200, 2000, TRUE)))
#' class(term_matrix(term_build(linpar(~ g), big)))
#'
#'
#' # Fitted. The data are simulated from a known truth, so the
#' # estimates below can be read against it.
#' if (requireNamespace("statmodels7", quietly = TRUE)) {
#'   set.seed(1)
#'   fd <- data.frame(x = runif(200, -2, 2),
#'                    g = factor(sample(c("a", "b"), 200, TRUE)))
#'   fd$y <- 1 + 2 * fd$x + 0.5 * (fd$g == "b") + rnorm(200, sd = 0.4)
#'   ft <- statmodels7::statmod(y ~ x + g,
#'                              distributions7::gaussian1_distrib(), fd)
#'   # truth: intercept 1, slope 2, and 0.5 for the second level
#'   round(coef(ft)$mu, 2)
#' }
#' @export
linpar <- function(formula, label = "", sparse = NULL, contrasts = NULL) {
  if (!inherits(formula, "formula")) {
    stop("'formula' must be a formula.", call. = FALSE)
  }
  if (length(formula) != 2L) {
    stop("'formula' must be one-sided, e.g. ~ x1 + x2.", call. = FALSE)
  }
  if (!is.character(label) || length(label) != 1L || is.na(label)) {
    stop("'label' must be a single character string.", call. = FALSE)
  }
  o <- .check_design_opts(sparse, contrasts, "linpar")
  LinparTerm(label = label, formula = formula, sparse = o$sparse,
             contrasts = o$contrasts,
             X = NULL, coef_names = character(0),
             blueprint = list(), penalty = NULL)
}

#' The Model Matrix of a Formula, in Either Storage
#'
#' @description
#' Runs [stats::model.matrix()] or [Matrix::sparse.model.matrix()] on the same
#' terms object and model frame, strips the `assign` and `contrasts`
#' attributes, and returns the block beside the contrasts that were used.
#'
#' @details
#' The two routes give identical numbers; what differs is the storage and the
#' cost of producing it. The sparse route builds the matrix sparse and never
#' forms a dense intermediate: at 20000 rows and a factor of 1000 levels,
#' 0.007 s and 1.8 MB against 0.164 s and 161.3 MB.
#'
#' The contrasts come back beside the block instead of staying on it: the block
#' is what a consumer reads, and the contrasts are what the blueprint records.
#'
#' @param tt A terms object, from [stats::model.frame()] or
#'   [stats::delete.response()].
#' @param mf The model frame `tt` was built from, or one built against the same
#'   levels.
#' @param contrasts A named list of contrasts, or `NULL` for the session's.
#' @param sparse `TRUE` for a `dgCMatrix`, `FALSE` for a base matrix.
#'
#' @return A list of two: `X`, the block, and `contrasts`, the named list
#'   [stats::model.matrix()] recorded, which is `NULL` when no factor was
#'   coded.
#'
#' @seealso [linpar()], whose build and prediction both call this;
#'   [.resolve_sparse()] for the choice of storage.
#'
#' @keywords internal
.design_matrix <- function(tt, mf, contrasts = NULL, sparse = FALSE) {
  X <- if (isTRUE(sparse)) {
    Matrix::sparse.model.matrix(tt, mf, contrasts.arg = contrasts)
  } else {
    stats::model.matrix(tt, mf, contrasts.arg = contrasts)
  }
  ctr <- attr(X, "contrasts")
  attr(X, "assign") <- NULL
  attr(X, "contrasts") <- NULL
  list(X = X, contrasts = ctr)
}


#' @title Check a Term's Storage and Contrasts
#'
#' @description
#' Validates the `sparse` and `contrasts` arguments the design-building
#' constructors share, and normalizes `contrasts = NULL` to an empty list so
#' that the class property is always a list. Called from the constructor, so a
#' mistake is reported where it was written.
#'
#' @param sparse What the constructor was given: `TRUE`, `FALSE`, or `NULL` for
#'   the storage to be settled at build. Anything else throws, the message
#'   naming `what`.
#' @param contrasts What the constructor was given: a named list, or `NULL`.
#'   Anything else throws. The names are not checked against the formula's
#'   factors; [stats::model.matrix()] does that at build.
#' @param what The constructor's name, used in the two messages.
#'
#' @return A list of two: `sparse` unchanged, and `contrasts` as given or an
#'   empty list.
#'
#' @seealso [linpar()] and [penalized_terms()], the constructors that call it.
#'
#' @keywords internal
.check_design_opts <- function(sparse, contrasts, what = "this term") {
  if (!is.null(sparse) &&
      (!is.logical(sparse) || length(sparse) != 1L || is.na(sparse))) {
    stop(sprintf(paste0("'sparse' in '%s' must be TRUE, FALSE, or NULL to",
                        " settle it from the design."), what), call. = FALSE)
  }
  if (!is.null(contrasts) && !is.list(contrasts)) {
    stop(sprintf(paste0("'contrasts' in '%s' must be a named list, one entry",
                        " per factor, or NULL."), what), call. = FALSE)
  }
  list(sparse = sparse,
       contrasts = if (is.null(contrasts)) list() else contrasts)
}


#' @title How Many Columns of a Model Matrix Come From Factors
#'
#' @description
#' Counts the columns a terms object contributes through indicator variables,
#' read from the model frame without building the matrix. It is the second
#' factor of the product [.resolve_sparse()] compares against its threshold.
#'
#' @details
#' A term contributes the product of its variables' level counts, a numeric
#' variable counting one, and a term carrying no factor at all is skipped: its
#' columns are dense whatever the storage, so they say nothing about the
#' choice. Character and logical columns are counted by their number of
#' distinct values, since [stats::model.matrix()] will code them as factors.
#'
#' The count is an **upper bound**: contrasts drop one level per factor, and an
#' interaction of two factors is counted at the product of their full level
#' counts. That is the right side to err on, the number being read against a
#' threshold below which the sparse route loses little and above which it wins
#' by orders.
#'
#' @param tt A terms object.
#' @param mf The model frame it was built from.
#'
#' @return A single number, `0` when no term carries a factor.
#'
#' @seealso [.resolve_sparse()], which reads it; [linpar()].
#'
#' @keywords internal
.indicator_cols <- function(tt, mf) {
  fac <- attr(tt, "factors")
  if (is.null(fac) || !length(fac) || !ncol(fac)) return(0)
  nlev <- vapply(rownames(fac), function(v) {
    x <- mf[[v]]
    if (is.null(x)) 1
    else if (is.factor(x)) as.numeric(nlevels(x))
    else if (is.character(x) || is.logical(x)) as.numeric(length(unique(x)))
    else 1
  }, numeric(1))
  tot <- 0
  for (j in seq_len(ncol(fac))) {
    inn <- fac[, j] > 0
    if (!any(inn) || !any(nlev[inn] > 1)) next
    tot <- tot + prod(nlev[inn])
  }
  tot
}


#' @title Settle Whether a Block Is Built Sparse
#'
#' @description
#' Passes an explicit `TRUE` or `FALSE` through, and where the caller left
#' `NULL` decides from the size of the indicator part: sparse when
#' `n * ncol_ind` exceeds \eqn{10^5}.
#'
#' @details
#' The dense indicator part holds `n * ncol_ind` cells where the sparse one
#' holds one non-zero per row, so that product is what separates the two
#' routes, and the threshold is read off it rather than off a count of levels.
#'
#' Measured over fifteen combinations of sample size and level count, building
#' the block and forming its crossproduct, the routes cross between
#' \eqn{10^5} and \eqn{3 \times 10^5} cells. Below that the sparse route loses
#' a little (0.6 to 0.9 times the dense route); above it the gap opens quickly:
#' at 20000 rows the ratios are 1.4 at 15 levels, 12.3 at 25, 14.4 at 60,
#' 139 at 200 and 445 at 400.
#'
#' A design carrying no factor has no indicator part, so the product is zero
#' and the block is built dense, which is the right answer for a purely
#' continuous design.
#'
#' @param sparse `TRUE`, `FALSE`, or `NULL` to decide here.
#' @param n The number of rows.
#' @param ncol_ind The columns coming from indicators, from
#'   [.indicator_cols()].
#'
#' @return A single logical, never `NA`. A non-finite `n` or `ncol_ind` gives
#'   `FALSE`.
#'
#' @seealso [.indicator_cols()] for the count it reads, [linpar()] for the
#'   argument it settles.
#'
#' @keywords internal
.resolve_sparse <- function(sparse, n, ncol_ind) {
  if (!is.null(sparse)) return(isTRUE(sparse))
  isTRUE(is.finite(n) && is.finite(ncol_ind) && n * ncol_ind > 1e5)
}


#' @title Build a Parametric Block
#' @name term_build.LinparTerm
#'
#' @description
#' Builds the model matrix of a [linpar()] term's formula against `data`,
#' prefixes the coefficient names with the term's label, and records the
#' blueprint [term_predict()] will reapply: the terms object with the response
#' deleted, the factor levels, the contrasts used, and the storage settled on.
#'
#' @details
#' # The model frame
#'
#' [stats::model.frame()] is called with `na.action = na.pass`, so a row with a
#' missing covariate keeps its place and the block stays aligned with the
#' response, and with `drop.unused.levels = FALSE`, so a factor level present
#' in the data but used by no row still gets a column. Both choices are about
#' alignment: a block whose rows have been silently dropped no longer matches
#' the response it is fitted against.
#'
#' # The storage
#'
#' Where the constructor left `sparse = NULL`, [.resolve_sparse()] settles it
#' from `n` times the indicator column count. The value settled on goes into
#' `blueprint$sparse` and the `sparse` property is left as the caller wrote it,
#' because new rows may be any number and a prediction deciding again could
#' build a block of a different kind from the fitted one.
#'
#' @param term An unbuilt or built [LinparTerm()].
#' @param data A data frame carrying every variable the formula names.
#' @param ... Unused.
#'
#' @return The term with `X`, `coef_names` and `blueprint` filled, so that
#'   [term_is_built()] is `TRUE`. The block has `nrow(data)` rows, and its
#'   column names are the coefficient names.
#'
#' @seealso [linpar()], [term_predict.LinparTerm()], [.resolve_sparse()].
#'
#' @examples
#' dd <- data.frame(x = 1:8, g = factor(rep(c("a", "b", "c", "d"), 2)))
#' b <- term_build(linpar(~ x + g), dd)
#' term_matrix(b)
#' names(b@blueprint)
#'
#' # A missing value keeps its row.
#' term_matrix(term_build(linpar(~ x), data.frame(x = c(1, NA, 3))))
#'
#' @keywords internal
S7::method(term_build, LinparTerm) <- function(term, data, ...) {
  mf <- stats::model.frame(term@formula, data,
                           na.action = stats::na.pass,
                           drop.unused.levels = FALSE)
  tt <- attr(mf, "terms")
  # the block is a plain numeric matrix, or a dgCMatrix where the caller
  # asked for one or where the design is large enough that one pays: the
  # model.matrix bookkeeping lives in the blueprint, not on the result
  sp <- .resolve_sparse(term@sparse, nrow(mf), .indicator_cols(tt, mf))
  b <- .design_matrix(tt, mf,
                      if (length(term@contrasts)) term@contrasts else NULL,
                      sp)
  X <- b$X
  cn <- colnames(X)
  if (nzchar(term@label)) cn <- paste(term@label, cn, sep = ".")
  colnames(X) <- cn
  term@X <- X
  term@coef_names <- cn
  term@blueprint <- list(
    terms = stats::delete.response(tt),
    xlev = stats::.getXlevels(tt, mf),
    contrasts = b$contrasts,
    # the SETTLED storage, not what the constructor was given: new data may
    # carry any number of rows, and a prediction that decided again could
    # build a block of a different kind from the one that was fitted
    sparse = sp
  )
  term
}

#' @title A Parametric Block at New Rows
#' @name term_predict.LinparTerm
#'
#' @description
#' Rebuilds the model frame at `newdata` against the levels recorded at build
#' time, then the model matrix with the recorded contrasts and in the recorded
#' storage, and labels the columns with the term's coefficient names. Nothing
#' is re-derived from the new rows.
#'
#' @details
#' The levels come from `blueprint$xlev`, so a factor in `newdata` need carry
#' only the levels its own rows use and still gets the full set of columns. A
#' level the blueprint does not know is rejected by [stats::model.frame()] with
#' `"factor g has new levels zz"`, which is the right answer: a coefficient was
#' never fitted for it.
#'
#' The storage comes from `blueprint$sparse` rather than being decided again,
#' so a prediction does not spend at new rows what the build was careful not
#' to. `na.action = na.pass` again keeps every row.
#'
#' @param term A built [LinparTerm()]. An unbuilt one throws
#'   `"the term has not been built; call term_build(term, data) first."`.
#' @param newdata A data frame carrying every variable the formula names.
#' @param ... Unused.
#'
#' @return A block of `nrow(newdata)` rows and [term_npar()] columns, in the
#'   storage the build settled on, with the term's coefficient names as column
#'   names.
#'
#' @seealso [term_predict()] for the generic and the identity it satisfies,
#'   [term_build.LinparTerm()] for what recorded the blueprint.
#'
#' @examples
#' dd <- data.frame(x = 1:8, g = factor(rep(c("a", "b", "c", "d"), 2)))
#' b <- term_build(linpar(~ x + g), dd)
#'
#' # On the fitting data it returns the block itself.
#' all.equal(term_predict(b, dd), term_matrix(b))
#'
#' # A subset that drops two levels keeps all five columns.
#' nd <- droplevels(dd[dd$g %in% c("a", "b"), ])
#' levels(nd$g)
#' dim(term_predict(b, nd))
#'
#' # A level the fit never saw is refused.
#' bad <- dd
#' levels(bad$g) <- c("a", "b", "c", "zz")
#' try(term_predict(b, bad))
#'
#' @keywords internal
S7::method(term_predict, LinparTerm) <- function(term, newdata, ...) {
  .assert_built(term)
  bp <- term@blueprint
  mf <- stats::model.frame(bp$terms, newdata,
                           na.action = stats::na.pass,
                           xlev = bp$xlev)
  # the STORAGE is part of the blueprint: a prediction that densified would
  # spend at new data what the build was careful not to
  X <- .design_matrix(bp$terms, mf, bp$contrasts, isTRUE(bp$sparse))$X
  colnames(X) <- term@coef_names
  X
}
