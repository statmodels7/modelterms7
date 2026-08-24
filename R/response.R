#' @include term_classes.R
NULL

#' @title S7 Class for a Censored Response
#' @name censored_response
#'
#' @description
#' Holds a response whose observations are not all exact: the values, a lower
#' and an upper bound for each of them, and a status saying which of the four
#' cases the observation is. [cens()] is the constructor to call, and it derives
#' the statuses from the values and the bounds. The raw S7 constructor
#' documented here takes all four vectors and checks only that they agree, so it
#' is the form to use when the statuses are already known.
#'
#' @details
#' # The four statuses and what each asserts
#'
#' Writing \eqn{Y} for the unobserved response, \eqn{L} for `lwr` and \eqn{U}
#' for `upr`, the four values of `status` assert
#'
#' | status | what is known | likelihood contribution |
#' | --- | --- | --- |
#' | `"observed"` | \eqn{Y = y} | \eqn{f(y)} |
#' | `"left"` | \eqn{Y \le L} | \eqn{F(L)} |
#' | `"right"` | \eqn{Y \ge U} | \eqn{1 - F(U)} |
#' | `"interval"` | \eqn{L \le Y \le U} | \eqn{F(U) - F(L)} |
#'
#' with \eqn{f} and \eqn{F} the density and the distribution function of the
#' fitted family. `y` is `NA` for an interval-censored observation, its value
#' being unknown; every other status carries a number.
#'
#' # What reads it today
#'
#' [interpret_formula()] accepts `cens(...)` on the left of a formula and
#' returns the object as the `response` element of its result. Nothing then
#' assembles the four contributions above: `statmodels7::statmod()` stops with
#' a message naming the gap, and no other function in the toolkit reads the
#' class. The pieces exist: `distributions7::distrib_grad_cdf()` and
#' `distrib_hess_cdf()` carry the derivatives of \eqn{F} in the parameters. The
#' assembler that would use them is not written, so the class records what is
#' known about each observation and no more.
#'
#' # What the validator enforces
#'
#' `y`, `lwr`, `upr` and `status` must be the same length; every lower bound
#' must be **strictly** below its upper bound, so `lwr = upr` is rejected;
#' and every status must be one of the four names. A failure throws with the
#' offending rule quoted.
#'
#' @param y A numeric vector of response values, one per observation, `NA`
#'   where the observation is interval-censored.
#' @param lwr,upr Numeric vectors of censoring bounds, of the same length as
#'   `y`. Infinite entries are allowed and mean no censoring on that side. The
#'   validator rejects any pair with `lwr >= upr`.
#' @param status A character vector of the same length as `y`, each entry one
#'   of `"observed"`, `"left"`, `"right"` or `"interval"`. Any other string
#'   throws, naming the first one it meets.
#'
#' @return An S7 object of class `censored_response` with properties `y`,
#'   `lwr`, `upr` (numeric, all of one length) and `status` (character, the
#'   same length). It carries no methods beyond [print()].
#'
#' @seealso [cens()], which derives the statuses instead of taking them;
#'   [interpret_formula()], which accepts it on the left of a formula.
#'
#' @examples
#' # Built through cens(), which is the ordinary route.
#' r <- cens(c(0, 1, 5, NA), lwr = 0, upr = 5)
#' r@status
#' r@y                       # NA survives on the interval-censored row
#' cbind(lwr = r@lwr, upr = r@upr)
#'
#' # The raw constructor takes statuses already known, and checks them.
#' censored_response(y = c(2, 9), lwr = c(-Inf, -Inf), upr = c(Inf, 9),
#'                   status = c("observed", "right"))
#'
#' # Three ways to fail the validator.
#' try(censored_response(y = 1, lwr = 2, upr = 1, status = "observed"))
#' try(censored_response(y = 1, lwr = 0, upr = 2, status = "cut"))
#' try(censored_response(y = c(1, 2), lwr = 0, upr = 2, status = "observed"))
#'
#' @export
censored_response <- S7::new_class(
  name = "censored_response",
  properties = list(
    y = S7::class_numeric,
    lwr = S7::class_numeric,
    upr = S7::class_numeric,
    status = S7::class_character
  ),
  validator = function(self) {
    n <- length(self@y)
    if (length(self@lwr) != n || length(self@upr) != n ||
        length(self@status) != n) {
      return("y, lwr, upr and status must have the same length")
    }
    if (any(self@lwr >= self@upr)) {
      return("every lower bound must be strictly below its upper bound")
    }
    bad <- setdiff(unique(self@status),
                   c("observed", "left", "right", "interval"))
    if (length(bad)) {
      return(sprintf("unknown status '%s'", bad[1L]))
    }
    NULL
  }
)

#' Mark a Response as Censored
#'
#' @description
#' Builds a [censored_response()] from a vector of values and one or two
#' censoring bounds, deriving each observation's status from where its value
#' falls. `cens(y, lwr = 0)` on the left of a model formula declares that any
#' value at or below zero is left-censored there; `cens(y, upr = 8)` does the
#' same at the top; giving both bounds allows either, and an `NA` value between
#' two finite bounds is interval-censored.
#'
#' @details
#' # The rule that assigns a status
#'
#' The bounds are recycled to `length(y)` first, so a scalar applies to every
#' observation and a vector gives one bound per observation. Each entry is then
#' classified by where it sits, with the bounds inclusive:
#'
#' | condition | status | what is known |
#' | --- | --- | --- |
#' | `y <= lwr` | `"left"` | \eqn{Y \le} `lwr` |
#' | `y >= upr` | `"right"` | \eqn{Y \ge} `upr` |
#' | `lwr < y < upr` | `"observed"` | \eqn{Y = y} |
#' | `is.na(y)`, both bounds finite | `"interval"` | \eqn{Y \in [}`lwr`, `upr`\eqn{]} |
#'
#' The tests are applied in that order, so a value at or below `lwr` is left-
#' censored even when it is also at or above `upr`; the validator forbids
#' `lwr >= upr`, so the two cannot both bind. An `NA` value without two finite
#' bounds says nothing about \eqn{Y} at all and throws.
#'
#' At the defaults `lwr = -Inf` and `upr = Inf` no value can reach a bound,
#' every status comes back `"observed"`, and the object carries the same
#' information the bare vector does.
#'
#' # Where it can be used
#'
#' On the left of a formula passed to [interpret_formula()], which returns the
#' object as its `response` element. `statmodels7::statmod()` refuses that
#' response with a message naming the gap: the toolkit marks censoring and does
#' not yet assemble a censored likelihood from it. See [censored_response()]
#' for the four contributions such an assembler would need.
#'
#' @param y A numeric vector of responses, `NA` at an interval-censored
#'   observation. Coerced with [as.numeric()], so an integer vector is accepted
#'   and a factor is not.
#' @param lwr A numeric lower bound, of length 1 or `length(y)`, recycled to
#'   `length(y)`. `-Inf` by default, which is no left censoring. Any other
#'   length throws `"'lwr' must have length 1 or length(y)."`, and an `NA`
#'   entry throws `"'lwr' must not contain NA."`.
#' @param upr A numeric upper bound, on the same terms, `Inf` by default. Every
#'   `upr` must be strictly above its `lwr`; equal bounds throw from the class
#'   validator.
#'
#' @return A [censored_response()] object of length `length(y)`, carrying `y`
#'   unchanged (`NA` included), the recycled `lwr` and `upr`, and the derived
#'   `status`.
#'
#' @seealso [censored_response()] for the class and the four likelihood
#'   contributions; [interpret_formula()] for the formula it goes into;
#'   `distributions7::distrib_grad_cdf()` for the distribution-function
#'   derivatives a censored likelihood is built from.
#'
#' @examples
#' # Left censoring at zero: the first value is at the bound, so it is
#' # censored there and the other two are exact.
#' r <- cens(c(0, 0.7, 2.4), lwr = 0)
#' r@status
#'
#' # One of each status. The NA is interval-censored between the bounds.
#' r2 <- cens(c(0, 1, 5, NA), lwr = 0, upr = 5)
#' r2
#' data.frame(y = r2@y, lwr = r2@lwr, upr = r2@upr, status = r2@status)
#'
#' # Per-observation bounds: only the rows whose own bound binds are censored.
#' cens(c(1, 2, 3), lwr = c(-Inf, 2, -Inf), upr = c(Inf, Inf, 3))@status
#'
#' # With no bounds given, every observation is exact.
#' all(cens(rnorm(20))@status == "observed")
#'
#' # An NA with no finite pair of bounds carries no information.
#' try(cens(c(1, NA)))
#'
#' # It goes on the left of a formula.
#' d <- data.frame(t = c(1, 5, 9, 2), x = c(1, 2, 3, 4))
#' names(interpret_formula(cens(t, upr = 8) ~ x, d))
#'
#' @export
cens <- function(y, lwr = -Inf, upr = Inf) {
  y <- as.numeric(y)
  n <- length(y)
  recycle <- function(v, nm) {
    v <- as.numeric(v)
    if (length(v) == 1L) v <- rep(v, n)
    if (length(v) != n) {
      stop(sprintf("'%s' must have length 1 or length(y).", nm),
           call. = FALSE)
    }
    if (anyNA(v)) {
      stop(sprintf("'%s' must not contain NA.", nm), call. = FALSE)
    }
    v
  }
  lwr <- recycle(lwr, "lwr")
  upr <- recycle(upr, "upr")

  status <- rep("observed", n)
  obs <- !is.na(y)
  status[obs & y <= lwr] <- "left"
  status[obs & y >= upr] <- "right"
  na <- !obs
  if (any(na & !(is.finite(lwr) & is.finite(upr)))) {
    stop("an NA response is interval-censored and needs finite 'lwr' and 'upr'.",
         call. = FALSE)
  }
  status[na] <- "interval"

  censored_response(y = y, lwr = lwr, upr = upr, status = status)
}

#' @title Print a Censored Response
#' @name print.censored_response
#'
#' @description
#' Prints one line giving the number of observations and how many carry each
#' status, in the fixed order `observed`, `left`, `right`, `interval`. A status
#' no observation has is left out, so a response with no censoring at all
#' prints as `<censored_response> 20 observations: 20 observed`. The values and
#' the bounds are not shown; read them from the `y`, `lwr` and `upr`
#' properties.
#'
#' @param x A [censored_response()] object.
#' @param ... Unused, and accepted so that the signature matches [print()]'s.
#'
#' @return `x`, invisibly. Called for the line it writes to the console.
#'
#' @seealso [cens()], which builds the object and assigns the statuses.
#'
#' @examples
#' cens(c(0, 1, 5, NA), lwr = 0, upr = 5)
#'
#' # Only the statuses present are listed.
#' cens(c(1, 2, 3))
#'
#' @keywords internal
S7::method(print, censored_response) <- function(x, ...) {
  n <- length(x@y)
  counts <- table(factor(x@status,
                         levels = c("observed", "left", "right", "interval")))
  shown <- counts[counts > 0L]
  cat(sprintf("<censored_response> %d observations: %s\n", n,
              paste(sprintf("%d %s", shown, names(shown)), collapse = ", ")))
  invisible(x)
}
