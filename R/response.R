#' @include term_classes.R
NULL

#' @title S7 Class for a Censored Response
#' @name censored_response
#'
#' @description
#' The response object \code{\link{cens}} constructs: the observed values,
#' the per-observation censoring bounds, and the status each observation
#' carries (\code{"observed"}, \code{"left"}, \code{"right"} or
#' \code{"interval"}). The likelihood assembler of the model layer consumes
#' it, contributing a density where the observation is exact and a
#' difference of distribution functions where it is censored.
#'
#' @param y The numeric response values (\code{NA} for an
#'   interval-censored observation).
#' @param lwr,upr The numeric censoring bounds, one value per observation.
#' @param status The character vector of per-observation statuses.
#'
#' @return An object of class \code{censored_response}.
#'
#' @seealso \code{\link{cens}}
#' @examples
#' S7::S7_inherits(cens(c(0, 1.2), lwr = 0), censored_response)
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

#' Censored Response Constructor
#'
#' @description
#' Marks a response as censored, for the left-hand side of a model formula:
#' \code{cens(y, lwr = 0)} in a formula declares that values at or below
#' the bound are left-censored there.
#'
#' @details
#' The bounds are recycled to the length of \code{y}, so a scalar bound
#' applies to every observation and a vector gives per-observation bounds.
#' The status of each observation follows from the values: an observation
#' with \code{y <= lwr} is left-censored (all that is known is
#' \eqn{Y \le lwr}), one with \code{y >= upr} is right-censored, one with
#' \code{y} strictly inside the bounds is observed exactly, and one with
#' \code{y = NA} and both bounds finite is interval-censored
#' (\eqn{Y \in [lwr, upr]}). An \code{NA} value without two finite bounds
#' carries no information and is rejected.
#'
#' @param y A numeric vector; \code{NA} for interval-censored
#'   observations.
#' @param lwr A numeric vector of lower bounds, length 1 or
#'   \code{length(y)}. Defaults to \code{-Inf} (no left censoring).
#' @param upr A numeric vector of upper bounds, length 1 or
#'   \code{length(y)}. Defaults to \code{Inf} (no right censoring).
#'
#' @return An object of class \code{\link{censored_response}}.
#'
#' @examples
#' r <- cens(c(0, 0.7, 2.4), lwr = 0)
#' r@status
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

S7::method(print, censored_response) <- function(x, ...) {
  n <- length(x@y)
  counts <- table(factor(x@status,
                         levels = c("observed", "left", "right", "interval")))
  shown <- counts[counts > 0L]
  cat(sprintf("<censored_response> %d observations: %s\n", n,
              paste(sprintf("%d %s", shown, names(shown)), collapse = ", ")))
  invisible(x)
}
