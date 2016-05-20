#' Summary results from autoSEM.
#'
#' @param object An object from autoSEM.
#' @param ... Other arguments.
#' @export
#' @method summary autoSEM


summary.autoSEM <- function(object, ...){


  ret <- list(call = object$call,
              fit = object$fit,
              solution = object$solution)

  #ret <- object$fit


  class(ret) <- "summary.autoSEM"
  return(ret)
}
