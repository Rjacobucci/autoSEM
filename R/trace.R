#'
#'
#' Calculates the trace of a matrix
#' @param A matrix
#' @keywords trace
#' @export
#' @examples
#' trace()

trace = function(A){
  return(sum(diag(A)))
}
