#'
#'
#' Calculates the fit functions
#' @param df degrees of freedom
#' @param N sample size
#' @param chisq
#' @param ncp non-centrality parameter
#' @keywords fit chisq ncp rmsea
#' @export
#' @examples
#' d()
#' rmsea()


#null.df = function(k) k*(k-1)/2

#observations = function(k) k*(k+1)/2

# NCP
d = function(chisq,df,N) max(0,(chisq -df)/(N-1))

# RMSEA
rmsea = function(ncp,df) sqrt(ncp/df)

