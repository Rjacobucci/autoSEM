#' Implementation of ant colony optimization by ross jacobucci
#'
#'
#' @param size Tabu list size. This is the number of variables * factors.
#' @param iters Maximum number of iterations to allow.
#' @param fitness Fitness function to use.
#' @param criterion Criterion to use.
#' @keywords aco_rj
#' @export
#' @examples
#' \dontrun{
#' aco_rj()
#'}



aco_rj <- function(size, iters,fitness,criterion){

  string = matrix(NA,1,size)
  for(i in 1:size) string[i] = rbinom(1,1,prob=0.5)

  out = list()
  best.sol = 999999999
  first.bic = NA
  best.string <- rep(NA,length(string))
  count=0
  pheromone.1 <- rep(1,size)
  pheromone.0 <- rep(1,size)
  samp.wgt <- rep(1,size)

while(count < iters){

  count=count+1


  current.fit = fitness(string)

  if(criterion =="BIC" & is.na(first.bic) == T){
    first.bic = current.fit
  }

  if(is.numeric(current.fit)==TRUE){
    if(criterion == "BIC"){
      val.add = round(pnorm(current.fit-first.bic)*3,0)
    }else if(criterion=="RMSEA"){
      val.add = round(dchisq(current.fit,1),0)
    }else if(criterion=="NCP"){
      val.add = round(100/(current.fit+1),0)
    }
  }else{
    val.add=0
  }


  pheromone.1[string==1] <- pheromone.1[string==1] + val.add
  pheromone.0[string==0] <- pheromone.0[string==0] + val.add


  if(current.fit < best.sol){
    best.sol = current.fit
    best.string = string
  }

  samp.wgt = pheromone.1/pheromone.0
  for(i in 1:size) string[i] = rbinom(1,1,prob=pnorm(samp.wgt[i]))


}
out$solution = best.string
out$value = best.sol
out$samp.wgt = samp.wgt
out

}

