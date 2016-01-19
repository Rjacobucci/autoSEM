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

  if(criterion == "BIC"){
    stop("BIC is not working for aco")
  }

  string = matrix(NA,1,size)
  for(i in 1:size) string[i] = rbinom(1,1,prob=0.5)

  out = list()
  best.sol = 1e10
  first.bic = NA
  best.string <- rep(NA,length(string))

  count=0
  count2 = 1
  pheromone.1 <- rep(1,size)
  pheromone.0 <- rep(1,size)
  samp.wgt <- rep(1,size)
  first.fit2 = NA
  first.fit = fitness(string)
  fit.avg = first.fit

while(count < iters){

  count=count+1

  if(is.na(first.fit)==T & is.na(first.fit2)==T & first.fit < 1e10){
    first.fit2 = 0
  }else{
    first.fit2 = fitness(string)
  }


    current.fit = fitness(string)




  if(criterion =="BIC" & is.na(first.fit) == T){
    first.bic = current.fit
  }

  if(is.numeric(current.fit)==TRUE & current.fit < 1e8){
    if(criterion == "BIC"){
      val.add = round(fit.avg-current.fit,0)
      #val.add = rbinom(1,3,0.5)
    }else if(criterion=="RMSEA"){
      val.add = round(dchisq(current.fit,1)*2,0)
    }else if(criterion=="NCP"){
      #val.add = round(100/(current.fit+1),0)
      val.add = round((fit.avg - current.fit)*3,0)+1
      if(val.add < 0) val.add = 0
    }
  }else{
    val.add=0
  }


    #if(is.numeric(val.add)==F){
     # stop(val.add)
    #}

    if(is.na(val.add)==T){ # val.add > 50
      val.add = 0
    }


  pheromone.1[string==1] <- pheromone.1[string==1] + val.add
  pheromone.0[string==0] <- pheromone.0[string==0] + val.add


  if(is.na(current.fit)==FALSE & current.fit < 1e8){
    if(current.fit < best.sol){
      best.sol = current.fit
      best.string = string
    }
  }


  #samp.wgt = pheromone.1/pheromone.0
  #for(i in 1:size) string[i] = rbinom(1,1,prob=samp.wgt[i])
  for(i in 1:size) {
    pher.samp = c(rep(1,pheromone.1[i]),rep(0,pheromone.0[i]))
    #string[i] = sample(pher.samp,1)
    string[i] = sample(c(sample(pher.samp,1),rbinom(1,1,0.5)),1,prob=c(0.9,.1))
  }


  if(is.na(current.fit)==FALSE & current.fit < 1e8 & val.add != 0){
    fit.avg = (fit.avg*count2 +current.fit)/(count2+1)
    count2 = count2+1
  }


}
out$solution = best.string
out$value = best.sol
out$val.iters = count2
out$samp.wgt = pheromone.1/pheromone.0
out$pheromone.1 = pheromone.1
out$pheromone.0 = pheromone.0
out

}

