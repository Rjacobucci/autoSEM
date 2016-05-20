
aco <- function(size, iters,fitness,criterion,min.improve,seed){

 # if(criterion == "BIC" | criterion == "BIC2" | criterion == "AIC"){
  #  stop("BIC is not working for aco")
  #}

  string = matrix(NA,1,size)
  set.seed(seed)
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
  iter.no.improve = 0

while(count < iters){

  count=count+1


    current.fit = fitness(string)




 # if(criterion =="BIC" | criterion =="BIC2" |
    #    criterion =="AIC" & is.na(first.bic) == T){
    #   first.bic = current.fit
    # }

  if(is.numeric(current.fit)==TRUE & current.fit < 1e8){
    if(criterion == "BIC" | criterion == "AIC" | criterion == "BIC2"){
      if(count == 1){
        val.add = 4
      }else{
        val.add = round(best.sol-current.fit,0)
      }
      if(val.add < 0 | is.na(val.add)==TRUE) val.add = 0
      val.add = min(val.add,50)
      #val.add = rbinom(1,3,0.5)
    }else if(criterion=="RMSEA"){
      val.add = round(dchisq(current.fit,1)*2,0)
      val.add = min(val.add,50)
    }else if(criterion=="NCP"){
      #val.add = round(100/(current.fit+1),0)
      val.add = round((best.sol - current.fit)*3,0)+1
      if(val.add < 0) val.add = 0
      val.add = min(val.add,50)
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
      iter.no.improve = 0
    }else{
      iter.no.improve = iter.no.improve + 1
    }
  }

    if(iter.no.improve == min.improve){
      break
    }

  if(any(pheromone.1 <1) | any(pheromone.0 < 1)){
    stop(c(max(c(pheromone.1,pheromone.0)),val.add))
  }

  #samp.wgt = pheromone.1/pheromone.0
  #for(i in 1:size) string[i] = rbinom(1,1,prob=samp.wgt[i])
  for(i in 1:size) {
    pher.samp = c(rep(1,pheromone.1[i]),rep(0,pheromone.0[i]))
    #string[i] = sample(pher.samp,1)
    string[i] = sample(c(sample(pher.samp,1),rbinom(1,1,0.5)),1,prob=c(0.85,.1))
  }


  if(is.na(current.fit)==FALSE & current.fit < 1e8){ # & val.add != 0){
#    fit.avg = (fit.avg*count2 +current.fit)/(count2+1)
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

