# Author : Nicolas S. Müller (nicolas.muller@unige.ch)
#
# This is an example of how to run the tabusearch

TABU <- read.csv(file.choose(), sep=",", header=TRUE)
 tabusearch("v1", names(TABU)[2:18], TABU$y, TABU, 10)


#### ARGUMENTS ####
# K : starting model
# vars : vector with the variable names
# dep : dependent variable
# dataset : name of the dataset
# tabusize : size of the tabu list
# type : type of model ("lm" for a linear model, "cox" for a proportional risk model)
# test : measure used to compare models. At the moment, only "Ftest" is supported for "lm", and "AIC" or "wald" for "cox"
###################
#### OUTPUT #####
# The output is a list with 3 elements. The "K" element is a vector with the names of the retained variables, the "measure" one is the value of the test for the best model,
# and the "bestmodel" is the "lm" output object for the best model found.
#################

tabusearch <- function(K, vars, dep, dataset, tabusize, type="lm", test="Ftest") {
 

  tabulist <- vector("list", length(vars))
  names(tabulist) <- vars
  minaic <- vector("numeric")
  for(i in 1:length(tabulist)) {
    tabulist[[i]] <- 1000
  }

  iter <- 1
  stop.it <- 0

  bestscore <- 1
  if(test=="AIC") {
    bestscore <- 100000
  }
  
  Kbest <- NULL
  bestmodel <- NULL
  message("TABU search, type=", type, ", test=", test)
  while(stop.it < 30) {
    neigh <- NULL
    neigh <- neighbour(K, vars, tabulist, iter, tabusize)
    measurelist <- NULL
    models <- list()
    if(length(neigh)>1 && !is.null(neigh)) {
   
      for(i in 1:length(neigh)) {
        tdata <- subset(dataset,,neigh[[i]])
  

        if(type=="lm") {
          if(test=="Ftest") {
            model <- lm(dep ~., data=tdata)
            smmodel <- summary(model)
            measurelist[i] <- 1-pf(smmodel$fstatistic[1], smmodel$fstatistic[2], smmodel$fstatistic[3])
          }
          else {
            stop("Test not supported")
          }
        }
        else if(type=="cox") {
          if(test=="wald") {
            model <- coxph(Surv(dep[[1]], dep[[2]], dep[[3]]) ~ ., data=tdata)
            measurelist[i] <- summary(model)$waldtest[3]
          }
          else if(test=="AIC") {
            model <- coxph(Surv(dep[[1]], dep[[2]], dep[[3]]) ~ ., data=tdata)
            measurelist[i] <- extractAIC(model)[2]
          }
          
          else {
            stop("Test not supported")
          }
          
        }
        else {
          stop("Method not implemented")
        }
        models[[i]]<-model
      }

      if (min(measurelist) < bestscore) {
        bestscore <- min(measurelist)
        Kbest <- neigh[[which.min(measurelist)]]
        bestmodel <- models[[which.min(measurelist)]]
        for(i in 1:length(tabulist)) {
          tabulist[[i]] <- 1000
        }
      }
      
      Ktbest <- neigh[[which.min(measurelist)]]
      
    
      if(length(Ktbest)<length(K)) {
        tb <- K[which(!K%in%Ktbest)]
      }
      else {
        tb <- Ktbest[which(!Ktbest%in%K)]
      }

      tabulist[[tb]] <- iter

    }
    

    iter <- iter+1
    if(setequal(Ktbest,K)) {
      stop.it <- stop.it+1
    }
    else {
      stop.it <- 0
    }
    
    
    K <- Ktbest
    
  }

  message("Number of iterations : ", iter)
  message("Selected variables : ")
  print(Kbest)

  if(test=="Ftest") { message("F-score : ", summary(bestmodel)$fstatistic[1], " df = ", summary(bestmodel)$fstatistic[2], " p-value = ", bestscore) }
  else { 
	message("Best value of ", test, " is : ", bestscore)
	}
  return(list("K"=Kbest, "measure"=bestscore, "bestmodel"=bestmodel))

}


neighbour <- function(K, vars, tabulist, iter, tabusize) {
  #
  newK <- list()
  #newK[[1]] <- K
  n <- 1
  ## add variables
  for(i in 1:length(vars)) {
   # print(vars[i])
   # print(!(vars[i]%in%K) && !(vars[i]%in%tabulist))
    if(!(vars[i]%in%K) && ((tabulist[[vars[i]]] - iter) > tabusize)) {
     # print(paste("adding ", vars[i]))
      newK[[n]] <- c(K, vars[i])
      n<-n+1
    }
  }
  ## remove variables
  for(i in 1:length(K)) {
    if( ((tabulist[[vars[i]]] - iter) > tabusize) && (length(K) > 1)) {
      newK[[n]] <- K[-i]
      n<-n+1
    }
  }
  ## swap variables
  for(i in 1:length(K)) {
    for(j in 1:length(vars)) {
      if(!(vars[j]%in%K) && ( (tabulist[[vars[j]]] - iter) > tabusize ) && ( (tabulist[[K[i]]] - iter) > tabusize)  ) {
        #print(paste("deleting ", 
        newK[[n]] <- c(K[-i], vars[j])
        n<-n+1
      }

    }
  }
  return(newK)
}


tabusearch(c("bmi","map","hdl","ltg"), names(diab.train)[2:11], diab.train$Y, diab.train, 10)

lm.tabu <- lm(Y ~ bmi + map + hdl + ltg + age,diab.train)
summary(lm.tabu)

pred.tabu = predict(lm.tabu,diab.test)
cor(pred.tabu,diab.test$Y)**2


# from other TabuSearch program in R

lm.tabu2 <- lm(Y ~ sex + bmi + map + tc + ldl + ltg + glu,diab.train)
summary(lm.tabu2)

pred.tabu2 = predict(lm.tabu2,diab.test)
cor(pred.tabu2,diab.test$Y)**2

# compare to lasso

YY.train <- as.matrix(diab.train$Y)
XX.train <- as.matrix(diab.train[,2:11])

YY.test <- as.matrix(diab.test$Y)
XX.test <- as.matrix(diab.test[,2:11])

lasso.out9 <- cv.glmnet(XX.train,YY.train,family="gaussian",alpha=1)
(lmin <- lasso.out9$lambda.1se)

lasso.out99 <- glmnet(XX.train,YY.train,family="gaussian",alpha=1,lambda=lminSE)
coef(lasso.out99)

lm.lasso <- lm(Y ~ sex + bmi + map + hdl + ltg,diab.train)
pred.lasso2 = predict(lm.lasso,diab.test)
cor(pred.lasso2,diab.test$Y)**2


# and again compare to original regression

lm.full <- lm(Y ~ .,diab.train)
summary(lm.full)

pred.full = predict(lm.full,diab.test)
cor(pred.full,diab.test$Y)**2
