#' This function houses a number of different heuristic optimization algorithms
#'     for specification search.
#'
#'
#' @param method which optimization algorithm to use. Currently, it is only
#'        recommended to use "GA" for the genetic algorithm from the GA
#'        package, "aco", an implementation of the ant colony
#'        algorithm by Ross Jacobucci, and "tabu", an implementation of
#'        the Tabu search procedure by Ross Jacobucci. The latter two
#'        algorithms are based on the book chapter by Marcoulides &
#'        Leite, 2013.
#'
#' @param data a required dataset to search with.
#' @param nfac the number of factors to test.
#' @param varList list containing the names of the
#'        variables to use from the dataset.
#' @param criterion The fit index to use as a criterion for
#'        choosing the best model. Current options are "NCP",
#'        "RMSEA","AIC", "BIC", and "BIC2", which is the sample
#'        size adjusted BIC.
#' @param minInd The minimum number of indicators per factor.
#' @param niter The maximum number of iterations to use. "default" changes the number
#'        of iterations based on the algorithm used.
#' @param parallel Whether to use the snowfall package for parallelization.
#'        Note that this is only applicable for the GA package at this time.
#' @param missing Argument to be passed to cfa() as to what to do with missing
#'        values. Note: missing="fiml" can't be paired with CV=TRUE
#' @param CV Whether to use cross-validation for choosing the best model. The
#'        default is to use fit indices without CV.It is currently recommended to either
#'        use FALSE or "boot". Note that "boot" will take significantly longer.
#' @param R If using bootstrap, how many samples to take? Default is 100
#' @param min.improve Number of iterations to wait for improvement
#'        before breaking.
#' @param seed random seed number.
#' @param std.lv Defaults to true. So lavaan uses all variables for each factor
#' @param ... Additional arguments to pass to cfa(). An example is
#'        is setting orth=FALSE,std.lv=TRUE.
#' @return fit the fit index
#' @return solution the solution with the best fit
#' @return out returned object from optimization algorithm
#' @keywords autoSEM
#' @export
#' @import lavaan
#' @import snowfall
#' @import GA
#' @importFrom stats cov dchisq rbinom runif
#' @examples
#' library(autoSEM)
#' myData =  HolzingerSwineford1939[,7:15]
#'
#' f1.vars <- c("x1","x2","x3","x4","x5","x6","x7","x8","x9")
#'
#' out = autoSEM(method="GA",data=myData,nfac=1,
#'              varList=list(f1.vars),CV=FALSE,
#'              criterion="RMSEA",minInd=3,niter=1)
#' summary(out)
#'
autoSEM <- function(method="GA",
                    data=NULL,
                    nfac=NULL,
                    varList=NULL,
                    criterion="BIC",
                    minInd=3,
                    niter="default",
                    parallel="no",
                    missing="listwise",
                    CV="boot",
                    R=100,
                    min.improve=niter,
                    seed=NULL,
                    std.lv=TRUE,
                    ...){
  ret <- list()
  options(warn=2)


  if(missing == "fiml" & CV == TRUE){
    stop("Can't pair fiml with cross-validation at this time")
  }


  if(niter == "default"){
    if(method=="tabu"){
      niter=length(varList)*20
    }else if(method =="aco"){
      niter=100
    }else if(method=="GA"){
      niter=20
    }
  }

  if(is.null(seed)==TRUE) seed = round(runif(1,0,1)*10000000,0)

 if(method != "GA" & method != "tabu" & method != "aco"){
    stop("Only GA, tabu, and aco are currently implemented")
  }

  if(CV==T){
    ids = sample(nrow(data),nrow(data)/2)
    data_train = data[ids,]
    data_test = data[-ids,]
  }else if(CV==FALSE){
    data_train = data
  }else if(CV=="boot"){
    data_train = data
  }else{
    stop("Wrong Assignment to CV")
  }
#  }else if(replaceSamp==TRUE){
#    ids = sample(nrow(data),nrow(data),replace=TRUE)
#    ids2 = sample(nrow(data),nrow(data),replace=TRUE)
#    data_train = data[ids,]
#    data_test = data[ids2,]
#  }


  fitness <- function(string) {

    lll = varList


    jjj = list()
    for(i in 1:nfac){
      jjj[[i]] = lll
    }

    jjj =  list()
    for(i in 1:nfac){
      jjj[[i]] =  string[1:length(lll[[1]])]
      string = string[-(1:length(lll[[1]]))]
    }

        for(i in 1:length(jjj)){
          if(sum(jjj[[i]]) < minInd){
            if(method=="GA"){
              -44
            }else if(method=="tabu" | method=="aco"){
              1e10
            }
          }
        }





  ooo =  list()

    for(i in 1:nfac){
          facc = paste("f",i,sep="")
          uu = gsub("1","start(1)*",jjj[[i]])
          uu = gsub("0","0*",uu)
          uu2 = paste(uu,lll[[1]],sep="")
          ooo[[i]] =  paste0(uu2,collapse="+")
    }

    for(i in 1:nfac){
      facc = paste("f",i,sep="")
      ooo[[i]] = paste(paste(facc," =~ "), ooo[[i]])
    }


   # ppp <- list()
   # for(i in 1:nfac){
    #  ppp[[i]] = gsub("x","NA*x",ooo[[i]])
   # }


   # for(uu in 1:length(start)){
   #   ppp[[1]] = gsub("NA",start[uu],ppp[[1]])
    #  ppp[[1]]
    #}



    fmld <- ""
    for(jj in 1:nfac){
      fmld <- paste(fmld,ooo[[jj]],sep="\n")
    }

    p = length(unique(unlist(varList)))


    d = function(chisq,df,N) max(0,(chisq -df)/(N-1))
    rmsea = function(ncp,df) sqrt(ncp/df)


    outt = try(lavaan::cfa(fmld,data_train,std.lv=std.lv,missing=missing,...),silent=TRUE)
    #summary(outt)
    if(inherits(outt, "try-error")) {
      if(method=="GA"){
        #return(-99999999)
        -1e10
      }else if(method=="tabu" | method=="aco"){
        1e10
      }
      }else{
      if(CV==F){
        fits.lav = lavaan::fitMeasures(outt)
        bic = fits.lav["bic"]
        bic2 = fits.lav["bic2"]
        aic = fits.lav["aic"]
        RMSEA = fits.lav["rmsea"]
        NCP = d(fits.lav["chisq"],fits.lav["df"],fits.lav["ntotal"])
      }else if(CV==T){

        if(criterion=="BIC" | criterion=="BIC2" | criterion=="AIC"){
          stop("Only use CV=F with BIC, BIC2, or AIC")
        }

        df=outt@Fit@test[[1]]$df
        cov.order = outt@Data@ov.names
        cov.test = cov(data_test[,unlist(cov.order)])
        impcov = fitted(outt)$cov
        N=nrow(data_test)
        fit.test = 0.5*(log(det(impcov)) + trace(cov.test %*% solve(impcov)) - log(det(cov.test))  - p)
        chisq.test = N*fit.test
        NCP = d(chisq.test,df,N)
        RMSEA = rmsea(NCP,df)
      }else if(CV == "boot"){

       # if(criterion=="BIC" | criterion=="BIC2" | criterion=="AIC"){
        #  stop("Only use CV=F with BIC, BIC2, or AIC")
       # }

       # RMSEA.rep <- rep(NA,R)
        #NCP.rep <- rep(NA,R)
        #impcov = fitted(outt)$cov
        N=nrow(data_train)
        df=outt@Fit@test[[1]]$df
        #for(i in 1:100){
        #  ids <- sample(nrow(data_train),nrow(data_train),replace=TRUE)
        #  new.dat <- data_train[ids,]
        #  cov.boot <- cov(new.dat[,outt@pta$vnames$ov.nox[[1]]])
        #  fit.boot = 0.5*(log(det(impcov)) + trace(cov.boot %*% solve(impcov)) - log(det(cov.boot))  - p)
        #  chisq.boot = N*fit.boot
        #  NCP.rep[i] = d(chisq.boot,df,N)
        #  RMSEA.rep[i] = rmsea(NCP.rep[i],df)
        #}
        chisq.boot = rep(NA,R)
        NCP.rep = rep(NA,R)
        #RMSEA.rep = rep(NA,R)

        fit.boot = try(lavaan::bootstrapLavaan(outt,type="yuan",R=R,FUN=fitmeasures,
                                       fit.measures=c("chisq","rmsea","bic","bic2","aic")),silent=TRUE)
        #print(fit.boot)
        if(inherits(fit.boot, "try-error")) {
          NCP = 9987
          RMSEA = 9987
          bic <- 1e10
          bic2 = 1e10
          aic = 1e10
        }else{
          RMSEA <- mean(fit.boot[,2])
          for(i in 1:R){
            NCP.rep[i] = d(fit.boot[i,1],df,N)
          }
          NCP <- mean(NCP.rep)
          bic <- mean(fit.boot[,3])
          bic2 <- mean(fit.boot[,4])
          aic <- mean(fit.boot[,5])
        }


      }

    if(method== "GA"){
      if(criterion=="BIC"){
        return_val = 100 /bic
      }else if(criterion=="BIC2"){
        return_val = 100 /bic2
      }else if(criterion=="AIC"){
        return_val = 100 /aic
      }else if(criterion=="RMSEA"){
        return_val= 100/RMSEA
      }else if(criterion=="NCP"){
        return_val = 100/NCP
      }
    }else if(method=="tabu" | method=="aco"){
      if(criterion=="BIC"){
        return_val = bic
      }else if(criterion=="BIC2"){
        return_val = bic2
      }else if(criterion=="AIC"){
        return_val = aic
      }else if(criterion=="RMSEA"){
        return_val= RMSEA
      }else if(criterion=="NCP"){
        return_val=NCP
      }
    }

#10
      if(method=="GA"){
        return(return_val)
        #10
      }else if(method=="tabu" | method=="aco"){
        return(return_val)
      }

    }
   # 10
  }



  p_length = length(unlist(varList))

  if(method=="GA"){
   # if(parallel=="no"){
      out = GA::ga("binary", fitness = fitness, nBits = p_length*nfac,monitor=T,maxiter=niter,run=10)
      # }else if(parallel=="yes"){
      #   out = GA::ga("binary", fitness = fitness, nBits = p_length*nfac,monitor=T,maxiter=niter,parallel=TRUE)
      # }
  }else if(method=="tabu"){
    out = tabu(size=p_length*nfac,iters=niter,fitness=fitness,min.improve=min.improve,seed=seed)
  }else if(method=="aco"){
    out = aco(size=p_length*nfac,iters=niter,fitness=fitness,
                 criterion=criterion,min.improve=min.improve,seed=seed)
  }

  fit= NULL
  if(method=="tabu" | method=="aco"){
    if(criterion=="BIC" | criterion=="BIC2" | criterion=="AIC"){
      ret$fit = out$value
    }else if(criterion=="RMSEA"){
      ret$fit = out$value
    }else if(criterion=="NCP"){
      ret$fit = out$value
    }
  }

  if(method == "GA" | method=="aco" | method=="tabu"){
    ret$out = out
  }

  if(method=="GA"){
    ret$solution = out@solution
    if(criterion=="BIC"){
      ret$fit <- 100/out@fitnessValue
    }else{
      ret$fit <- 100/out@fitnessValue
    }

  }else if(method=="tabu" | method=="aco"){
    ret$solution = out$solution
  }


  ret$call <- match.call()
  class(ret) <- "autoSEM"
  return(ret)
}
