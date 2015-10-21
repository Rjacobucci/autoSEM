#'
#'
#' This function houses a number of different heuristic optimization algorithms for specification search.
#' @param method which optimization algorithm to use
#' @keywords autoSEM
#' @export
#' @examples
#' autoSEM()
#'
#'
autoSEM <- function(method="tabuSearch",
                    data=NULL,
                    nfac=NULL,
                    varList=NULL,
                    criterion="BIC",
                    minInd=3,
                    stdlv=TRUE,
                    orth=TRUE,
                    niter=30,
                    parallel="no",
                    replaceSamp=FALSE,
                    boot=FALSE){

  options(warn=2)

  if(replaceSamp==FALSE){
    ids = sample(nrow(data),nrow(data)/2)
    data_train = data[ids,]
    data_test = data[-ids,]
  }else if(replaceSamp==TRUE){
    ids = sample(nrow(data),nrow(data),replace=TRUE)
    ids2 = sample(nrow(data),nrow(data),replace=TRUE)
    data_train = data[ids,]
    data_test = data[ids2,]
  }


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

#    for(i in 1:length(jjj)){
#      if(sum(jjj[[i]]) < minInd){
#        if(method=="GA"){
#          return(-99999999)
#        }else if(method=="tabuSearch"){
#          return(0)
#        }
#      }
#    }





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


    outt = try(lavaan::cfa(fmld,data_train,orthogonal=orth,std.lv=stdlv),silent=TRUE)

    if(inherits(outt, "try-error")) {
      if(method=="GA"){
        return(-99999999)
      }else if(method=="tabuSearch"){
        return(0)
      }
      }else{
      bic = fitMeasures(outt)["bic"]
      df=outt@Fit@test[[1]]$df
      cov.order = outt@Data@ov.names
      cov.test = cov(data_test[,unlist(cov.order)])
      impcov = fitted(outt)$cov
      N=nrow(data_test)
      fit.test = 0.5*(log(det(impcov)) + trace(cov.test %*% solve(impcov)) - log(det(cov.test))  - p)
      chisq.test = N*fit.test
      ncp.test = d(chisq.test,df,N)
      RMSEA.test = rmsea(ncp.test,df)

      if(boot == TRUE){
        RMSEA.rep <- rep(NA,100)
        for(i in 1:100){
          ids <- sample(nrow(myData),nrow(myData),replace=TRUE)
          new.dat <- myData[ids,]
          cov.boot <- cov(new.dat)
          fit.boot = 0.5*(log(det(impcov)) + trace(cov.boot %*% solve(impcov)) - log(det(cov.boot))  - p)
          chisq.boot = N*fit.boot
          ncp.boot = d(chisq.boot,df,N)
          RMSEA.rep[i] = rmsea(ncp.boot,df)
        }
      }
      RMSEA.boot <- mean(RMSEA.rep)

      if(criterion=="BIC"){
        return_val = 1/bic
      }else if(criterion=="RMSEA"){
        return_val= -RMSEA.test
      }else if(criterion=="RMSEA.boot"){
        return_val= -RMSEA.boot
      }


      if(method=="GA"){
        return(return_val + 1)
      }else if(method=="tabuSearch"){
        return(return_val + 1)
      }

    }
  }
  p_length = length(unlist(varList))

  if(method=="GA"){
    if(parallel=="no"){
      out = ga("binary", fitness = fitness, nBits = p_length*nfac,monitor=T,maxiter=niter,run=10)
    }else if(parallel=="yes"){
      out = ga("binary", fitness = fitness, nBits = p_length*nfac,monitor=T,maxiter=niter,parallel=TRUE)
    }

  }else if(method=="tabuSearch"){
    out = tabuSearch(size = p_length*nfac, iters = niter,objFunc = fitness)
  }

  out
}
