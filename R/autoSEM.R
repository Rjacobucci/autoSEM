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
                    replaceSamp=FALSE){

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

    uuu = list()
    for(i in 1:length(lll)){
      uuu[[i]] = length(lll[[i]])
    }

    jjj =  list()
    for(i in 1:nfac){
      jjj[[i]] =  string[1:uuu[[i]]]
      string = string[-(1:uuu[[i]])]
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

    ooo = list()
  #  for(i in 1:nfac){
  #    facc = paste("f",i,sep="")
   #   ooo[[i]] = paste(paste(facc," =~ "), paste(lll[[i]][jjj[[i]]==1], collapse= "+"))
   # }

    for(i in 1:nfac){
          facc = paste("f",i,sep="")
         ooo[[i]] = paste(paste(facc," =~ "), paste(lll[[i]], collapse= "+"))
    }

    ppp <- list()
    for(i in 1:nfac){
      ppp[[i]] = gsub("x","NA*x",ooo[[i]])
    }


    for(i in 1:nfac){
      ppp[[i]] = gsub("NA",string,ooo[[i]])
    }

    for(i in 1:9){
      gsub("NA",ppp[[1]]
    }


    fmld <- ""
    for(jj in 1:nfac){
      fmld <- paste(fmld,ppp[[jj]],sep="\n")
    }

    p = length(unique(unlist(varList)))

    outt = lavaan::cfa(fmld,data_train,orthogonal=orth,std.lv=stdlv)

    if(inspect(outt,"converged")==F | any(eigen(inspect(outt,"cov.lv"))$values < 0)){

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

      if(criterion=="BIC"){
        return_val = 1/bic
      }else if(criterion=="RMSEA"){
        return_val= -RMSEA.test
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
      out = ga("binary", fitness = fitness, nBits = p_length,monitor=T,maxiter=niter,run=10)
    }else if(parallel=="yes"){
      out = ga("binary", fitness = fitness, nBits = p_length,monitor=T,maxiter=niter,parallel=TRUE)
    }

  }else if(method=="tabuSearch"){
    out = tabuSearch(size = p_length, iters = niter,objFunc = fitness)
  }

  out
}
