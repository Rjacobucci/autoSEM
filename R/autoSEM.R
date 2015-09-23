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
autoSEM <- function(method="tabuSearch",nfac=NULL,varList=NULL,
                    criterion="BIC",minInd=3,stdlv=TRUE,orth=TRUE,
                    niter=30){

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

    for(i in 1:length(jjj)){
      if(sum(jjj[[i]]) < minInd){
        if(method=="GA"){
          return(-99999999)
        }else if(method=="tabuSearch"){
          return(1)
        }
      }
    }

    ooo = list()
    for(i in 1:nfac){
      facc = paste("f",i,sep="")
      ooo[[i]] = paste(paste(facc," =~ "), paste(lll[[i]][jjj[[i]]==1], collapse= "+"))
    }

    fmld <- ""
    for(jj in 1:nfac){
      fmld <- paste(fmld,ooo[[jj]],sep="\n")
    }

    p = length(unique(unlist(varList)))

    outt = lavaan::cfa(fmld,myData,orthogonal=orth,std.lv=stdlv)

    if(inspect(outt,"converged")==F | any(eigen(inspect(outt,"cov.lv"))$values < 0)){

      if(method=="GA"){
        return(-99999999)
      }else if(method=="tabuSearch"){
        return(1)
      }

    }else{
      bic = fitMeasures(outt)["bic"]
     # df=outt@Fit@test[[1]]$df
      #cov.order = outt@Data@ov.names
      #cov.test = cov(myData.test[,unlist(cov.order)])
     # impcov = fitted(out)$cov
      #fit.test = 0.5*(log(det(impcov)) + trace(cov.test %*% solve(impcov)) - log(det(cov.test))  - p)
      #chisq.test = N*fit.test
     # ncp.test = d(chisq.test,df,N)
     # RMSEA.test = rmsea(ncp.test,df)

      if(criterion=="BIC"){
        return_val = -bic
      }else if(criterion=="RMSEA"){
        return_val= -RMSEA.test
      }


      if(method=="GA"){
        return(return_val + 1000000000000)
      }else if(method=="tabuSearch"){
        return(return_val + 1000000000000)
      }

    }
  }
  p_length = length(unlist(varList))

  if(method=="GA"){
    out = ga("binary", fitness = fitness, nBits = p_length,monitor=T,maxiter=niter)
  }else if(method=="tabuSearch"){
    out = tabuSearch(size = p_length, iters = niter,objFunc = fitness)
  }

  out
}
