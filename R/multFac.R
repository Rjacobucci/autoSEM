#'
#'
#' Tests multiple factors
#' @param method
#' @keywords multFac
#' @export
#' @examples
#' multFac()
#'
#'


multFac <- function(facList,parallel="no",ncore=1,
                    method="tabuSearch",
                    data=NULL,
                    varList=NULL,
                    criterion="BIC",
                    minInd=3,
                    stdlv=TRUE,
                    orth=TRUE,
                    niter=30,
                    CV=FALSE){

  if(parallel=="yes"){
    #library(snowfall)
    snowfall::sfInit(T,ncore)
    snowfall::sfExport("data","facList","criterion","CV","orth",
                       "stdlv","minInd","niter","varList","method")
    snowfall::sfLibrary(autoSEM); snowfall::sfLibrary(lavaan);
    snowfall::sfLibrary(GA);snowfall::sfLibrary(tabuSearch)

    ret.auto <- function(facs){

      ret = autoSEM(method=method,data=data,nfac=facs,orth=orth,CV=CV,
                    varList=varList,criterion=criterion,minInd=minInd,niter=niter)
      ret
    }


    out = sfLapply(facList,ret.auto)
    snowfall::sfStop()
  }else if(parallel=="no"){
    out = list()

    for(y in 1:length(facList)){
      out[[y]] = autoSEM(method=method,data=data,nfac=facList[y],orth=orth,CV=CV,
                         varList=varList,criterion=criterion,minInd=minInd,niter=niter)
    }
  }

out
}


