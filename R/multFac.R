#' Tests multiple factors
#'
#'
#' @param facList A vector containing the number of factors to test. Ex: ll = c(1,2,3)
#' @param parallel Whether to use the snowfall package for parallelization. Note that
#'        this is different than in autoSEM. Parallelization with multFac runs the
#'        different factor models separately, not in the actual search algorithm.
#' @param ncore Number of cores to use.
#' @param method which optimization algorithm to use. Currently, it is only
#'        recommended to use "GA" for the genetic algorithm from the GA
#'        package, "aco", an implementation of the ant colony
#'        algorithm by Ross Jacobucci, and "tabu", an implementation of
#'        the Tabu search procedure by Ross Jacobucci. The latter two
#'        algorithms are based on the book chapter by Marcoulides &
#'        Leite, 2013.
#' @param missing Argument to be passed to cfa() as to what to do with missing
#'        values. Note: missing="fiml" can't be paired with CV=TRUE
#' @param data a required dataset to search with.
#' @param varList list containing the names of the
#'        variables to use from the dataset.
#' @param criterion The fit index to use as a criterion for
#'        choosing the best model. Current options are "NCP",
#'        "RMSEA", and "BIC".
#' @param minInd The minimum number of indicators per factor.
#' @param niter The maximum number of iterations to use. "default" changes the number
#'        of iterations based on the algorithm used.
#' @param CV Whether to use cross-validation for choosing the best model. The
#'        default is to use fit indices without CV.
#' @param min.improve Number of iterations to wait for improvement
#'        before breaking.
#' @param seed random seed number.
#' @param std.lv Defaults to true. So lavaan uses all variables for each factor
#' @param ... Additional arguments to pass to cfa(). An example is
#'        is setting orth=FALSE,std.lv=TRUE.
#' @keywords multFac
#' @export
#' @examples
#' \dontrun{
#' library(autoSEM)
#' myData =  HolzingerSwineford1939[,7:15]
#'
#' f1.vars <- c("x1","x2","x3","x4","x5","x6","x7","x8","x9")
#' rrr = list(f1.vars)
#' facs <- 1:4
#'
#' out = multFac(facList=facs,parallel="yes",ncore=4,method="GA",
#'             data=myData,orth=FALSE,CV=FALSE,std.lv=TRUE,
#'             varList=rrr,criterion="RMSEA",niter="default")
#'}



multFac <- function(facList,
                    parallel="no",
                    ncore=1,
                    method="GA",
                    missing="listwise",
                    data=NULL,
                    varList=NULL,
                    criterion="BIC",
                    minInd=3,
                    niter="default",
                    CV=FALSE,
                    min.improve=niter,
                    seed=NULL,
                    std.lv=TRUE,
                    ...){

  if(niter == "default"){
    if(method=="tabu"){
      niter=30
    }else if(method =="aco"){
      niter=100
    }else if(method=="GA"){
      niter=20
    }
  }


  if(length(facList) < ncore){
    ncore = length(facList)
  }

  if(parallel=="yes"){

    if(Sys.info()[1] == "Windows"){

      snowfall::sfInit(parallel=TRUE, cpus=ncore)
      snowfall::sfExport("data","facList","criterion","CV",
                         "minInd","niter","varList","method")
      snowfall::sfLibrary(autoSEM); snowfall::sfLibrary(lavaan);
      snowfall::sfLibrary( "GA", character.only=TRUE )

      ret.auto <- function(facs){

        ret = autoSEM(method=method,data=data,nfac=facs,CV=CV,,std.lv=std.lv,missing=missing,
                      ...,
                      varList=varList,criterion=criterion,minInd=minInd,
                      niter=niter,min.improve=min.improve)
        ret
      }


      out = sfLapply(facList,ret.auto)
      snowfall::sfStop()

    }else if(Sys.info()[1]=="Darwin"){

    snowfall::sfStop()
    snowfall::sfInit(parallel=TRUE, cpus=ncore)
    snowfall::sfExport("data","facList","criterion","CV",
                       "minInd","niter","varList","method")
    snowfall::sfLibrary(autoSEM); snowfall::sfLibrary(lavaan);snowfall::sfLibrary( "GA", character.only=TRUE )

    ret.auto <- function(facs){

      ret = autoSEM(method=method,data=data,nfac=facs,CV=CV,,std.lv=std.lv,missing=missing,
                    ...,
                    varList=varList,criterion=criterion,minInd=minInd,
                    niter=niter,min.improve=min.improve)
      ret
    }


    out = sfLapply(facList,ret.auto)
    snowfall::sfStop()
    }
  }else if(parallel=="no"){
    out = list()

    for(y in 1:length(facList)){
      out[[y]] = autoSEM(method=method,data=data,nfac=facList[y],CV=CV,,std.lv=std.lv,missing=missing,
                         ...,
                         varList=varList,criterion=criterion,minInd=minInd,
                         niter=niter,min.improve=min.improve)
    }
  }

out
}


