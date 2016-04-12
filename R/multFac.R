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
#'        package, "aco_rj", an implementation of the ant colony
#'        algorithm by Ross Jacobucci, and "tabu_rj", an implementation of
#'        the Tabu search procedure by Ross Jacobucci. The latter two
#'        algorithms are based on the book chapter by Marcoulides &
#'        Leite, 2013. The other methods: "pso", "NMOF", "DEoptim",
#'        "tabuSearch", and "rgenoud" are all based on real-value
#'        optimization, not binary strings. This substantially increases
#'        the computation time and these methods are not currently
#'        recommended for use.
#' @param data a required dataset to search with.
#' @param varList list containing the names of the
#'        variables to use from the dataset.
#' @param criterion The fit index to use as a criterion for
#'        choosing the best model. Current options are "NCP",
#'        "RMSEA", and "BIC".
#' @param minInd The minimum number of indicators per factor.
#' @param lav.control Additional arguments to pass to cfa(). An example is
#'        is setting lav.control = list(orth=FALSE,std.lv=TRUE).
#' @param niter The maximum number of iterations to all.
#' @param CV Whether to use cross-validation for choosing the best model. The
#'        default is to use fit indices without CV.
#' @param min.improve Number of iterations to wait for improvement
#'        before breaking.
#' @param seed random seed number.
#' @keywords multFac
#' @export
#' @examples
#' \dontrun{
#' multFac()
#'}



multFac <- function(facList,
                    parallel="no",
                    ncore=1,
                    method="tabuSearch",
                    data=NULL,
                    varList=NULL,
                    criterion="BIC",
                    minInd=3,
                    lav.control = list(orth=TRUE,std.lv=TRUE),
                    niter=30,
                    CV=FALSE,
                    min.improve=niter,
                    seed=NULL){

  if(length(facList) < ncore){
    ncore = length(facList)
  }

  if(parallel=="yes"){

    if(Sys.info()[1] == "Windows"){

      snowfall::sfInit(parallel=TRUE, cpus=ncore)
      snowfall::sfExport("data","facList","criterion","CV","orth",
                         "stdlv","minInd","niter","varList","method")
      snowfall::sfLibrary(autoSEM); snowfall::sfLibrary(lavaan);
      snowfall::sfLibrary(GA);#snowfall::sfLibrary(tabuSearch);
      #snowfall::sfLibrary(rgenoud)

      ret.auto <- function(facs){

        ret = autoSEM(method=method,data=data,nfac=facs,CV=CV,
                      lav.control=lav.control,
                      varList=varList,criterion=criterion,minInd=minInd,
                      niter=niter,min.improve=min.improve)
        ret
      }


      out = sfLapply(facList,ret.auto)
      snowfall::sfStop()

    }else if(Sys.info()[1]=="Darwin"){

    snowfall::sfStop()
    snowfall::sfInit(parallel=TRUE, cpus=ncore)
    snowfall::sfExport("data","facList","criterion","CV","orth",
                       "stdlv","minInd","niter","varList","method")
    snowfall::sfLibrary(autoSEM); snowfall::sfLibrary(lavaan);#snowfall::sfLibrary(hydroPSO)
    snowfall::sfLibrary(GA);#snowfall::sfLibrary(tabuSearch);snowfall::sfLibrary(rgenoud)

    ret.auto <- function(facs){

      ret = autoSEM(method=method,data=data,nfac=facs,CV=CV,
                    lav.control=lav.control,
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
      out[[y]] = autoSEM(method=method,data=data,nfac=facList[y],CV=CV,
                         lav.control=lav.control,
                         varList=varList,criterion=criterion,minInd=minInd,
                         niter=niter,min.improve=min.improve)
    }
  }

out
}


