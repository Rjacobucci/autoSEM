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
autoSEM <- function(method="tabuSearch",nfac=NULL,varNames=NULL,criterion="BIC"){

  fitness <- function(string) {
    pt1 = string[1:6]
    pt2 = string[7:12]
    pt3 = string[13:18]

    if(sum(pt1) < 3 | sum(pt2) < 3 | sum(pt3) < 3){

      if(method=="GA"){
        return(-99999999)
      }else if(method="tabuSearch"){
        return(1)
      }
    }

    inc1 <- which(pt1 == 1)
    inc2 <- which(pt2 == 1)
    inc3 <- which(pt3 == 1)

    f1.vars <- c("x1","x2","x3","x4","x5","x6")
    f2.vars <- c("x4","x5","x6","x7","x8","x9")
    f3.vars <- c("x1","x2","x3","x7","x8","x9")

    fmld <- c(paste("factor1 =~ ", paste(f1.vars[inc1], collapse= "+")),
              paste("factor2 =~ ", paste(f2.vars[inc2], collapse= "+")),
              paste("factor3 =~ ", paste(f3.vars[inc3], collapse= "+")))


    #vars = unique(c(f1.vars[inc1],f2.vars[inc2],f3.vars[inc3]))
    #cov.test = cov(myData.test[vars])


    out = lavaan::cfa(fmld,myData,orthogonal=T,std.lv=T)

    if(inspect(out,"converged")==F | any(eigen(inspect(out,"cov.lv"))$values < 0)){

      if(method=="GA"){
        return(-99999999)
      }else if(method="tabuSearch"){
        return(1)
      }

    }else{
      bic = fitMeasures(out)["bic"]
      df=out@Fit@test[[1]]$df
      cov.order = out@Data@ov.names
      cov.test = cov(myData.test[,unlist(cov.order)])
      impcov = fitted(out)$cov
      fit.test = 0.5*(log(det(impcov)) + trace(cov.test %*% solve(impcov)) - log(det(cov.test))  - 9)
      chisq.test = N*fit.test
      ncp.test = d(chisq.test,df,N)
      RMSEA.test = rmsea(ncp.test,df)

      if(criterion=="BIC"){
        return_val = -bic
      }else if(criterion=="RMSEA"){
        return_val= -RMSEA.test
      }


      if(method=="GA"){
        return(return_val + 1000000000000)
      }else if(method="tabuSearch"){
        return(return_val + 1000000000000)
      }

    }
  }


  if(method=="GA"){
    out = ga("binary", fitness = fitness, nBits = p_length,monitor=T)
  }else if(method="tabuSearch"){
    out = tabuSearch(size = p_length, iters = 50,objFunc = fitness)
  }

  out
}
