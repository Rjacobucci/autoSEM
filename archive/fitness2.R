#'
#'
#' fitness2 just for tabu_rj
#' @param string
#' @keywords fitness2
#' @export
#' @examples
#' fitness2()
#'
#'
#'

fitness2 <- function(string,varList,method,nfac,data_train,orth,stdlv,CV,criterion) {

  lll = varList

  if(method=="rgenoud" | method=="pso" | method=="NMOF" | method=="DEoptim"){
    string = round(string)
  }

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
  #          -44
  #        }else if(method=="tabuSearch"){
  #          0
  #        }else if(method=="rgenoud"){
  #          99999999
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
      #return(-99999999)
      -99999999
    }else if(method=="tabuSearch"){
      #return(0)
      0
    }else if(method == "rgenoud" | method=="pso" | method=="NMOF" |
             method=="DEoptim" | method=="tabu_rj"){
      9999999999
    }
  }else{
    if(CV==F){
      bic = fitMeasures(outt)["bic"]
      RMSEA = fitMeasures(outt)["rmsea"]
    }else if(CV==T){
      df=outt@Fit@test[[1]]$df
      cov.order = outt@Data@ov.names
      cov.test = cov(data_test[,unlist(cov.order)])
      impcov = fitted(outt)$cov
      N=nrow(data_test)
      fit.test = 0.5*(log(det(impcov)) + trace(cov.test %*% solve(impcov)) - log(det(cov.test))  - p)
      chisq.test = N*fit.test
      ncp.test = d(chisq.test,df,N)
      RMSEA = rmsea(ncp.test,df)
    }else if(CV == "boot"){
      RMSEA.rep <- rep(NA,100)
      impcov = fitted(outt)$cov
      N=nrow(data_train)
      df=outt@Fit@test[[1]]$df
      for(i in 1:100){
        ids <- sample(nrow(data_train),nrow(data_train),replace=TRUE)
        new.dat <- data_train[ids,]
        cov.boot <- cov(new.dat[,outt@pta$vnames$ov.nox[[1]]])
        fit.boot = 0.5*(log(det(impcov)) + trace(cov.boot %*% solve(impcov)) - log(det(cov.boot))  - p)
        chisq.boot = N*fit.boot
        ncp.boot = d(chisq.boot,df,N)
        RMSEA.rep[i] = rmsea(ncp.boot,df)
      }
      RMSEA <- mean(RMSEA.rep)
    }

    if(method == "tabuSearch" | method== "GA"){
      if(criterion=="BIC"){
        return_val = 100 /bic
      }else if(criterion=="RMSEA"){
        return_val= 1 - RMSEA
      }
    }else if(method=="rgenoud" | method == "pso" |
             method == "NMOF" | method=="DEoptim" |
             method=="tabu_rj"){
      if(criterion=="BIC"){
        return_val = bic
      }else if(criterion=="RMSEA"){
        return_val= RMSEA
      }
    }

    #10
    if(method=="GA"){
      return(return_val)
      #10
    }else if(method=="tabuSearch"){
      return(return_val)
      #10
    }else if(method=="rgenoud" | method=="pso" | method=="NMOF" | method=="DEoptim"){
      return(return_val)
      #10
    }else if(method=="tabu_rj"){
      return(return_val)
    }
  }
  # 10
}


