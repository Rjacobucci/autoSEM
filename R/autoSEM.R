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
                    CV=FALSE,
                    aco.iters=100){
  ret <- list()
  options(warn=2)

  #if(method != "tabuSearch" & method != "GA"){
  #  stop("Only tabuSearch and GA are currently working well.")
  #}

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

        for(i in 1:length(jjj)){
          if(sum(jjj[[i]]) < minInd){
            if(method=="GA"){
              -44
            }else if(method=="tabuSearch"){
              0
            }else if(method=="rgenoud" | method=="tabu_rj" | method=="aco_rj"){
              99999999
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


    outt = try(lavaan::cfa(fmld,data_train,orthogonal=orth,std.lv=stdlv),silent=TRUE)

    if(inherits(outt, "try-error")) {
      if(method=="GA"){
        #return(-99999999)
        -99999999
      }else if(method=="tabuSearch"){
        #return(0)
        0
      }else if(method == "rgenoud" | method=="pso" | method=="NMOF" |
               method=="DEoptim" | method=="tabu_rj" | method=="aco_rj"){
        9999999999
      }
      }else{
      if(CV==F){
        fits.lav = fitMeasures(outt)
        bic = fits.lav["bic"]
        RMSEA = fits.lav["rmsea"]
        NCP = d(fits.lav["chisq"],fits.lav["df"],fits.lav["ntotal"])
      }else if(CV==T){
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
        RMSEA.rep <- rep(NA,100)
        NCP.rep <- rep(NA,100)
        impcov = fitted(outt)$cov
        N=nrow(data_train)
        df=outt@Fit@test[[1]]$df
        for(i in 1:100){
          ids <- sample(nrow(data_train),nrow(data_train),replace=TRUE)
          new.dat <- data_train[ids,]
          cov.boot <- cov(new.dat[,outt@pta$vnames$ov.nox[[1]]])
          fit.boot = 0.5*(log(det(impcov)) + trace(cov.boot %*% solve(impcov)) - log(det(cov.boot))  - p)
          chisq.boot = N*fit.boot
          NCP.rep[i] = d(chisq.boot,df,N)
          RMSEA.rep[i] = rmsea(NCP.rep[i],df)
        }
        RMSEA <- mean(RMSEA.rep)
        NCP <- mean(NCP.rep)

      }

    if(method == "tabuSearch" | method== "GA"){
      if(criterion=="BIC"){
        return_val = 100 /bic
      }else if(criterion=="RMSEA"){
        return_val= 1 - RMSEA
      }else if(criterion=="NCP"){
        return_val = 10/NCP
      }
    }else if(method=="rgenoud" | method == "pso" |
             method == "NMOF" | method=="DEoptim" |
             method=="tabu_rj" | method=="aco_rj"){
      if(criterion=="BIC"){
        return_val = bic
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
      }else if(method=="tabuSearch"){
        return(return_val)
        #10
      }else if(method=="rgenoud" | method=="pso" | method=="NMOF" | method=="DEoptim"){
        return(return_val)
        #10
      }else if(method=="tabu_rj" | method=="aco_rj"){
        return(return_val)
      }
    }
   # 10
  }



  p_length = length(unlist(varList))

  if(method=="GA"){
    if(parallel=="no"){
      out = ga("binary", fitness = fitness, nBits = p_length*nfac,monitor=T,maxiter=niter,run=10)
    }else if(parallel=="yes"){
      out = ga("binary", fitness = fitness, nBits = p_length*nfac,monitor=T,maxiter=niter,parallel=TRUE)
    }
  }else if(method=="tabuSearch"){
    out = tabuSearch(size = p_length*nfac, iters = niter,objFunc = fitness,listSize=5)
  }else if(method=="tabu_rj"){
    out = tabu_rj(size=p_length*nfac,iters=niter,fitness=fitness)
  }else if(method=="aco_rj"){
    out = aco_rj(size=p_length*nfac,iters=aco.iters,fitness=fitness,criterion=criterion)
  }else if(method=="rgenoud"){
    dom = cbind(rep(0,p_length*nfac),rep(1,p_length*nfac))
    out = rgenoud::genoud(fitness,nvars=p_length*nfac,Domains=dom,boundary=2,print.level=0)
  }else if(method=="pso"){

    out = hydroPSO::hydroPSO(rep(0.5,p_length*nfac),fitness,
                             lower=rep(0,p_length*nfac),upper=rep(1,p_length*nfac))

    #out = pso::psoptim(rep(0.5,p_length*nfac),fitness,
    #                         lower=0,upper=1)
  }else if(method=="NMOF"){
    fitness(rep(0.5,p_length*nfac))
    algo=list(nB=p_length*nfac)
    out = NMOF::GAopt(fitness,algo=algo)
  }else if(method=="DEoptim"){
    out = RcppDE::DEoptim(fitness,lower=rep(0,p_length*nfac),upper=rep(1,p_length*nfac),
                          control=DEoptim.control(NP=9000))
  }

  if(method=="GA"){
    if(criterion=="BIC"){
      ret$fit = 100/(summary(out)$fitness)
    }else if(criterion=="RMSEA"){
      ret$fit = 1- (summary(out)$fitness)
    }else if(criterion=="NCP"){
      ret$fit = 1/NCP
    }
  }else if(method=="rgenoud" | method=="pso" | method=="tabu_rj" | method=="aco_rj"){
    if(criterion=="BIC"){
      ret$fit = out$value
    }else if(criterion=="RMSEA"){
      ret$fit = out$value
    }else if(criterion=="NCP"){
      ret$fit = out$value
    }
  }else if(method=="tabuSearch"){
    if(criterion=="BIC"){
      ret$fit = 100/max(out$eUtilityKeep)
    }else if(criterion=="RMSEA"){
      ret$fit = 1-max(out$eUtilityKeep)
    }
  }else if(method=="NMOF"){
    if(criterion=="BIC"){
      ret$fit = out$OFvalue
    }else if(criterion=="RMSEA"){
      ret$fit = out$OFvalue
    }
  }

  if(method == "GA" | method == "rgenoud" | method=="pso" |
     method=="NMOF" | method=="aco_rj"){
    ret$out = out
  }

  if(method=="tabuSearch"){
    ret$solution = out$configKeep[out$eUtilityKeep == max(out$eUtilityKeep),]
  }else if(method=="GA"){
    ret$solution = summary(out)$solution
  }else if(method=="rgenoud" | method=="pso"){
    ret$solution = round(out$par)
  }else if(method=="NMOF"){
    ret$solution = out$xbest
  }else if(method=="tabu_rj" | method=="aco_rj"){
    ret$solution = out$solution
  }

  ret
}
