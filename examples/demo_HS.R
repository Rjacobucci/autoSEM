library(lavaan)
myData =  HolzingerSwineford1939[,7:15]



f1.vars <- c("x1","x2","x3","x4","x5","x6","x7","x8","x9")
rrr = list(f1.vars)
facs <- 1:4

ret.auto <- function(facs){
  res <- list()
  ret = autoSEM(method="GA",data=myData,nfac=facs,orth=FALSE,replaceSamp=TRUE,
                  varList=rrr,criterion="RMSEA",minInd=3,niter=60)
  #res$fit = 1/(summary(ret)$fitness-1)
  res$fit = -(summary(ret)$fitness-1)
  res$solution = summary(ret)$solution
  res
}

library(snowfall)
sfStop()
sfInit(T,4)
snowfall::sfExport("myData","rrr","facs")
sfLibrary(autoSEM); sfLibrary(lavaan); sfLibrary(GA)
(out = sfLapply(facs,ret.auto))
