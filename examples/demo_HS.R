library(lavaan)
myData =  HolzingerSwineford1939[,7:15]



f1.vars <- c("x1","x2","x3","x4","x5","x6","x7","x8","x9")
rrr = list(f1.vars)
facs <- 2:3


uu = multFac(facList=facs,parallel="yes",ncore=2,method="GA",data=myData,orth=FALSE,CV="boot",
             varList=rrr,criterion="RMSEA",minInd=3,niter=1)
uu
fits = c(uu[[1]]$fit,uu[[2]]$fit)
which(min(fits)==fits)
