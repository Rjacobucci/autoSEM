#library(devtools)

#install_github("RJacobucci/autoSEM")
#library(autoSEM)




library(lavaan)

mod = '
f1 =~ 1*x1 + 1*x2 + 1*x3 + 1*x4 + 1*x5 + 1*x6 + 0*x7 + 0*x8 + 0*x9
f1~~1*f1'

myData = simulateData(mod,sample.nobs=300,model.type="cfa")

f1.vars <- c("x1","x2","x3","x4","x5","x6","x7","x8","x9")
rrr = list(f1.vars)
facs <- 1:4


uu = multFac(facList=facs,parallel="yes",ncore=4,method="GA",
             data=myData,orth=FALSE,CV=FALSE,
             varList=rrr,criterion="RMSEA",niter=30)
uu

summary(uu[[1]]$out)
print(uu[[1]]$out)
max(uu[[1]]$out$eUtilityKeep)

fits = c(uu[[1]]$fit,uu[[2]]$fit)
which(min(fits)==fits)


out = autoSEM(method="aco_rj",data=myData,nfac=1,varList=list(f1.vars),orth=FALSE,CV=FALSE,
        criterion="BIC",minInd=3,aco.iters=100)
out
round(out$out$samp.wgt,3)
