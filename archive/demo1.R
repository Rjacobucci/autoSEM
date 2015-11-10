library(devtools)

install_github("RJacobucci/autoSEM")
library(autoSEM)
library(tabuSearch)
library(GA)


library(lavaan)

N=5000
population.model <- '
f1 =~ 1*x1 + 0.2*x2 + 1*x3
f2 =~ 1*x4 + 0.2*x5 + 1*x6
f3 =~ 1*x7 + 0.2*x8 + 1*x9
f1~~1*f1
f2~~1*f2
f3~~1*f3
f1~~0*f2
f2~~0*f3
f1~~0*f3
'


myData <- simulateData(population.model, sample.nobs=N)
myData.test <- simulateData(population.model, sample.nobs=N)


fa <-'
f1 =~ 1*x1 + x2 + x3 + x4 + x5 + x6
f2 =~ 1*x4 + x5 + x6 + x7 + x8 + x9
f3 =~ NA*x1 + x2 + x3 + 1*x7 + x8 + x9
f1~~1*f1
f2~~1*f2
f3~~1*f3
f1~~0*f2
f2~~0*f3
f1~~0*f3
'
sim.out = cfa(fa,myData)

nfac=3
f1.vars <- c("x1","x2","x3","x4","x5","x6")
f2.vars <- c("x4","x5","x6","x7","x8","x9")
f3.vars <- c("x1","x2","x3","x7","x8","x9")
rrr = list(f1.vars,f2.vars,f3.vars)


ret = autoSEM(method="tabuSearch",data=myData,nfac=3,varList=rrr,criterion="RMSEA",minInd=3,niter=20)


ret$configKeep[max(ret$eUtilityKeep)==ret$eUtilityKeep,]
