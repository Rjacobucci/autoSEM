library(lavaan)

N=5000
population.model <- '
f1 =~ 1*x1 + 1*x2 + 1*x3 + 1*x4 + 1*x5
f2 =~ 0.3*x6 + 0.3*x7 + 0.2*x8 + 0.1*x9
f1~~1*f1
f2~~1*f2
f1~~0*f2
'


myData <- simulateData(population.model, sample.nobs=N)
myData.test <- simulateData(population.model, sample.nobs=N)


fa <-'
f1 =~ NA*x1 + x2 + x3 + x4 + x5 + x6 + x7 + x8 + x9
f1~~1*f1
'
sim.out = cfa(fa,myData)
summary(sim.out)

nfac=1
f1.vars <- c("x1","x2","x3","x4","x5","x6","x7","x8","x9")
rrr = list(f1.vars)


ret = autoSEM(method="GA",data=myData,nfac=1,
              varList=rrr,criterion="BIC",minInd=3,niter=20)


ret$configKeep[max(ret$eUtilityKeep)==ret$eUtilityKeep,]
