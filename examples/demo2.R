library(lavaan)

N=1000
population.model <- '
f1 =~ 1*x1 + 1*x2 + 1*x3 + 1*x4 + 1*x5
f2 =~ 1*x6 + 1*x7 + 1*x8 + 1*x9
f1~~1*f1
f2~~1*f2
f1~~0*f2
'


population.model <- '
f1 =~ 1*x1 + 1*x2 + 1*x3 + 1*x4 + 0*x5 + 0*x6 + 0*x7 + 0.3*x8 + 0*x9
f1~~1*f1
'


myData <- simulateData(population.model, sample.nobs=N)
myData.test <- simulateData(population.model, sample.nobs=N)


fa <-'
f1 =~ NA*x1 + x2 + x3 + x4 + x5 + x6 + x7 + x8 + x9
f1~~1*f1
'

fa <-'
f1 =~ NA*x1 + x2
f1~~1*f1
'
options(warn=1)
sim.out = cfa(fa,myData)
summary(sim.out)
1/fitMeasures(sim.out)["bic"]

start= c(1,1,1,0,0,0,1,1,1)
nfac=1





f1.vars <- c("x1","x2","x3","x4","x5","x6","x7","x8","x9")
rrr = list(f1.vars)
facs <- 1:4

ret.auto <- function(facs){
    res <- list()
    ret = autoSEM(method="GA",data=myData,nfac=facs,orth=FALSE,replaceSamp=FALSE,boot=TRUE,
                  varList=rrr,criterion="RMSEA",minInd=3,niter=3)
    res$fit = 1/(summary(ret)$fitness-1)
    #res$fit = -(summary(ret)$fitness-1)
    res$solution = summary(ret)$solution
    res
}

library(snowfall)
sfStop()
sfInit(T,4)
snowfall::sfExport("myData","rrr","facs")
sfLibrary(autoSEM); sfLibrary(lavaan); sfLibrary(GA); sfLibrary(tabuSearch)
(out = sfLapply(facs,ret.auto))

