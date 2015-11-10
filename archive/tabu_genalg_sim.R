library(lavaan)

N=500
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

library(GA)


pars = parTable(sim.out)


f1.vars <- c("x1","x2","x3","x4","x5","x6")
f2.vars <- c("x4","x5","x6","x7","x8","x9")
f3.vars <- c("x1","x2","x3","x7","x8","x9")

c(1,1,1,0,0,0)
f1.vars[c(1,1,1,0,0,0) == 1]

(c(1,1,1,0,0,0) == 1)
mod <- lm(Y ~ ., diab.train)

x <- model.matrix(mod)[, -1]
y <- model.response(model.frame(mod))

start = c(1,1,1,0,0,0,1,1,1,0,0,0,0,0,0,1,1,1)

#cov.order = out@Data@ov.names
#cov.train = cov(myData[,unlist(cov.order)])







fitness <- function(string) {
  pt1 = string[1:6]
  pt2 = string[7:12]
  pt3 = string[13:18]

  if(sum(pt1) < 3 | sum(pt2) < 3 | sum(pt3) < 3) return(1) #return(-99999999999)

  inc1 <- which(pt1 == 1)
  inc2 <- which(pt2 == 1)
  inc3 <- which(pt3 == 1)

  f1.vars <- c("x1","x2","x3","x4","x5","x6")
  f2.vars <- c("x4","x5","x6","x7","x8","x9")
  f3.vars <- c("x1","x2","x3","x7","x8","x9")

  lll = list(f1.vars,f2.vars,f3.vars)

  fmld <- c(paste("factor1 =~ ", paste(f1.vars[inc1], collapse= "+")),
    paste("factor2 =~ ", paste(f2.vars[inc2], collapse= "+")),
    paste("factor3 =~ ", paste(f3.vars[inc3], collapse= "+")))


  #vars = unique(c(f1.vars[inc1],f2.vars[inc2],f3.vars[inc3]))
  #cov.test = cov(myData.test[vars])


  out = lavaan::cfa(fmld,myData,orthogonal=T,std.lv=T)

  if(inspect(out,"converged")==F | any(eigen(inspect(out,"cov.lv"))$values < 0)){
    #return(-99999999999999)
    return(1)
  }else{
    #fitMeasures(out)["bic"]
    df=out@Fit@test[[1]]$df
    cov.order = out@Data@ov.names
    cov.test = cov(myData.test[,unlist(cov.order)])
    impcov = fitted(out)$cov
    fit.test = 0.5*(log(det(impcov)) + trace(cov.test %*% solve(impcov)) - log(det(cov.test))  - 9)
    chisq.test = N*fit.test
    ncp.test = d(chisq.test,df,N)
    RMSEA.test = rmsea(ncp.test,df)
    -RMSEA.test + 100
  }
}

#GA <- ga("binary", fitness = fitness, nBits = 18,monitor=T)
library(tabuSearch)
result <- tabuSearch(size = 18, iters = 50,objFunc = fitness,config=start)

pss = which(max(result$eUtilityKeep) == result$eUtilityKeep)
result$configKeep[pss,]

plot(GA)
summary(GA)
summary(GA)$solution
