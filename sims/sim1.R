library(parallel)

mod1 <- '
f1 =~ 1*x1 + 1*x2 + 1*x3 + 1*x4 + 0*x5 + 0*x6 + 0*x7 + 0.3*x8 + 0*x9
f1~~1*f1
'


mod2 <- '
f1 =~ 1*x1 + 1*x2 + 1*x3 + 1*x4 + 1*x5
f2 =~ 1*x6 + 1*x7 + 1*x8 + 1*x9
f1~~1*f1
f2~~1*f2
f1~~0*f2
'

mod3 <- '
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

crit = c("BIC","RMSEA")
cv = c(FALSE,TRUE,"boot")



f1.vars <- c("x1","x2","x3","x4","x5","x6","x7","x8","x9")
rrr = list(f1.vars)
facs <- 1:4
iters=20
count=0
mat = matrix(NA,iters*3*2*3,4)
mods = c(mod1,mod2,mod3)

for(i in 1:iters){
  for(j in 1:length(mods)){
    for(o in 1:length(crit)){
      for(p in 1:length(cv)){
        count=count+1

        if(crit[o]=="BIC" & cv[p] == TRUE | crit[o]=="BIC" & cv[p] == "boot"){

          mat[count,1] = NA
          mat[count,2] = NA
          mat[count,3] = NA
          mat[count,4] = NA
        }else{


    myData = simulateData(mods[[j]],model.type="cfa",sample.nobs=500)


    uu = multFac(facList=facs,parallel="yes",ncore=4,method="tabuSearch",data=myData,orth=FALSE,CV=cv[p],
                 varList=rrr,criterion=crit[o],minInd=3,niter=30)

    fits = c(uu[[1]]$fit,uu[[2]]$fit,uu[[3]]$fit,uu[[4]]$fit)
    mat[count,1] = j
    mat[count,2] = paste(which(min(fits)==fits),collapse="")
    mat[count,3] = crit[o]
    mat[count,4] = cv[p]

        }
      }
    }
  }
  print(i)
}


mat
mat2 = na.omit(mat)
res = matrix(NA,nrow(mat2),2)
for(j in 1:nrow(mat2)){
  res[j,1] = mat2[j,1] == mat2[j,2]

  res[j,2] = paste(c(mat2[j,3],mat2[j,4]),collapse="")
}

res1 = res[res[,2] == "BICFALSE",]
sum(res1[,1]== "TRUE")/60

res2 = res[res[,2] == "RMSEAFALSE",]
sum(res2[,1]== "TRUE")/60

res3 = res[res[,2] == "RMSEATRUE",]
sum(res3[,1]== "TRUE")/60

res4 = res[res[,2] == "RMSEAboot",]
sum(res4[,1]== "TRUE")/60


