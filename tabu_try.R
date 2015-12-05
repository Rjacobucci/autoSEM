library(lavaan)

mod = '
f1 =~ 1*x1 + 1*x2 + 1*x3 + 1*x4 + 1*x5 + 1*x6 + 0*x7 + 0*x8 + 0*x9
f1~~1*f1'

myData = simulateData(mod,sample.nobs=300,model.type="cfa")

# have to load just string
f1.vars <- c("x1","x2","x3","x4","x5","x6","x7","x8","x9")
rrr = list(f1.vars)
#facList=facs
nfac=1
parallel="no"
method="tabu_rj"
data=myData
minInd=3
stdlv=TRUE
orth=TRUE
orth=FALSE
CV=FALSE
varList=rrr
criterion="BIC"
niter=30
length=9
string = matrix(NA,1,length)
for(i in 1:9) string[i] = rbinom(1,1,prob=0.5)
fitness(string)
fitness(new.taboo[300,])





tabu_rj <- function(string){

  out = list()
  K.best = 999999999
  best.string <- rep(NA,length(string))
  taboo.list <- matrix(NA,100,length(string))
  count=0

while(count < n.iter){

  count=count+1

  string.mat = matrix(NA,length(string),length(string))

  for(i in 1:nrow(string.mat)){
    new.string=string
    new.string[i] = ifelse(new.string[i] == 1,0,1)
    string.mat[i,] = new.string
  }


  rowmatch <- function(A,B) {
    # Rows in A that match the rows in B
    f <- function(...) paste(..., sep=":")
    if(!is.matrix(B)) B <- matrix(B, 1, length(B))
    a <- do.call("f", as.data.frame(A))
    b <- do.call("f", as.data.frame(B))
    match(b, a)
  }


  rww = rowmatch(taboo.list,string.mat)
  rww2 = rww[is.na(rww)==FALSE]
  if(length(rww2)==0){
    string.mat2 = string.mat
  }else{
    string.mat2 = string.mat[-rww2,]
  }

  if(nrow(string.mat2)==0){
    break
  }




  K.vec = rep(NA,nrow(string.mat2))
  for(i in 1:nrow(string.mat2)){
    K.vec[i] = fitness(string.mat2[i,])
  }

  pos = which(K.vec == min(K.vec))

  if(min(K.vec) < K.best){
    K.best=min(K.vec)
    best.string = string.mat2[pos,]
  }

  new.taboo = rbind(taboo.list,string.mat2)
  if(sum(is.na(taboo.list[,1])) >0){
    new.taboo = new.taboo[-c(1:sum(is.na(taboo.list[,1]))),]
  }

  if(nrow(new.taboo)>100){
    new.taboo2 = new.taboo[(nrow(new.taboo)-100):nrow(new.taboo),]
  }
  taboo.list=new.taboo2

  if(length(pos)>1){
    string=string.mat[1,]
  }else{
    string=string.mat[pos,]
  }

}
out$solution = best.string
out$value = K.best

out

}

