
tabu_rj <- function(size, iters,fitness,min.improve,seed){

  string = matrix(NA,1,size)
  set.seed(seed)
  for(i in 1:size) string[i] = rbinom(1,1,prob=0.5)

  n.models.run = 0
  n.models.run2 = 0
  iter.no.improve = 0
  out = list()
  K.best = 999999999
  best.string <- rep(NA,length(string))
  taboo.list <- matrix(NA,length(string)*5,length(string))
  count=0

while(count < iters){

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
    stop("no matches")
  }




  K.vec = rep(NA,nrow(string.mat2))
  for(i in 1:nrow(string.mat2)){
    K.vec[i] = fitness(string.mat2[i,])
  }

  n.models.run = n.models.run + sum(is.na(K.vec) == FALSE)
  n.models.run2 = n.models.run2 + sum(K.vec > -1 )

  pos = which(K.vec == min(K.vec))

  if(min(K.vec) < K.best){
    K.best=min(K.vec)
    best.string = string.mat2[pos,]
    iter.no.improve = 0
  }else{
    iter.no.improve = iter.no.improve + 1
  }

  if(iter.no.improve == min.improve){
    break
  }

  taboo.list = rbind(taboo.list[length(string):(length(string)*5),],string.mat)
#  if(sum(is.na(taboo.list[,1])) >0){
  #    #new.taboo = new.taboo[-c(1:sum(is.na(taboo.list[,1]))),]
  #    taboo.list=new.taboo
  #  }

  #  if(nrow(taboo.list)>100){
  #    new.taboo2 = new.taboo[(nrow(new.taboo)-100):nrow(new.taboo),]
  #    taboo.list=new.taboo2
  #  }

  if(length(pos)>1){
    string=string.mat2[pos[1],]
  }else{
    string=string.mat2[pos,]
  }

}
out$solution = best.string
out$value = K.best
out$n.models.run = n.models.run
out$n.models.run2 = n.models.run2
#out$tabu_list = taboo.list
out$count = count
out$iter.no.improve = iter.no.improve
out

}

