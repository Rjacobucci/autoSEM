
caretGA <- list(
  fit = function(x, y, ...) train(x, y, "qda", 
                                  metric = "ROC",
                                  trControl = trainControl(method = "repeatedcv", 
                                                           repeats = 1,
                                                           summaryFunction = twoClassSummary,
                                                           classProbs = TRUE)),
  fitness = function(object, x, y)
  {
    caret:::getTrainPerf(object)[, paste("Train", object$metric, sep = "")]
  },
  predict = function(object, x)
  {
    library(caret)
    predict(object, newdata = x)
  }
)


x2index <- function(x)
{  
  binary <- binomial()$linkinv(x)
  binary <- ifelse(binary >= .5, 1, 0) 
  apply(binary, 2, function(x) which(x == 1))
}

psofs <- function(x, y, iterations = 10, functions = caretGA, 
                  maximize = TRUE, verbose = FALSE, 
                  parallel = TRUE, 
                  convert = x2index, tx = NULL, ty = NULL, ...)
{
  call <- match.call()
  varNames <- colnames(x)
  numVar <- ncol(x)
  par <- rep(0, numVar)
  par[sample(seq(along = par), 10)] <- 1
  
  psoObj <- function(vars, x, y, tx, ty, maximize, func, ...)
  {
    fit <- func$fit(x[, vars, drop = FALSE], y, ...)
    out <- func$fitness(fit, tx[, vars, drop = FALSE], ty)
    out 
  }
  
  results <- binary_pso(par, 
                        fn = psoObj, 
                        parallel = parallel,
                        x = x,
                        y = y,
                        convert = convert,
                        tx = tx,
                        ty = ty,
                        maximize = maximize,
                        func = functions,
                        control = list(maxit = iterations),
                        ...)
  ## fit final model
  results$varNames <- varNames
  results$call <- call
  results$maximize <- maximize
  results$func <- functions
  results$fit <- functions$fit(x[, results$bestVars[[1]], drop = FALSE], y)
  class(results) <- "psofs"
  results
}

print.psofs <- function(x, top = 5, digits = max(3, getOption("digits") - 3), ...)
{
  topx <- x$par[1:min(top, length(x$pa))]
  topx <- x$varNames[topx]
  topx <- paste(topx, collapse = ", ")
  cat("\nParticle Swarm Feature Selection\n\n")
  cat("Iterations: ", x$iterations,
      "\nSwarm Size: ", x$numSwarm,
      "\nAv. Informants:", signif(x$numInform, digits),
      "\nExploitation Value:", signif(x$exploitation, digits),
      "\n\n")
  cat("Best fitness value: ", signif(x$fitness[order(-x$fitness[,1]),][1,1], digits),
      "\n", length(x$par), " Predictors (",
      topx, ifelse(top < length(x$varNames), ", ...)", ")"), "\n",
      sep = "")
  ## mean swarm size and fitness at end?
  
  
}


## This function is a modifed version of the psoptim funciton in the pso package. 
binary_pso <- function (par, fn, ..., 
                        convert = x2index,
                        maximize = TRUE,
                        parallel = FALSE,
                        verbose = FALSE,
                        control = list()) {
  lower <- -1 
  upper <- 1

  fn1 <- function(par) fn(par, ...)/p.fnscale
  mrunif <- function(n,m,lower,upper) {
    return(matrix(runif(n*m,0,1),nrow=n,ncol=m)*(upper-lower)+lower)
  }
  norm <- function(x) sqrt(sum(x*x))
  rsphere.unif <- function(n,r) {
    temp <- runif(n)
    return((runif(1,min=0,max=r)/norm(temp))*temp)
  }
  svect <- function(a,b,n,k) {
    temp <- rep(a,n)
    temp[k] <- b
    return(temp)
  }
  mrsphere.unif <- function(n,r) {
    m <- length(r)
    temp <- matrix(runif(n*m),n,m)
    return(temp%*%diag(runif(m,min=0,max=r)/apply(temp,2,norm)))
  }
  npar <- length(par)
  lower <- as.double(rep(lower, ,npar))
  upper <- as.double(rep(upper, ,npar))
  con <- list(trace = 0, fnscale = 1, maxit = 1000L, maxf = Inf,
              abstol = -Inf, reltol = 0, REPORT = 10,
              s = NA, k = 3, p = NA, w = 1/(2*log(2)),
              c.p = .5+log(2), c.g = .5+log(2), d = NA,
              v.max = NA, rand.order = TRUE, max.restart=Inf,
              maxit.stagnate = Inf,
              vectorize=FALSE, 
              trace.stats = FALSE, type = "SPSO2007")
  nmsC <- names(con)
  con[(namc <- names(control))] <- control
  if (length(noNms <- namc[!namc %in% nmsC])) 
    warning("unknown names in control: ", paste(noNms, collapse = ", "))
  ## Argument error checks
  if (any(upper==Inf | lower==-Inf))
    stop("fixed bounds must be provided")
  
  p.type <- pmatch(con[["type"]],c("SPSO2007","SPSO2011"))-1
  if (is.na(p.type)) stop("type should be one of \"SPSO2007\", \"SPSO2011\"")
  
  p.trace <- con[["trace"]]>0L # provide output on progress?
  p.fnscale <- con[["fnscale"]] # scale funcion by 1/fnscale
  p.maxit <- con[["maxit"]] # maximal number of iterations
  p.maxf <- con[["maxf"]] # maximal number of function evaluations
  p.abstol <- con[["abstol"]] # absolute tolerance for convergence
  p.reltol <- con[["reltol"]] # relative minimal tolerance for restarting
  p.report <- as.integer(con[["REPORT"]]) # output every REPORT iterations
  p.s <- ifelse(is.na(con[["s"]]),ifelse(p.type==0,floor(10+2*sqrt(npar)),40),
                con[["s"]]) # swarm size
  p.p <- ifelse(is.na(con[["p"]]),1-(1-1/p.s)^con[["k"]],con[["p"]]) # average % of informants
  p.w0 <- con[["w"]] # exploitation constant
  if (length(p.w0)>1) {
    p.w1 <- p.w0[2]
    p.w0 <- p.w0[1]
  } else {
    p.w1 <- p.w0
  }
  p.c.p <- con[["c.p"]] # local exploration constant
  p.c.g <- con[["c.g"]] # global exploration constant
  p.d <- ifelse(is.na(con[["d"]]),norm(upper-lower),con[["d"]]) # domain diameter
  p.vmax <- con[["v.max"]]*p.d # maximal velocity
  p.randorder <- as.logical(con[["rand.order"]]) # process particles in random order?
  p.maxrestart <- con[["max.restart"]] # maximal number of restarts
  p.maxstagnate <- con[["maxit.stagnate"]] # maximal number of iterations without improvement
  p.vectorize <- as.logical(con[["vectorize"]]) # vectorize?
  p.trace.stats <- as.logical(con[["trace.stats"]]) # collect detailed stats?
  
  if (p.trace) {
    message("S=",p.s,", K=",con[["k"]],", p=",signif(p.p,4),", w0=",
            signif(p.w0,4),", w1=",
            signif(p.w1,4),", c.p=",signif(p.c.p,4),
            ", c.g=",signif(p.c.g,4))
    message("v.max=",signif(con[["v.max"]],4),
            ", d=",signif(p.d,4))
    if (p.trace.stats) {
      stats.trace.it <- c()
      stats.trace.error <- c()
      stats.trace.f <- NULL
      stats.trace.x <- NULL
    }
  }
  ## Initialization
  outX <- outY <- vector(mode = "list", length = p.maxit)
  charIt <- format(1:p.maxit)
  if (p.reltol!=0) p.reltol <- p.reltol*p.d
  
  lowerM <- matrix(lower,nrow=npar,ncol=p.s)
  upperM <- matrix(upper,nrow=npar,ncol=p.s)
  
  X <- mrunif(npar,p.s,lower,upper)
  if (!any(is.na(par)) && all(par>=lower) && all(par<=upper)) X[,1] <- par
  if (p.type==0) {
    V <- (mrunif(npar,p.s,lower,upper)-X)/2
  } else { ## p.type==1
    V <- matrix(runif(npar*p.s,min=as.vector(lower-X),max=as.vector(upper-X)),npar,p.s)
    p.c.p2 <- p.c.p/2 # precompute constants
    p.c.p3 <- p.c.p/3
    p.c.g3 <- p.c.g/3
    p.c.pg3 <- p.c.p3+p.c.g3
  }
  if (!is.na(p.vmax)) { # scale to maximal velocity
    temp <- apply(V,2,norm)
    temp <- pmin.int(temp,p.vmax)/temp
    V <- V%*%diag(temp)
  }

  # first evaluations
  index <- convert(X)
  if(length(index) == 0) index <- list()
  if(parallel)
  {
    library(foreach)
    f.x <- foreach(part = seq(along = index), .combine = "rbind") %dopar% fn(index[[part]], ...)
    f.x <- as.data.frame(f.x)
    rownames(f.x) <- paste(1:nrow(f.x))
  } else {
    f.x <- lapply(index, fn, ...) 
    f.x <- do.call("rbind", f.x)
    f.x <- as.data.frame(f.x)
  }
  f.x$iter <- 1
  outX[[1]] <- index
  outY[[1]] <- f.x 
  f.x <- f.x[,1]
  if(maximize) f.x <- -f.x

  stats.feval <- p.s
  P <- X
  f.p <- f.x
  P.improved <- rep(FALSE,p.s)
  i.best <- which.min(f.p)
  error <- f.p[i.best]
  init.links <- TRUE
  if (verbose) {
    message(charIt[1], ": best = ",
            ifelse(maximize, -signif(error,4), signif(error,4)),
            " mean = ",
            ifelse(maximize, -signif(mean(f.p),4), signif(mean(f.p),4)),
            " av size = ",
            signif(mean(unlist(lapply(index, length))),4))
    if (p.trace.stats) {
      stats.trace.it <- c(stats.trace.it,1)
      stats.trace.error <- c(stats.trace.error,error)
      stats.trace.f <- c(stats.trace.f,list(f.x))
      stats.trace.x <- c(stats.trace.x,list(X))
    }
  }
  ## Iterations
  stats.iter <- 1
  stats.restart <- 0
  stats.stagnate <- 0
  
  while (stats.iter<p.maxit && stats.feval<p.maxf && error>p.abstol &&
           stats.restart<p.maxrestart && stats.stagnate<p.maxstagnate) {
    stats.iter <- stats.iter+1
    if (p.p!=1 && init.links) {
      links <- matrix(runif(p.s*p.s,0,1)<=p.p,p.s,p.s)
      diag(links) <- TRUE
    }
    ## The swarm moves
    
    if (p.p==1) j <- rep(i.best,p.s)
    else # best informant
      j <- sapply(1:p.s,function(i)
        which(links[,i])[which.min(f.p[links[,i]])]) 
    temp <- (p.w0+(p.w1-p.w0)*max(stats.iter/p.maxit,stats.feval/p.maxf))
    V <- temp*V # exploration tendency
    if (p.type==0) {
      V <- V+mrunif(npar,p.s,0,p.c.p)*(P-X) # exploitation
      temp <- j!=(1:p.s)
      V[,temp] <- V[,temp]+mrunif(npar,sum(temp),0,p.c.p)*(P[,j[temp]]-X[,temp])
    } else { # SPSO 2011
      temp <- j==(1:p.s)
      temp <- P%*%diag(svect(p.c.p3,p.c.p2,p.s,temp))+
        P[,j]%*%diag(svect(p.c.g3,0,p.s,temp))-
        X%*%diag(svect(p.c.pg3,p.c.p2,p.s,temp)) # G-X
      V <- V+temp+mrsphere.unif(npar,apply(temp,2,norm))
    }
    if (!is.na(p.vmax)) {
      temp <- apply(V,2,norm)
      temp <- pmin.int(temp,p.vmax)/temp
      V <- V%*%diag(temp)
    }
    X <- X+V
    ## Check bounds
    temp <- X<lowerM
    if (any(temp)) {
      X[temp] <- lowerM[temp] 
      V[temp] <- 0
    }
    temp <- X>upperM
    if (any(temp)) {
      X[temp] <- upperM[temp]
      V[temp] <- 0
    }
    ## Evaluate function
    index <- convert(X)
    if(length(index) == 0) index <- list()
    if(parallel)
    {
      library(foreach)
      f.x <- foreach(part = seq(along = index), .combine = "rbind") %dopar% fn(index[[part]], ...)
      f.x <- as.data.frame(f.x)
      rownames(f.x) <- paste(1:nrow(f.x))
    } else {
      f.x <- lapply(index, fn, ...) 
      f.x <- do.call("rbind", f.x)
      f.x <- as.data.frame(f.x)
    }
    f.x$iter <- stats.iter
    outX[[stats.iter]] <- index
    outY[[stats.iter]] <- f.x 
    f.x <- f.x[,1]
    if(maximize) f.x <- -f.x
    stats.feval <- stats.feval+p.s
    
    temp <- f.x<f.p
    if (any(temp)) { # improvement
      P[,temp] <- X[,temp]
      f.p[temp] <- f.x[temp]
      i.best <- which.min(f.p)
    }
    if (stats.feval>=p.maxf) break
    
    if (p.reltol!=0) {
      d <- X-P[,i.best]
      d <- sqrt(max(colSums(d*d)))
      if (d<p.reltol) {
        X <- mrunif(npar,p.s,lower,upper)
        V <- (mrunif(npar,p.s,lower,upper)-X)/2
        if (!is.na(p.vmax)) {
          temp <- apply(V,2,norm)
          temp <- pmin.int(temp,p.vmax)/temp
          V <- V%*%diag(temp)
        }
        stats.restart <- stats.restart+1
        if (verbose) message(charIt[stats.iter], ": restarting")
      }
    }
    init.links <- f.p[i.best]==error # if no overall improvement
    stats.stagnate <- ifelse(init.links,stats.stagnate+1,0)
    error <- f.p[i.best]
    if (verbose) {
      if (p.reltol!=0) 
        message(charIt[stats.iter], ": best = ",
                ifelse(maximize, -signif(error,4), signif(error,4)),
                ", swarm diam.=",signif(d,4))
      else
        message(charIt[stats.iter], ": best = ",
                ifelse(maximize, -signif(error,4), signif(error,4)),
                " mean = ",
                ifelse(maximize, -signif(mean(f.p),4), signif(mean(f.p),4)),
                " av size = ",
                signif(mean(unlist(lapply(index, length))),4))
      if (p.trace.stats) {
        stats.trace.it <- c(stats.trace.it,stats.iter)
        stats.trace.error <- c(stats.trace.error,error)
        stats.trace.f <- c(stats.trace.f,list(f.x))
        stats.trace.x <- c(stats.trace.x,list(X))
      }
    }
  }
  if (error<=p.abstol) {
    msg <- "Converged"
    msgcode <- 0
  } else if (stats.feval>=p.maxf) {
    msg <- "Maximal number of function evaluations reached"
    msgcode <- 1
  } else if (stats.iter>=p.maxit) {
    msg <- "Maximal number of iterations reached"
    msgcode <- 2
  } else if (stats.restart>=p.maxrestart) {
    msg <- "Maximal number of restarts reached"
    msgcode <- 3
  } else {
    msg <- "Maximal number of iterations without improvement reached"
    msgcode <- 4
  }
  if (verbose) message(msg)
  part <- unlist(lapply(outY, function(x) 1:nrow(x)))
  fitness <- do.call("rbind", outY)
  fitness$particle <- part
  ord <- if(maximize) order(-fitness[,1]) else order(fitness[,1])
  top <- 5
  bestFit <- head(fitness[ord,], top)
  bestVars <- vector(mode = "list", length = top)
  for(i in 1:top)
  {
    bestVars[[i]] <- outX[[bestFit$iter[i]]][[bestFit$particle[i]]]
  }
  
  
  vImp <- table(unlist(outX))/nrow(fitness)*100
  vImp <- data.frame(Variable = names(training)[as.integer(names(vImp))],
                     Percent = as.vector(vImp))
  vImp <- vImp[order(-vImp$Percent),]
  ## dotplot(reorder(Variable, Percent) ~ Percent, vImp[1:30, ])

  o <- list(par=convert(P[,i.best, drop = FALSE]),
            value= if(maximize) -f.p[i.best] else f.p[i.best],
            counts=c("function"=stats.feval,"iteration"=stats.iter,
                     "restarts"=stats.restart),
            convergence=msgcode,
            message=msg,
            options = con,
            exploitation = p.w0,
            numSwarm = p.s,
            numInform = p.p,
            iterations = p.maxit,
            index = outX,
            selectionPct  = vImp,
            bestVars = bestVars,
            bestFit = bestFit,
            fitness = fitness)
  return(o)
}
