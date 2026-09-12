require(lpSolveAPI)

predict.PVM <- function(object, ..., newdata, eta="adjtogrid", probepsilon="adjtogrid", trndt=NULL, retallprd=FALSE, newdataasKmat=FALSE)
{
  if (probepsilon=="adjtogrid") probepsilon <- min(object$grid)/2
  if ( !is.numeric(probepsilon) || length(probepsilon) > 1 || probepsilon[1] < 0 || probepsilon[1] > 1)
    stop("Wrong value for argument probepsilon\n")

  if (!is.matrix(newdata) && !is.data.frame(newdata))
  {
    if (!is.numeric(newdata) && is.integer(newdata)) stop("Wrong type for argument newdata\n")
    else newdata <- matrix(newdata,nrow=1)
  }
  n <- nrow(newdata)
  npoints <- length(object$svms)
  k <- object$nclasses
  if (eta=="adjtogrid") eta <- k*nrow(object$grid)
  Pprob   <- matrix(nrow=n,ncol=k)
  if (!newdataasKmat) {

    if (is.data.frame(newdata) || is.matrix(newdata)) {
      if (!is.null(object$center)) {
        if (!is.null(object$scale)) newdata <- scale(newdata,center=object$center,scale=object$scale)
        else newdata <- scale(newdata,center=object$center,scale=FALSE)
      } else if (!is.null(object$scale)) newdata <- scale(newdata,center=FALSE,scale=object$scale)
      if (is.data.frame(newdata))  newdata <- as.matrix(newdata)
    } else if (is.numeric(newdata) || is.integer(newdata)) {
      if (!is.null(object$center)) {
        if (!is.null(object$scale)) newdata <- (newdata-object$center)/object$scale
        else newdata <- newdata-object$center
      } else if (!is.null(object$scale)) newdata <- newdata/object$scale
    }

    if (is.character(trndt) && trndt=="newdata") {
      K <- matrix(nrow=n,ncol=n)
      K[lower.tri(K,diag=TRUE)] <- makeKMat0(newdata, kernel=object$kernel, kpar=object$kpar)
      K[upper.tri(K,diag=FALSE)] <- t(K)[upper.tri(K,diag=FALSE)]

    } else {
     if (!is.null(object$x) && !is.null(trndt))
       warning("The value of argument trndt will be ignored, since the training data was already suplied by the x component of the value of argument object\n")
     if (is.null(object$x)) {
       if (is.null(trndt))
         stop("Argument trndt (training data) needs to be provided when argument object does not have an x component and argument newdataasKmat is not set to TRUE\n")

       if ( !is.null(dim(trndt)) ) {
          if (!is.null(object$center)) {
            if (!is.null(object$scale))  trndt <- scale(trndt,center=object$center,scale=object$scale)
            else trndt <- scale(trndt,center=object$center,scale=FALSE)
          } else {
            if (!is.null(object$scale)) trndt <- scale(trndt,center=FALSE,scale=object$scale)
          }
       }
       if (is.data.frame(newdata) || is.matrix(newdata)) K <- makeKMat0(trndt, newdata, kernel=object$kernel, kpar=object$kpar)
       else K <- makeKMat0(trndt,  t(as.matrix(newdata)), kernel=object$kernel, kpar=object$kpar)
      } else {
        if (is.data.frame(newdata) || is.matrix(newdata)) K <- makeKMat0(object$x, newdata, kernel=object$kernel, kpar=object$kpar)
        else K <- makeKMat0(object$x,  t(as.matrix(newdata)), kernel=object$kernel, kpar=object$kpar)
      }
    }
    if (class(object$svms)[1]!="kernelksvms") classpred <- predict(object$svms,newdata=NULL,KMat=K)
    else classpred <- predict(object$svms,newdata=NULL,KMat=K)

  } else {
    if (!is.matrix(newdata))  stop("Wrong value for argument newdata\n")
    if (!is.null(trndt))
      warning("The value of argument trndt is ignored when the new data is provided as an already computed kernel matrix\n")
    if (class(object$svms)[1]!="kernelksvms") classpred <- predict(object$svms,newdata=NULL,KMat=newdata)
    else classpred <- predict(object$svms,newdata=NULL,KMat=newdata)
  }

  if (is.factor(classpred)) {
    if (k>2) alpha <- LPdual(classpred,object$grid,eta,probepsilon)
    else alpha <- k2Prob(classpred,object$grid[,1],eta)
  } else if (is.data.frame(classpred)) {
     for (obs in 1:n) {
       if (k>2) alpha <- LPdual(classpred[obs,],object$grid,eta,probepsilon)
       else alpha <- k2Prob(classpred[obs,,drop=FALSE],object$grid[,1],eta)
       if (is.null(alpha)) warning("Class probability estimation failed for observation",obs,"\n")
       else Pprob[obs,] <- alpha
     }
  } else stop("Wrong class for classpred object\n")
  if (!is.null(rownames(newdata))) rownames(Pprob) <- rownames(newdata)
  else rownames(Pprob) <- 1:nrow(newdata)
  colnames(Pprob) <- object$svms$grplvls

  class(Pprob) <- "ClassProb"

  if (retallprd) {
    return(list(Pprob=Pprob,classpred=classpred))
  }

  Pprob
}

k2Prob <- function(Yhat,PI1g,epsilon=PI1g[1]/2,eta=2*length(PI1g))
{
  Yhat <- unclass(Yhat)
  PI2g <- 1 - PI1g
  p1 <- PI1g
  p2 <- 1- p1
  f <- 0.
  for (i in 1:length(PI1g)) {
    if (Yhat[i]==1) negdif <- PI2g[i]*p2 - PI1g[i]*p1
    else if (Yhat[i]==2) negdif <- PI1g[i]*p1 - PI2g[i]*p2
    f <- f + ifelse(negdif>0,eta*negdif,negdif)
  }
  p1hat <- max(PI1g[which.min(f)],epsilon)
  c(p1hat,1-p1hat)
}

