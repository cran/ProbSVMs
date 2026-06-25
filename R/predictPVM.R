require(lpSolveAPI)

predict.PVM <- function(object, ..., newdata, eta=15., probepsilon="adjtogrid", trndt=NULL, retallprd=FALSE, newdataasKmat=FALSE)
{
  if (probepsilon=="adjtogrid") probepsilon <- min(object$grid)/2
  if ( !is.numeric(probepsilon) || length(probepsilon) > 1 || probepsilon[1] < 0 || probepsilon[1] > 1)
    stop("Wrong value for argument probepsilon\n")

  n <- nrow(newdata)
  npoints <- length(object$svms)
  k <- object$nclasses
  Pprob   <- matrix(nrow=n,ncol=k)
  if (!newdataasKmat) {
    if (is.null(trndt)) stop("Argument trndt (training data) needs to be provided when argument newdataasKmat is not set to TRUE\n")
    if (!is.null(dim(trndt)) || trndt!="newdata") K <- makeKMat0(trndt, newdata, kernel=object$kernel, kpar=object$kpar)
    else {
      K <- matrix(nrow=n,ncol=n)
      K[lower.tri(K,diag=TRUE)] <- makeKMat0(newdata, kernel=object$kernel, kpar=object$kpar)
      K[upper.tri(K,diag=FALSE)] <- t(K)[upper.tri(K,diag=FALSE)]
    }
    if (class(object$svms)[1]!="kernelksvms") classpred <- predict(object$svms,newdata=NULL,KMat=K)
    else classpred <- predict(object$svms,newdata=NULL,KMat=K,trndt=object$x)
  } else {
    if (!is.matrix(newdata))  stop("Wrong value for argument newdata\n")
#    if (ncol(newdata)!=n) stop("Wrong dimensions for the K matrix given by argument newdata\n")
    if (class(object$svms)[1]!="kernelksvms") classpred <- predict(object$svms,newdata=NULL,KMat=newdata)
    else classpred <- predict(object$svms,newdata=NULL,KMat=newdata,trndt=object$x)
  }

  for (obs in 1:n) {
    if (k>2) alpha <- LPdual(factor(as.matrix(classpred)[obs,]),object$grid,eta,probepsilon)
    else {
#      if (!withLP) alpha <- k2Prob(factor(as.matrix(classpred)[obs,]),object$grid[,1])
#      else alpha <- LPdual(factor(as.matrix(classpred)[obs,]),object$grid,eta)
      alpha <- LPdual(factor(as.matrix(classpred)[obs,]),object$grid,eta,probepsilon)
    }
    if (is.null(alpha)) warning("Class probability estimation failed for observation",obs,"\n")
    else Pprob[obs,] <- alpha
  }
  if (!is.null(rownames(newdata))) rownames(Pprob) <- rownames(newdata)
  else rownames(Pprob) <- 1:nrow(newdata)
  colnames(Pprob) <- object$svms$grplvls

  class(Pprob) <- "ClassProb"

  if (retallprd) {
    return(list(Pprob=Pprob,classpred=classpred))
  }

  Pprob
}

k2Prob <- function(Yhat,frstcgrid,epsilon=frstcgrid[1]/2)
{
   Yhat <- unclass(Yhat)
   alphalob <- which.min(Yhat==2)
   alphaupb <- which.max(Yhat==1)
   p1 <- max((frstcgrid[alphalob]+frstcgrid[alphaupb])/2,epsilon)
   c(1-p1,p1)
}

