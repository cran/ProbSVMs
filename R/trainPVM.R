require(kernlab)

trainPVM <- function(x=NULL, y, scaled=TRUE, K=NULL, loss=c("WW","LLW"), withbias=FALSE,
                    kernel=c("rbfdot","vanilladot","polydot"), kpar=list(sigma="d2median"),
                    C=0.25, lambda=NULL, tunex=x, tuney=y, tuneK=K,
                    grid=NULL, dpiinv=ceiling(sqrt(length(y))/0.2), keepdt=TRUE, ... )
{
  maxn <- 10000    # Review this limite later, particularly when cached implementations become available

  loss <- match.arg(loss)
  if (is.null(K)) kernel <- match.arg(kernel)
  else kernel <- NULL
  if (!is.factor(y)) stop("Argument y is not a factor\n")
  n <- length(y)
  if (n>maxn)  stop(paste("Currently trainPVM can only handle problems with up to",maxn,"observations (examples)\n"))
  k <- length(levels(y))
  if (!is.null(x) && nrow(x)!=n) stop("Length of argument y (",n,") does not agree wit the number of rows (",nrow(x),") in argument x\n")

  if (!is.null(lambda)) {
    if (length(lambda)!=1 ) stop("Wrong value for the lambda argument\n")
      if (lambda=="tuneit") C <- "tuneit"
      else {
        if (!is.numeric(lambda)) stop("Wrong value for the lambda argument\n")
        C <- 1/(2*n*lambda)
      }
  } else {
    if (length(C)!=1 ) stop("Wrong value for the C argument\n")
    if (C=="tuneit") lambda <- "tuneit"
    else {
      if (!is.numeric(C)) stop("Wrong value for the C argument\n")
      lambda <- 1/(2*n*C)
    }
  }

  if (is.null(K)) {
    if (!is.matrix(x) &&  !is.data.frame(x)) stop("Argument x is not a matrix neither a data frame\n")
    if (all(apply(x,2,is.numeric))) {
      if (scaled) x <- scale(x,center=FALSE)
    } else {
      if (!is.data.frame(x)) stop("Argument x inludes categorical variables but is not a data frame\n")
      if (scaled) {
        numvar <- NULL                        #Note: apply(x,2,is.numeric)) does not work properly if some x variables are non-numeric
        for (i in 1:ncol(x)) if (is.numeric(x[,i])) numvar <- c(numvar,i)
        if (length(numvar)>0) {
          numxscld <- scale(x[,numvar])
          x[,numvar] <- numxscld
        }
      }
      x <- model.matrix(~ .,data=x)[,-1]  # Converting factors into dummy variables
      if (scaled) {
        p <- ncol(x)
        attr(x,"scaled:scale") <- rep(1.,p)
        if (length(numvar)>0) {
          attr(x,"scaled:scale")[numvar]  <- attr(numxscld,"scaled:scale")
        }
      }
    }
    if (nrow(x)!=n ) stop("Number of rows of argument x do not agree with length of argument y\n")
    if (!is.null(kpar) && !is.list(kpar))  stop("Wrong type for argument kpar")
    if (kernel=="Gauss") {
      if (is.null(kpar)) stop("The argument kpar must alwyas be specified for Gaussian kernels")
      if (is.null(kpar$sigma)) stop("The list provided by the kpar argument does not have a component named sigma as required\n")
      if (kpar$sigma=="d2median" || kpar$sigmar=="d2q01q09mean" || kpar$sigma=="d2q01q09hmean")   #Note: hmean stands for harmonic mean
        kpar <- GetrbfdotSigPar(x,kpar$sigma)
      if (!is.numeric(kpar$sigma) || length(kpar$sigma)!=1) stop("Wrong value for the sigma component of list provided by the kpar argument\n")
    }
    K <- makeKMat0(x, kernel=kernel, kpar=kpar) # Note: For efficiency the elements of K will be passed as a vector and not a matrix !!
  } else {
    if (!is.matrix(K)) {
      if (!is.numeric(K) || length(K)!=n*(n+1)/2) stop("Wrong value for argument K\n")
    } else {
      if (nrow(K)!=n || ncol(K)!=n)
        stop("Dimensions of the matrix given by argument K do not agree with length of argument y\n")
      if (!isTRUE(all.equal(K,t(K))))
        stop("Matrix given by argument K is not symmetric as it should\n")
      K <- K[lower.tri(K,diag=TRUE)]
    }
  }

  if ( !is.numeric(lambda) ) {

    if ( length(lambda)!=1 || lambda!="tuneit" )
      stop("Argument lambda is not a real number (of class numeric), neither the string 'tuneit' as it should\n")
    if (!is.factor(tuney)) stop("Argument tuney is not a factor as it should\n")
    if ( is.null(tunex) && is.null(tuneK)) stop("Arguments tunex and tuneK cannot be simultaneously NULL\n")
    n1 <- length(tuney)

    if ( !is.null(tuneK) && !is.matrix(tuneK) ) {
      tuneK <- matrix(0.,nrow=n,ncol=n)
      tuneK[lower.tri(tuneK,diag=TRUE)] <- K
      tuneK[upper.tri(tuneK,diag=FALSE)] <- tuneK[lower.tri(tuneK,diag=FALSE)]
    } else  if (is.null(tuneK)) {
      if (!is.matrix(tunex) && !is.data.frame(tunex)) stop("Argument tunex is not a matrix neither a data frame\n")
      if (all(apply(tunex,2,is.numeric))) {
        if (scaled) tunex <- scale(tunex,scale=attr(x,"scaled:scale"))
      } else {
        if (!is.data.frame(tunex)) stop("Argument tunex inludes categorical variables but is not a data frame\n")
        if (scaled) {
          numvar <- NULL            # Note: apply(tunex,2,is.numeric)) does not work properly if some tunex variables are non-numeric
          for (i in 1:ncol(tunex)) if (is.numeric(tunex[,i])) numvar <- c(numvar,i)
          if (length(numvar)>0) {
            numxscld <- scale(tunex[,numvar],scale=attr(x[,numvar],"scaled:scale"))
            tunex[,numvar] <- numxscld
          }
        }
        tunex <- model.matrix(~ .,data=tunex)[,-1]  # Converting factors into dummy variables
      }

      if (nrow(tunex)!=n1 ) stop("Number of rows of argument tunex do not agree with length of argument tuney\n")
      tuneK <- makeKMat0(x, tunex, kernel=kernel, kpar=kpar)
    } else {
      if (!is.matrix(tuneK)) stop("Argument tuneK is not a matrix as it should\n")
      if (nrow(tuneK)!=n || ncol(tuneK)!=n1)
        stop("Dimensions of the matrix given by argument tuneK do not agree with length of argument tuney\n")
    }

  }
  if (is.null(grid)) grid <- buildgrid(k,dpiinv)

  if (lambda=="tuneit") {
    bestC <- Csearch(zereoexp="center",algorithm=trainandevalPVM,K=K,grp=y,testK=tuneK,testgrp=tuney,grid=grid,loss=loss)
    lambda <- 1/(2*n*bestC)
  }

#  PVM <- list(nclasses=k,kernel=kernel,kpar=kpar,grid=grid)
  PVM <- list(nclasses=k,kernel=kernel,kpar=kpar,grid=grid,scale=attr(x,"scaled:scale"))
  if (!withbias) PVM$svms <- trainSVM(y=y, K=K, lambda=lambda, loss=loss, class.weights=grid, scaled=scaled, kernel=NULL, kpar=NULL,  keepdt=FALSE, ...)
  else {
    if (k!=2) stop("Currrently PVM with bias coefficients only works for classification problems with 2 classes\n")
    Kmat <- matrix(nrow=n,ncol=n)
    Kmat[lower.tri(Kmat,diag=TRUE)] <- K
    Kmat[upper.tri(Kmat,diag=FALSE)] <- t(Kmat[lower.tri(Kmat,diag=FALSE)])
    Kmat <- as.kernelMatrix(Kmat)
    cllvls <- levels(y)
    colnames(grid) <- cllvls
    npoints <- nrow(grid)
    PVM$svms <- list(alpha=NULL,bias=NULL,ksvms=vector("list",npoints),grplvls=cllvls,lambda=lambda,C=C,
                     x=NULL,scale=NULL,kernel=NULL,kpar=NULL)
    for (pnt in 1:npoints) PVM$svms$ksvms[[pnt]] <- ksvm(Kmat, y=y, C=C, class.weights=grid[pnt,], ...)
    class(PVM$svms) <-c("kernelksvms","kernelSVMs","kernelSVM")
    keepdt <- TRUE
  }
  if (keepdt) PVM$x <- x
  class(PVM) <- "PVM"
  PVM
}

trainandevalPVM <- function(K, grp, testK, testgrp, C, grid, loss, evalcrt=c("loglik","errest"), minloglik=-1e12, ...)
{
  evalcrt <- match.arg(evalcrt)
  npoints <- nrow(grid)
  k <- ncol(grid)
  N <- nrow(testK)
  n <- length(grp)

  trnres <- trainPVM(K=K, y=grp, loss=loss, scaled=FALSE, lambda=1/(2*n*C), grid=grid, ...) #Note: Try one version with K, y and grid on the parent environment for efficiency!!
  testPprob <- predict(trnres,newdata=testK,newdataasKmat=TRUE)
  evalres <- evalPP(testPprob,testgrp)
  if (evalcrt=="loglik") {
    if (evalres$loglik > minloglik) return(evalres$loglik)
    else return(minloglik-evalres$errest)
  }
  else if (evalcrt=="errest") return(-evalres$errest)
}
