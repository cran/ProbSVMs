trainSVM <- function(x=NULL, y, class.weights=rep(1.,length(levels(y))), scaled=TRUE, K=NULL, loss=c("WW","LLW"),
                     kernel=c("rbfdot","vanilladot","polydot"), kpar=list(sigma="d2median"),
                     C=1., lambda=NULL, keepdt=TRUE, retotpst=FALSE, OptCntrl=SetUpOptPar(), TunCntrl=SetUpTunPar() )
{
  valid.weights <- function(w) isTRUE(all.equal(w,rep(1.,length(w)))) || isTRUE(all.equal(sum(w),1.))

  maxn <- 10000    # Review this limite later, particularly when cached implementations becoma available

  if (!is.factor(y)) stop("Argument y is not a factor as it should\n")
  ylvls <- levels(y)
  n <- length(y)
  if (n>maxn)  stop(paste("Currently trainSVM can only handle problems with up to",maxn,"observations (examples)\n"))
  k <- length(ylvls)
  OptCntrl <- CompOptPar(OptCntrl)
  TunCntrl <- CompTunPar(TunCntrl)
  if (is.null(TunCntrl$tuneK)) { }
  loss <- match.arg(loss)
  if (is.null(K)) kernel <- match.arg(kernel)
  else kernel <- NULL

  if (missing(C) && !is.null(lambda)) {
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

  if ( !is.numeric(class.weights) && !is.matrix(class.weights) ) stop("Wrong value for argument class.weights\n")
  if (!is.matrix(class.weights)) {
    if (length(class.weights)!=k) stop("Incorrect number of elements in the vector given by argument class.weights\n")
    if ( !valid.weights(class.weights))
      stop("Sum of the elements in the vector given by argument class.weights is not equal to one as it should\n")
  }  else  {
    if (ncol(class.weights)!=k) stop("Incorrect number of columns in the matrix given by argument class.weights\n")
    if ( !all(apply(class.weights,1,valid.weights)))
      stop("At least one row of the matrix given by argument class.weights is not a valid set of weight specifications\n")
  }
  if ( is.null(x) && is.null(K)) stop("Arguments x and K cannot be simultaneously NULL\n")
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
      K <- K[lower.tri(K,diag=TRUE)]  # Note: For efficiency the elements of K will be passed as a vector and not a matrix !!
    }
  }

  TrTunInd <- list()
  if ( !is.numeric(lambda) ) {

    if ( length(lambda)!=1 || lambda!="tuneit" )
      stop("Argument lambda is not a real number (of class numeric), neither the string 'tuneit' as it should\n")

    if (TunCntrl$crossval) {
      TunCntrl$tunex <- TunCntrl$tuneK <- matrix(nrow=0,ncol=0)
      if (!is.null(TunCntrl$tuney)) warning("Value of argument TunCntrl$tuney will be ignored, since training data will be used for tuning\n")
      TunCntrl$tuney <- factor(integer())
      n1 <- 0
      TrTunInd <- GenCrossVFolds(y,TunCntrl$crossvalpar$Strfolds,TunCntrl$crossvalpar$kfold,TunCntrl$crossvalpar$CVrep)
    } else {
      if ( is.null(TunCntrl$tunex) && is.null(TunCntrl$tuneK)) {
        TunCntrl$tunex <- x
        TunCntrl$tuney <- y
        if (is.null(K)) TunCntrl$tuneK <- makeKMat0(x, TunCntrl$tunex, kernel=kernel, kpar=kpar)
        else {
          TunCntrl$tuneK <- matrix(nrow=n,ncol=n)
          TunCntrl$tuneK[lower.tri(TunCntrl$tuneK,diag=TRUE)] <- K
          TunCntrl$tuneK[upper.tri(TunCntrl$tuneK,diag=FALSE)] <- t(TunCntrl$tuneK)[upper.tri(TunCntrl$tuneK,diag=FALSE)]
        }
      } else {
        if (!is.factor(TunCntrl$tuney)) stop("Argument TunCntrl$tuney is not a factor as it should\n")
        n1 <- length(TunCntrl$tuney)
        if (is.null(TunCntrl$tuneK)) {
          if (!is.matrix(TunCntrl$tunex) && !is.data.frame(TunCntrl$tunex)) stop("Argument TunCntrl$tunex is not a matrix neither a data frame\n")
          if (all(apply(TunCntrl$tunex,2,is.numeric))) {
            if (scaled) TunCntrl$tunex <- scale(TunCntrl$tunex,scale=attr(x,"scaled:scale"))
          } else {
            if (!is.data.frame(TunCntrl$tunex)) stop("Argument TunCntrl$tunex inludes categorical variables but is not a data frame\n")
            if (scaled) {
              numvar <- NULL            # Note: apply(TunCntrl$tunex,2,is.numeric)) does not work properly if some TunCntrl$tunex variables are non-numeric
              for (i in 1:ncol(TunCntrl$tunex)) if (is.numeric(TunCntrl$tunex[,i])) numvar <- c(numvar,i)
              if (length(numvar)>0) {
                numxscld <- scale(TunCntrl$tunex[,numvar],scale=attr(x[,numvar],"scaled:scale"))
                TunCntrl$tunex[,numvar] <- numxscld
              }
            }
            TunCntrl$tunex <- model.matrix(~ .,data=TunCntrl$tunex)[,-1]  # Converting factors into dummy variables
          }
          if (nrow(TunCntrl$tunex)!=n1 ) stop("Number of rows of argument TunCntrl$tunex do not agree with length of argument TunCntrl$tuney\n")
          TunCntrl$tuneK <- makeKMat0(x, TunCntrl$tunex, kernel=kernel, kpar=kpar)
        } else {
          if (!is.matrix(TunCntrl$tuneK)) stop("Argument TunCntrl$tuneK is not a matrix as it should\n")
          if (nrow(TunCntrl$tuneK)!=n || ncol(TunCntrl$tuneK)!=n1)
            stop("Dimensions of the matrix given by argument TunCntrl$tuneK do not agree with length of argument TunCntrl$tuney\n")
        }
      }
    }
  }

  if (loss=="LLW")  model <- 0
  else if (loss=="WW")  model <- 1
       else stop("Optimization of SVMs using other than the LLW and WW losses were not implemented yet\n")

  if (OptCntrl$start=="compubwithlb") OptCntrlstart <- 0
  else  if (OptCntrl$start=="allub") OptCntrlstart <- 1
        else  if (OptCntrl$start=="alllb") OptCntrlstart <- 2
              else  if (OptCntrl$start=="warmstarts") OptCntrlstart <- 3
                    else stop("Wrong value for component start of the OptCntrl list\n")

  if (is.null(x)) scale <- NULL 
  else scale <- attr(x,"scaled:scale")
  if (!keepdt) x <- NULL
  
  if (!is.matrix(class.weights)) class.weights <- matrix(class.weights,nrow=1)
  nsvms <- nrow(class.weights)
  if (is.numeric(lambda)) {

    tmpres <- .Call("TrSVMs", unclass(y)-1, lambda, K, class.weights,
                    OptCntrl$epsilon, OptCntrl$maxiter, OptCntrl$tol,
                    model, retotpst, OptCntrlstart, OptCntrl$alpha0, PACKAGE = "ProbSVMs")
    if (is.null(tmpres)) return(NULL)

    if (nsvms==1) {

      alphaMat <- alphacnv(tmpres$alpha,y,loss)
      res <- list( alpha=alphaMat, grplvls=ylvls, lambda=lambda, C=C,
#                    x=x, scale=attr(x,"scaled:scale"), kernel=kernel, kpar=kpar
                    x=x, scale=scale, kernel=kernel, kpar=kpar
                  )
      if (retotpst) res$optlist <- list(optvalue=tmpres$fval,iterations=tmpres$iter,hitrate=NULL)
      class(res) <-"kernelSVM"

    } else {
      alphaArray <- array(dim=c(n,k,nsvms))
      for (svm in 1:nsvms)  alphaArray[,,svm] <- alphacnv(tmpres[svm,],y,loss)
      res <- list( alpha=alphaArray, grplvls=ylvls, lambda=lambda, C=C,
#                  x=x, scale=attr(x,"scaled:scale"), kernel=kernel, kpar=kpar
                  x=x, scale=scale, kernel=kernel, kpar=kpar
                  )
      if (retotpst) res$optlist <- list(optvalue=NULL,iterations=NULL,hitrate=NULL)
      class(res) <-c("kernelSVMs","kernelSVM")
    }

  }  else if (lambda=="tuneit") {

    if (nsvms!=1) stop("Currently, only one weight specification can be specifed when argument lambda equals the string 'tuneit'\n")
    tmpres <- NULL

    tmpres <- .Call("TunSVMs", unclass(y)-1, unclass(TunCntrl$tuney)-1, K, TunCntrl$tuneK, class.weights, TrTunInd,
                    OptCntrl$epsilon, OptCntrl$maxiter, OptCntrl$tol,
                    model, Csrchpar=TunCntrl$Csrchpar, retotpst, OptCntrlstart, OptCntrl$alpha0, PACKAGE = "ProbSVMs")

    if (is.null(tmpres)) return(NULL)                  # To do: change this to exit with an warning when called on its own, and exit gracefully otherwise !!!
    alphaMat <- matrix(tmpres$alpha,nrow=n,ncol=k)
    res <- list( alpha=alphaMat, grplvls=ylvls, lambda=1/(2*n*tmpres$bestC), C=tmpres$bestC,
#                 x=x, scale=attr(x,"scaled:scale"), kernel=kernel, kpar=kpar
                 x=x, scale=scale, kernel=kernel, kpar=kpar
                )
    if (retotpst) res$optlist <-  list(optvalue=tmpres$fval,iterations=tmpres$iter,hitrate=tmpres$hitrate)

    class(res) <-"kernelSVM"

  }  else stop("Argument lambda is not a real number (of class numeric), neither the string 'tuneit' as it should\n")

  res

}

alphacnv <- function(svmalpha,grouping,loss)
{
  if (!is.numeric(svmalpha)) stop("Wrong type for argumente svmalpha\n")
  if (!is.factor(grouping)) stop("Wrong type for argument grouping\n")
  m <- length(svmalpha)
  n <- length(grouping)
  k <- length(levels(grouping))
  if (loss=="TwoG") {
    if (k!=2) stop("TwoG losses can only be used in problems with two classes\n")
    return(matrix(c(svmalpha,-svmalpha),ncol=2))
  } else if (k != as.integer(m/n + 1))
    stop("Lengths of svmalpha (",m,"), and grouping (",n,") do not agree with the number of grouping classes (",k,")\n")
  y <- unclass(grouping)
  alphayind <- k*(0:(n-1)) + y
  if (any(alphayind < 1) || any(y > n*k)) stop("Some values of alphayind are outside their allowable range\n")

  talpha <- matrix(0.,nrow=k,ncol=n)
  if (loss=="LLW") talpha[-alphayind] <- -svmalpha
  else if (loss=="WW") {
    talpha[-alphayind] <- -svmalpha/2
    for (obs in 1:n) {
      baseind <- (obs-1)*(k-1)
      talpha[y[obs],obs] <- sum(svmalpha[baseind+(1:(k-1))])/2
    }
  } else stop("Currently only TwoG, LLW and WW loses are implemented\n")

  t(talpha)
}
