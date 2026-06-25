SetUpOptPar <- function(start=c("allub","alllb","compubwithlb","warmstarts"), alpha0=NULL,
                        epsilon=1e-12, maxiter=100000, tol=1.5e-12 )
{
#  bselcrit <- match.arg(bselcrit)
  start <- match.arg(start)
  if (start!="warmstarts" && !is.null(alpha0) && length(alpha0)>0)
  {
    alpha0 <- NULL
    warning("Value of the argument alpha0 is ignored when argument start is not set to 'warmstarts'\n")
  }
  else if ( start=="warmstarts" && is.null(alpha0) )
    stop("The argument alpha0 needs to be set to a vector of initial values when argument start is set to 'warmstarts'\n")
  if ( is.null(alpha0) )  alpha0 <- numeric()
  if (!is.numeric(alpha0)) stop("Wrong value for the argument alpha0\n")

  list(start=start,alpha0=alpha0,epsilon=epsilon,maxiter=maxiter,tol=tol)
}

CompOptPar <- function(parlist)
{
  if (is.null(parlist$start)) parlist$start <- "allub"
  if (parlist$start!="warmstarts" && !is.null(parlist$alpha0) && length(parlist$alpha0)>0)
  {
    parlist$alpha0 <- NULL
    warning("Value of the OptCntrl list component alpha0 is ignored when component start is not set to 'warmstarts'\n")
  }
  else if ( parlist$start=="warmstarts" && is.null(parlist$alpha0) )
    stop("The OptCntrl list component alpha0 needs to be set to a vector of initial values when component start is set to 'warmstarts'\n")
  if ( is.null(parlist$alpha0) )  parlist$alpha0 <- numeric()
  if (!is.numeric(parlist$alpha0)) stop("Wrong value for the OptCntrl list component alpha0\n")

  if (is.null(parlist$start)) parlist$start <- "compubwithlb"
  if (is.null(parlist$alpha0)) parlist$alpha0 <- numeric()
  if (is.null(parlist$epsilon)) parlist$epsilon <- 1e-12
  if (is.null(parlist$maxiter)) parlist$maxiter <- 100000
  if (is.null(parlist$tol)) parlist$tol <- 1.5e-12

  parlist
}
