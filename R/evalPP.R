evalPP <- function(ePP,trcl)
{
  if (!is.matrix(ePP)) stop("Argument epp is not a matrix\n")
  if (!is.factor(trcl)) stop("Argument trcl is not a factor\n")
  N <- nrow(ePP)
  k <- ncol(ePP)
  if (length(trcl)!=N) 
    stop("Dimensions of argument ePP (",dim(ePP),") does not agree with length argument trcl (",length(trcl),")\n")

  grplvls <- levels(trcl)
  prdclass <- grplvls[apply(ePP,1,which.max)]
  errcnt <- length(which(prdclass!=trcl))
  loglik <- 0.
  for (g in 1:k) {
    obsing <- which(trcl==grplvls[g])
    loglik <- loglik + sum(log(ePP[obsing,g]))
  }  
  list(loglik=loglik,errest=100*errcnt/N)
}
