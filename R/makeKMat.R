makeKMat <- function(dt1, dt2=NULL, kernel=c("rbfdot","vanilladot","polydot"), kpar=list(sigma = "d2median"))
{
  K0 <- makeKMat0(dt1,dt2,kernel,kpar)
  if (!is.null(dt2)) return(K0)
  n <- nrow(dt1)
  K <- matrix(nrow=n,ncol=n)
  K[lower.tri(K,diag=TRUE)] <- K0
  K[upper.tri(K,diag=FALSE)] <- t(K)[upper.tri(K,diag=FALSE)]

  K
}

makeKMat0 <- function(dt1, dt2=NULL, kernel=c("rbfdot","vanilladot","polydot"), kpar=list(sigma = "d2median"))
{

  kernel <- match.arg(kernel)
  if (!is.matrix(dt1) && !is.data.frame(dt1)) stop("Wrong type for argument dt1")
  if (!is.null(dt2) && !is.matrix(dt2) && !is.data.frame(dt2)) stop("Wrong type for argument dt2")
  if (!is.matrix(dt1))  dt1 <- as.matrix(dt1)
  if (!is.null(dt2) && !is.matrix(dt2))  dt2 <- as.matrix(dt2)

  if  (kernel=="rbfdot") {
    if (is.null(kpar)) stop("The argument kpar must alwyas be specified for Gaussian kernels")
    if (!is.list(kpar))  stop("Wrong type for argument kpar")
    if (is.null(kpar$sigma)) stop("The list provided by the kpar argument does not have a component named sigma as required\n")
    if (kpar$sigma=="d2median" || kpar$sigmar=="d2q01q09mean" || kpar$sigma=="d2q01q09hmean")   #Note: hmean stands for harmonic mean
 #     if (!is.null(dt2)) stop("Wrong type for argument kpar\n")
      kpar <- GetrbfdotSigPar(dt1,kpar$sigma)
    if (!is.numeric(kpar$sigma) || length(kpar$sigma)!=1) stop("Wrong value for the sigma component of list provided by the kpar argument\n")
    if (is.null(dt2)) .Call("makeGauKMat", dt1,  kpar$sigma, PACKAGE = "ProbSVMs" )
    else .Call("makeGauKMat2", dt1, dt2,  kpar$sigma, PACKAGE = "ProbSVMs" )
  } else if  (kernel=="vanilladot") {
    warning("As the linear kernel is not stricly positive definite, the generalization properties of the resulting non-bias kernel SVM cannot be guaranteed\n")
    if (is.null(dt2)) .Call("makeLinKMat", dt1,  PACKAGE = "ProbSVMs" )
    else .Call("makeLinKMat2", dt1, dt2,  PACKAGE = "ProbSVMs" )
  } else if  (kernel=="polydot") {
    warning("As the polynomial kernel is not stricly positive definite, the generalization properties of the resulting non-bias kernel SVM cannot be guaranteed\n")
    if (is.null(kpar)) kpar <- list(degree=1,scale=1.,offset=1.)
    if (!is.list(kpar))  stop("Wrong type for argument kpar")
    if (is.null(kpar$degree))  kpar$degree <- 1
    if (is.null(kpar$scale))  kpar$scale <- 1.
    if (is.null(kpar$degree))  kpar$offset <- 1.
    if (is.null(dt2)) .Call("makePolKMat", dt1, kpar$degree, kpar$scale, kpar$offset,  PACKAGE = "ProbSVMs" )
    else .Call("makePolKMat", dt1, dt2, kpar$degree, kpar$scale, kpar$offset,   PACKAGE = "ProbSVMs" )
  }

}
