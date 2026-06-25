l2norm2 <- function(x) sum(x^2)
GausKern <- function(x,y,sigma) exp(-l2norm2(x-y)/sigma^2)

MedianDist1Set <- function(x)
{
  n <- nrow(x)
  dist <- numeric(n*(n+1)/2)
  ind <- 0
  for (i in 1:n) for (j in 1:i) {
    ind <- ind+1
    dist[ind] <- l2norm2(x[i,]-x[j,])
  }
  sqrt(median(dist))
}

MedianDist2Sets <- function(xpos,xneg)
{
  npos <- nrow(xpos)
  nneg <- nrow(xneg)
  dist <- numeric(npos*nneg)
  ind <- 0
  for (i in 1:npos) for (j in 1:nneg) {
    ind <- ind+1
    dist[ind] <- l2norm2(xpos[i,]-xneg[j,])
  }
  sqrt(median(dist))
}

GetrbfdotSigPar <- function(x,kpar=c("d2median","d2q01q09mean","d2q01q09hmean"))
{
  kpar <- match.arg(kpar)
  d2srange <- sigest1(x,scaled=FALSE)    # Note: sigest1, unlike kernlab sigest, by default uses the whole training data, and not a sample of it
  if (kpar=="d2median") kpar <- list(sigma=1./d2srange[2])
  else if (kpar=="d2q01q09mean") kpar <- list(sigma=2./(d2srange[1]+d2srange[3]))
  else if (kpar=="d2q01q09hmean") kpar <-list(sigma=(1./d2srange[1]+1./d2srange[3])/2)

  kpar
}

