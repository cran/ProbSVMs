GenCrossVFolds <- function(grouping, Strfolds=TRUE, kfold=10, CVrep=3, retind=c("aslists","asmatrices"))
{
   fold <- function(n,kfold,fi) {
    lb <- floor((fi-1)*n/kfold) + 1
    ub <- floor(fi*n/kfold)
    if (lb>ub) return(NULL)
    lb:ub  # return(lb:ub)
  }
  
  retind <- match.arg(retind)
  if (retind=="asmatrices") trainncol <- tunncol <- 0

  if (!is.factor(grouping)) {
    grouping <- factor(grouping)
  }
  codes <- levels(grouping)
  nk <- table(grouping)
  k <- nrow(nk)
  if (k<2) stop("Factor grouping should have at least two different levels\n")
  nk <- as.numeric(nk)
  n <- sum(nk)
  if (Strfolds) permut <- vector("list",k)
  trep <- kfold*CVrep
  trainfolds <- vector("list",trep)
  tunfolds <- vector("list",trep)
  for (i in 1:CVrep)  {
    if (Strfolds) for (grp in 1:k) permut[[grp]] <- sort.int(runif(nk[grp]),index.return=TRUE)$ix
    else permut <- sort.int(runif(n),index.return=TRUE)$ix
    for (j in 1:kfold) {
      rep <- (i-1)*kfold + j
      if (Strfolds) {
        out <- which(grouping==codes[1])[permut[[1]]][fold(nk[1],kfold,j)]
        for (grp in 2:k) out <- c(out,which(grouping==codes[grp])[permut[[grp]]][fold(nk[grp],kfold,j)])
      } else {
        out <- permut[fold(n,kfold,j)]
      }
      if (length(out)==0) next
      trainfolds[[rep]] <- setdiff(1:n,out)
      tunfolds[[rep]] <- out
      if (retind=="asmatrices") {
        nout <-  length(out)
        tunncol <- max(tunncol,nout)
        trainncol <- max(trainncol,n-nout)
      }  
    }
  }
  if (retind=="asmatrices") {
    trfoldsasmat <- matrix(0,nrow=trep,ncol=trainncol)
    for (rep in 1:trep) trfoldsasmat[rep,1:length(trainfolds[[rep]])] <- trainfolds[[rep]]  
    tufoldsasmat <- matrix(0,nrow=trep,ncol=tunncol)
    for (rep in 1:trep) tufoldsasmat[rep,1:length(tunfolds[[rep]])] <- tunfolds[[rep]]  
    return(list(trIndL=trfoldsasmat,tunIndL=tufoldsasmat))
  }  
  list(trIndL=trainfolds,tunIndL=tunfolds)
}
