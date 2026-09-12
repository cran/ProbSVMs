require(lpSolveAPI)

LPdual <- function(Yhat,grid,eta=15.,epsilon=min(grid)/2)
{
  if ( !is.numeric(epsilon) || length(epsilon) > 1 || epsilon[1] < 0 || epsilon[1] > 1)
    stop("Wrong value for argument epsilon\n")
  if (!is.matrix(grid)) {
    if  (!is.data.frame(grid)) stop("Argument grid is not a matrix neither a data frame\n")
    grid <- as.matrix(grid)
  }
#  if (!is.factor(Yhat)) stop("Argument Yhat is not a factor as it should be\n")
  Yhat <- unclass(Yhat)

  k <- ncol(grid)
  npoints <- nrow(grid)
  npconst <- (k-1)*npoints;
  objcoef <- c(1.,-1,rep(epsilon,k),rep(0,npconst))
  M <- matrix(0.,nrow=k,ncol=2+k+npconst)
  dir <- rep("=",k)
  rhs <- numeric(k)
  ubs <- c(rep(Inf,2+k),rep(eta-1,npconst))

  M[1:k,1] <- 1.
  M[1:k,2] <- -1.
  for (c in 1:k) {
    M[c,c+2] <- 1.
    colib <- k+2
    sumPiy <- sumPic <- 0.
    for (pnt in 1:npoints) {
      gridval <- grid[pnt,c]
      Y <- Yhat[pnt]
      if (c==Y) {
        M[c,colib+1:(k-1)] <- gridval
        sumPiy <- sumPiy + (k-1) * gridval
      }
      else {
        if (c<Y) M[c,colib+c] <- -gridval
        else M[c,colib+c-1] <- -gridval
        sumPic <- sumPic + gridval
      }
      colib <- colib + k-1
    }
    rhs[c] <- sumPic - sumPiy
  }

  model <- make.lp(ncol=ncol(M))
  set.objfn(model,-objcoef)           # Note: lpSolveAPI always minimizes, therefore the sign change for the obj coefficients !!!
  for (c in 1:nrow(M)) add.constraint(model,M[c,],dir[c],rhs[c])
  set.bounds(model,upper=ubs)
#  nulllarg <- 0
#  status <- solve.lpExtPtr(model,nulllarg)
#  status <- solve.lpExtPtr(model,0)
  status <- solve.lpExtPtr(model,NULL)
  if (status!=0) return(NULL)

  -get.dual.solution(model)[2:(k+1)]   # First dual variable is the second element returned by lpSolveAPI,
                                        # and the minus sign is because we minimized the symmetric of the (max) objective function!!
}
