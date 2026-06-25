buildgrid <- function(k,dpiinv)
{
  if (!is.wholenumber(dpiinv)) stop("Argument dpiinv needs to be an integer\n")
  dpi <- 1./dpiinv
  np1 <- dpiinv - 1
  if (k>2) npoints <- np1*k
  else npoints <- np1
  grid  <- matrix(nrow=npoints,ncol=k)

  grdig <- dpi
  for (i in 1:np1) {
    if (k==2) {
      grid[i,1] <- grdig
      grid[i,2] <- 1. - grdig
    } else {
      for (g in 1:k) {
        pind <- (g-1)*np1 + i
        grid[pind,g] <- grdig
        u <- runif(k-1)
        grid[pind,-g] = u*(1.-grdig)/sum(u);
      }
    }
    grdig <- grdig+dpi
  }

  grid
}
