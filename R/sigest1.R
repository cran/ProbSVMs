sigest1 <- function (x, frac=0.5, scaled=TRUE, na.action=na.omit, sample=FALSE, smpwithrep=TRUE)
{
  if (!is.matrix(x) || !all(apply(x,2,is.numeric)))
    stop("Argument x is not a matrix of real numbers as it should\n")
  if(!sample && !smpwithrep)
    warning("Since sample argument is set to its FALSE default value (no sampling) the smpwithrep will be ignored\n")

  x <- na.action(x)

  if (length(scaled) == 1)
    scaled <- rep(scaled, ncol(x))
    if (any(scaled)) {
      co <- !apply(x[,scaled, drop = FALSE], 2, var)
      if (any(co)) {
        scaled <- rep(FALSE, ncol(x))
          warning(paste("Variable(s)",
            paste("`",colnames(x[,scaled, drop = FALSE])[co],"'", sep="", collapse=" and "),"constant. Cannot scale data.")
                 )
        } else {
          xtmp <- scale(x[,scaled])
          x[,scaled] <- xtmp
      }
    }

    m <- nrow(x)
    if (sample) {
      n <- floor(frac*m)
      if (smpwithrep) {
        index <- sample(1:m, n, replace = TRUE)
        index2 <- sample(1:m, n, replace = TRUE)
      } else {
        index <- sample(1:m, n, replace = FALSE)
        index2 <- sample(1:m, n, replace = FALSE)
      }
      temp <- x[index,, drop=FALSE] - x[index2,,drop=FALSE]
      edist <- rowSums(temp^2)
    } else edist <- as.vector(dist(x))^2

    res <- 1. / quantile(edist[edist!=0],probs=c(0.9,0.5,0.1))
    attr(res,"names") <- NULL

    res
}
