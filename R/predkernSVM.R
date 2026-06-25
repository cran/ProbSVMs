predict.kernelSVM <- function(object, ..., newdata=NULL, KMat=NULL, trndt=NULL)
{
  if (is.null(KMat)) {
    if (is.null(object$x) && is.null(trndt))
      stop("No training data was provided (neither by the classifier given by the argument object nor directly by the argument trndt\n")
    if (!is.null(object$x)) ntrndtcol <- ncol(object$x)
    else ntrndtcol <- ncol(trndt)
    if (is.data.frame(newdata) || is.matrix(newdata)) {
      for (v in 1:ncol(newdata)) if (!is.numeric(newdata[,v])) stop("Not all variables of newdata are numeric as they should\n")
      if (ntrndtcol!=ncol(newdata))
        stop("Dimension of the newdata vectors (=",ncol(newdata),") does not agree with the dimension of the training data vectors (=",ntrndtcol,")\n")
      if (!is.null(object$center)) {
        if (!is.null(object$scale)) {
          newdata <- scale(newdata,center=object$center,scale=object$scale)
        } else {
          newdata <- scale(newdata,center=object$center,scale=FALSE)
        }
      } else {
        if (!is.null(object$scale)) {
          newdata <- scale(newdata,center=FALSE,scale=object$scale)
        }
      }
      if (is.data.frame(newdata))  newdata <- as.matrix(newdata)
      if (!is.null(object$x)) KMat <- makeKMat0(object$x, newdata, kernel=object$kernel, kpar=object$kpar)
      else KMat <- makeKMat0(trndt, newdata, kernel=object$kernel, kpar=object$kpar)


    } else if (is.numeric(newdata) || is.integer(newdata)) {
      if (ntrndtcol!=length(newdata))
        stop("Dimension of the newdata vectors (=",length(newdata),") does not agree with the dimension of the training data vectors (=",ntrndtcol,")\n")
      if (!is.null(object$center)) {
        if (!is.null(object$scale)) {
          newdata <- (newdata-object$center)/object$scale
        } else {
          newdata <- newdata-object$center
        }
      } else {
        if (!is.null(object$scale)) {
          newdata <- newdata/object$scale
        }
      }
      if (!is.null(object$x)) KMat <- makeKMat0(object$x,  t(as.matrix(newdata)), kernel=object$kernel, kpar=object$kpar)
      else KMat <- makeKMat0(trndt,  t(as.matrix(newdata)), kernel=object$kernel, kpar=object$kpar)

    } else {
      stop("Argument newdata is not a matrix, a data.frame nor an numeric or integer vector\n")
    }

  }

  if (class(object)[1]=="kernelSVM") {
    factor(object$grplvls[apply(KMat%*%object$alpha,1,which.max)],levels=object$grplvls)
  } else if (class(object)[1]=="kernelSVMwbias") {
      factor(object$grplvls[apply(KMat%*%object$alpha+object$bias,1,which.max)],levels=object$grplvls)
  } else if (class(object)[1]=="kernelksvm") {
    predict(object$ksvms,as.kernelMatrix(KMat[,SVindex(object$ksvms),drop=FALSE]))

#  } else if (class(object)[1]=="kernelSVMs") {
  } else if (class(object)[1]=="kernelSVMs" || class(object)[1]=="kernelSVMswbias" || class(object)[1]=="kernelksvms" ) {
    if (!is.null(object$alpha)) nsvms <- dim(object$alpha)[3]
    if (class(object)[1]=="kernelSVMs") {
      res <-  data.frame(factor(object$grplvls[apply(KMat%*%object$alpha[,,1],1,which.max)],levels=object$grplvls))
      for (svm in 2:nsvms)  res <- cbind(res,factor(object$grplvls[apply(KMat%*%object$alpha[,,svm],1,which.max)],levels=object$grplvls))
    } else if (class(object)[1]=="kernelSVMswbias") {
        res <-  data.frame(factor(object$grplvls[apply(KMat%*%object$alpha[,,1]+object$bias[1,],1,which.max)],levels=object$grplvls))
        for (svm in 2:nsvms)  res <- cbind(res,factor(object$grplvls[apply(KMat%*%object$alpha[,,svm]+object$bias[svm,],1,which.max)],levels=object$grplvls))
    } else if (class(object)[1]=="kernelksvms") {
      nsvms <- length(object$ksvms)
      res <-  predict(object$ksvms[[1]],as.kernelMatrix(KMat[,SVindex(object$ksvms[[1]]),drop=FALSE]))
      for (svm in 2:nsvms)  res <- cbind(res,predict(object$ksvms[[svm]],as.kernelMatrix(KMat[,SVindex(object$ksvms[[svm]]),drop=FALSE])))
    }

    names(res) <- paste("SVM",1:nsvms,sep="")
    res
  } else stop("Wrong class of object argument\n")

}
