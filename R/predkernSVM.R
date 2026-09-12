predict.kernelSVM <- function(object, ..., newdata=NULL, KMat=NULL, trndt=NULL)
{
  if (!is.null(object$x) && !is.null(trndt)) 
    warning("The value of argument trndt will be ignored, since the training data was already suplied by the x component of the value of argument object\n")  

  if (is.null(KMat)) {
    if (is.data.frame(newdata) || is.matrix(newdata)) {
      if (!is.null(object$center)) {
        if (!is.null(object$scale)) newdata <- scale(newdata,center=object$center,scale=object$scale)
        else newdata <- scale(newdata,center=object$center,scale=FALSE)       
      } else if (!is.null(object$scale)) newdata <- scale(newdata,center=FALSE,scale=object$scale)
      if (is.data.frame(newdata))  newdata <- as.matrix(newdata)      
    } else if (is.numeric(newdata) || is.integer(newdata)) {
      if (!is.null(object$center)) {
        if (!is.null(object$scale)) newdata <- (newdata-object$center)/object$scale
        else newdata <- newdata-object$center        
      } else if (!is.null(object$scale)) newdata <- newdata/object$scale
    } else  stop("Argument newdata is not a matrix, a data.frame nor an numeric or integer vector\n")
  
    if (is.character(trndt) && trndt=="newdata") {
      n <- nrow(newdata)
      KMat <- matrix(nrow=n,ncol=n)
      KMat[lower.tri(KMat,diag=TRUE)] <- makeKMat0(newdata, kernel=object$kernel, kpar=object$kpar)
      KMat[upper.tri(KMat,diag=FALSE)] <- t(KMat)[upper.tri(KMat,diag=FALSE)]
      
    } else {
     if (!is.null(object$x) && !is.null(trndt)) 
       warning("The value of argument trndt will be ignored, since the training data was already suplied by the x component of the value of argument object\n")  
     if (is.null(object$x)) {
       if (is.null(trndt)) 
         stop("Argument trndt (training data) needs to be provided when argument object does not have an x component and argument newdataasKmat is not set to TRUE\n")
              
       if ( !is.null(dim(trndt)) ) {
          if (!is.null(object$center)) {
            if (!is.null(object$scale))  trndt <- scale(trndt,center=object$center,scale=object$scale)
            else trndt <- scale(trndt,center=object$center,scale=FALSE)
          } else {
            if (!is.null(object$scale)) trndt <- scale(trndt,center=FALSE,scale=object$scale)      
          }
       }
       if (is.data.frame(newdata) || is.matrix(newdata)) KMat <- makeKMat0(trndt, newdata, kernel=object$kernel, kpar=object$kpar)
       else KMat <- makeKMat0(trndt,  t(as.matrix(newdata)), kernel=object$kernel, kpar=object$kpar)       
      } else {
        if (is.data.frame(newdata) || is.matrix(newdata)) KMat <- makeKMat0(object$x, newdata, kernel=object$kernel, kpar=object$kpar)
        else KMat <- makeKMat0(object$x,  t(as.matrix(newdata)), kernel=object$kernel, kpar=object$kpar)       
      }
    }
    
  } else {
    if (!is.matrix(KMat))  stop("Wrong value for argument KMat\n")
    if (!is.null(trndt)) 
      warning("The value of argument trndt is ignored when the new data is provided as an already computed kernel matrix\n")  
  }

  if (class(object)[1]=="kernelSVM") {
    factor(object$grplvls[apply(KMat%*%object$alpha,1,which.max)],levels=object$grplvls)
  } else if (class(object)[1]=="kernelSVMwbias") {
      factor(object$grplvls[apply(KMat%*%object$alpha+object$bias,1,which.max)],levels=object$grplvls)
  } else if (class(object)[1]=="kernelksvm") {
#    predict(object$ksvms,as.kernelMatrix(KMat[,SVindex(object$ksvms),drop=FALSE]))
    factor(predict(object$ksvms,as.kernelMatrix(KMat[,SVindex(object$ksvms),drop=FALSE])),levels=object$grplvls)

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
#      res <-  predict(object$ksvms[[1]],as.kernelMatrix(KMat[,SVindex(object$ksvms[[1]]),drop=FALSE]))
#      for (svm in 2:nsvms)  res <- cbind(res,predict(object$ksvms[[svm]],as.kernelMatrix(KMat[,SVindex(object$ksvms[[svm]]),drop=FALSE])))
      res <-  data.frame(factor(predict(object$ksvms[[1]],as.kernelMatrix(KMat[,SVindex(object$ksvms[[1]]),drop=FALSE])),levels=object$grplvls))
      for (svm in 2:nsvms)
        res <- cbind(res,factor(predict(object$ksvms[[svm]],as.kernelMatrix(KMat[,SVindex(object$ksvms[[svm]]),drop=FALSE])),levels=object$grplvls))
    }

    names(res) <- paste("SVM",1:nsvms,sep="")
    res
  } else stop("Wrong class of object argument\n")

}
