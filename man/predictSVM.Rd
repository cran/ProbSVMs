\name{predict.kernelSVM}

\alias{predict.kernelSVM}

\title{S3 predict method for class 'kernelSVM'}

\description{predict classes based on object of class 'kernelSVM'.}


\usage{

\method{predict}{kernelSVM}(object, ..., newdata=NULL, KMat=NULL, trndt=NULL)
}


\arguments{

  \item{object}{Object of class inheriting from \dQuote{kernelSVM}.
    This could either contain a single trained SVM (if the \code{object}
    class is \dQuote{kernelSVM}) or multiple trained SVMs (if the
    \code{object} class is \dQuote{kernelSVMs}.}

   \item{\dots}{further arguments passed to or from other methods.}

  \item{newdata}{a matrix or data frame with the new data
    (with observations in rows and variables in columns)
    to be predicted by the SVM(s) given by \code{object}.
    It is ignored when argument \code{KMat} is not NULL.}

  \item{KMat}{a Kernel matrix with the values resulting from applying
    the kernel function to all pairs between newdata and training data
    observations. The \code{KMat} rows are associated with the new
    observations, while  its columns correspond to training data
    observations. When not NULL overrides the \code{newdata} argument.}

  \item{trndt}{either a matrix or data frame with the training data (with
    observations in rows and variables in columns) or the string
    \dQuote{newdata}, if (and only if) SVM predictions are to be found
    for the same data as the one used for training.}

}

\value{For predictions based on single SVM objects (classe \dQuote{kernelSVM}),
 a factor with the class predictions for each new data observation. For predictions
  based on multiple SVM objects (class \dQuote{kernelSVMs}), a data frame of factors
  where which each column corresponds to a different SVM and each row corresponds
  to a new data observation.}

\seealso{\code{\link{trainSVM}}}


\examples{

## The following examples are based on the MASS data set "crabs".
# This data records physical measurements on 200 specimens of
# Leptograpsus variegatus crabs observed on the shores
## of Western Australia. The crabs are classified by two factors, sex and sp
## (crab species as defined by its colour: blue or orange), with two levels
## each. The measurement variables include the natural logarithms of carapace length (CL),
## the carapace width (CW), the size of the frontal lobe (FL) and the size of
## the rear width (RW). In the  analysis provided, we created four classes
# by crossing the sex and sp levels.

library(MASS)
data(crabs)

crabs$grp <- interaction(crabs$sex,crabs$sp)

crabs$lnFL <- log(crabs$FL)
crabs$lnRW <- log(crabs$RW)
crabs$lnCL <- log(crabs$CL)
crabs$lnCW <- log(crabs$CW)
crabs$lnBD <- log(crabs$BD)

# train the WW SVM, with its default setings, or the crabs data


WWsvm <- trainSVM(crabs[,10:14],crabs$grp)

# Get in-sample classification results

WWpred <- predict(WWsvm,newdata=crabs[,10:14])
WWpred

# Compare classifications with true assignments

cat("Original classes:\n")
print(crabs$grp)
print(WWpred==crabs$grp)

}

