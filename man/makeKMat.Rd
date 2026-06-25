\name{makeKMat}

\alias{makeKMat}

\title{making kernal matrices}

\description{
makeKMat creates kernel matrices from row data. In can be used to create the Gram matrix
from a single data set, or a general kernel matrix from two different data sets.
}

\usage{

makeKMat(dt1, dt2=NULL, kernel=c("rbfdot","vanilladot","polydot"),
  kpar=list(sigma="d2median"))

}

\arguments{
  \item{dt1}{a data matrix or data frame, with observations in rows and variables in columns.}

  \item{dt2}{a data matrix or data frame, with observations in rows and variables in columns.
       If equal to NULL the kernel matrix to be created will consist of kernel values between
       all pairs of the \code{dt1} elements. Otherwise, it will consist of the kernel values between
       the elements of \code{dt1} (in rows) and \code{dt2} (in columns).}

  \item{kernel}{the kernel function used for training and prediction.
    Currentlty this argument can be set to one of the following strings:

    \itemize{
      \item \code{rbfdot} Radial Basis kernel \dQuote{Gaussian}

      \item \code{vanilladot} Linear kernel

      \item \code{polydot} Polynomial kernel
    }
  }

  \item{kpar}{the list of kernel hyper-parameters.
    This is a list which contains the parameters to be used with the
    kernel function. Valid parameters for existing kernels are :

    \itemize{
      \item \code{sigma}: inverse kernel width for the Radial Basis kernel function \dQuote{rbfdot}.
        Either a numeric value or one of the strings \dQuote{d2median} (default), \dQuote{d2q01q09mean"},
        or \dQuote{"d2q01q09hmean"}.  In the later case this value is computed automatically, from the
        distribution of all euclidean distances between \code{dt1}  observation pairs. 
        See \code{\link{GetrbfdotSigPar}} for details.

      \item \code{degree, scale, offset}: for the Polynomial kernel \dQuote{polydot}
    }
  }

}

\value{A kernel matrix. If argument \code{dt2} is NULL, a symmetrix matrix consisting of the kernel values
     between all pairs of \code{dt1} elements. If argument \code{dt2} is not NULL, a matrix with the kernel
     values between the elements of the \code{dt1} (in rows) and \code{dt2} (in columns).}

\seealso{\code{\link{GetrbfdotSigPar}}}

