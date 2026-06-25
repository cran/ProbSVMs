\name{as.ClassProb}

\alias{as.ClassProb}

\title{Coerce matrix and data.frames to a \code{ClassProb} object}

\description{conversion function to create \code{ClassProb} objects form matrices and data.frames.}


\usage{

as.ClassProb(x, \dots)

}

\arguments{

  \item{x}{a matrix or data frame (with observations in rows and classes in columns) containing class probability estimates
   in a supervised classification problem.}

 \item{\dots}{further arguments passed to or from other methods.}

}

\value{an object of class \code{ClassProb} containing class probability estimates, for which there is a specialized plot method.
The internal structure of objects inheriting from class \code{ClassProb} is simply the original matrix or data frame plus its
class attribute.}

\seealso{ \code{\link{plot.ClassProb}}}
