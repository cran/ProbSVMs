\name{GetrbfdotSigPar}

\alias{GetrbfdotSigPar}

\title{rbfdot Sigma hyper-parameter} 

\description{
\code{GetrbfdotSigPar} returns sensible values for the sigma parameter of the rbfdot (\dQuote{Gaussian}) kernel function.
}

\usage{

GetrbfdotSigPar(x, kpar=c("d2median","d2q01q09mean","d2q01q09hmean"))
}

\arguments{
        
  \item{x}{a matrix or data frame of training data (observations in rows, variables in columns).}
    
  \item{kpar}{a string specifying how sigma should be computed. Valid alternatives are:

    \itemize{
      \item \dQuote{d2median}: the median ofj the the euclidean  distances between all observations pairs in the data given by argument \code{x}.

      \item \dQuote{d2q01q09mean}: the arithmetic mean of the 10-th and 90-th percentiles of the the euclidean  distances between all observations pairs in the data given by argument \code{x}. 

      \item \dQuote{d2q01q09hmean}: the harmonic mean of the 10-th and 90-th percentiles of the the euclidean  distances between all observations pairs in the data given by argument \code{x}.
    }
  }  
   
}

\details{

\code{GetrbfdotSigPar} returns a sensible value for the sigma parameter of the Radial Basis (\dQuote{Gaussian}) kernel function.
According to Caputo  et. al (reference [1]) any sigma values between the 10-th and 90-th percentiles of the distribution of the 
pariwise euclidean distances between observation pairs, \eqn{d_{ij}^2 = ||x_i - x_j||^2}, lead to reasonable Support Vector Machine 
classification results. \code{GetrbfdotSigPar} returns one of three alternatives for the sigma value: (i) the median of the 
\eqn{d_{ij}^2} distribution (default). (ii) the arithmetic mean of the 10-th and 90-th percentiles
of the \eqn{d_{ij}^2} distribution. (iii) the harmonic mean of the 10-th and 90-th percentiles of the \eqn{d_{ij}^2} distribution.  
             
}

\value{a scalar with a sensible value for sigma hyper-parameter of the Radial Basis (\dQuote{Gaussian}) kernel function.}

\seealso{\code{\link{trainPVM}}, \code{\link{trainSVM}}.}

\references{

[1] Caputo, B.; Sim, K.; Furesjo, F. and Smola, A. (2002). Appearance-based
object recognition using svms: which kernel should i use? In \emph{Proceedings of NIPS
workshop on Statistical methods for computational experiments in visual
processing and computer vision}, Whistler.

}
