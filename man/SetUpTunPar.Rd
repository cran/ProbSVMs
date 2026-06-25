\name{SetUpTunPar}

\alias{SetUpTunPar}

\title{Setting up SVM tuning parameters.} 

\description{\code{SetUpTunPar} sets up several control parameters in the 
procedure used for tuning SVM regularization hyper-parameters.}


\usage{

SetUpTunPar(tunex=NULL, tuney=NULL, tuneK=NULL, 
Csrchpar=list(Cpowerbase=2.,Cgridinlev=7,Cnloops=2), crossval=FALSE, 
crossvalpar=list(Strfolds=TRUE,kfold=10,CVrep=3))

}

\arguments{

  \item{tunex}{a matrix or data frame containing the data (observations in rows
    and variables in columns) used for tuning SVM hyper-parameters. It is ignored 
    if argument \code{tuneK} is not equal to NULL.}
     
   \item{tuney}{a factor describing the response vector (with one label by observation) 
    used or tuning SVM hyper-parameters.}  

  \item{tuneK}{a kernel matrix based on the trainining data (columns), and the data used for tuning 
   SVM hyper-parameters (rows). When not NULL, overrides the value of argument \code{tunex}.}

   \item{Csrchpar}{list of control parameters used for the search of the tuning parameters.
    It is formed by the following components:    
    
    \itemize{
    
      \item Cpowerbase: base of sucessive powers defining the finite set of values to to searched.

  \item Cgridinlev: number of powers of \code{Cpowerbase} in each loop.

  \item Cnloops: number of times a search over \code{Cgridinlev} points is performed. Each aditional
    loop moves from wider searches to more compact searches.

    }
  }

  \item{crossval}{a logical flag flag indicating if the tuning procedure should be 
    based on a (statistically more sound, but more time intensitive) cross-validation
    strategy.}  

  \item{crossvalpar}{list of control parameters used for the cross-validation strategy.
    It is ignored when crossval is set to FALSE (default), and is formed by the following components:
    
    
    \itemize{
    
      \item Strfolds: a logical flag indicating if the folds should be stratified according to the original class
proportions (default), or randomly generated from the whole training sample, ignoring class membership.

  \item kfold: the number of training sample folds to be created in each replication.

  \item CVrep: the number of replications to be performed.

    }
  }
 
}

\details{
  \code{SetUpTunPar} sets up several control parameters in the  procedure used for tuning SVM regularization hyper-parameters.

  The tuning procedure is based on a search for the optimal value of an performance criterion over a finite grid of 
  \code{Csrchpar$Cgridinlev} different powers of \code{Csrchpar$Cpowerbase}, centered at zero. 
  
  When argument \code{Csrchpar$Cnloops}  is set to one, these powers are consecutive integers. When \code{Csrchpar$Cnloops} is higher 
  than one, there is an initial search over a wider grid, which is then refined around the best value found in the previous loop. 
  This procedure is repeated \code{Csrchpar$Cnloops}  times, in such way that in the last loop there are \code{Csrchpar$Cgridinlev} 
  consecutive integers.         

  The evaluation criterion to be optimized can be an SVM error rate on its training data (if arguments \code{tunex} and \code{tuneK} are both
  NULL), on a tuning data set defined by arguments \code{tunex} and \code{tuney} or \code{tuneK}  (if argument \code{crossval} is FALSE), or a 
  cross-validated estimate of the SVM error rate (if argument \code{crossval} is TRUE) as defined by argument \code{crossvalpar}.
}

\value{a list with components \code{tunex}, \code{tuney}, \code{tuneK}, \code{Csrchpar}, \code{crossval} and \code{crossvalpar} .}

\seealso{\code{\link{trainSVM}}}
  
