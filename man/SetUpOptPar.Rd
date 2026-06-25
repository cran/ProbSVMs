\name{SetUpOptPar}

\alias{SetUpOptPar}

\title{Set up optimization parameters.} 

\description{\code{SetUpOptPar} sets up several control parameters for the optimization of the mathematical
 programming models used for SVM training.}


\usage{

SetUpOptPar(start=c("allub","alllb","compubwithlb","warmstarts"), alpha0=NULL,
                        epsilon=1e-12, maxiter=100000, tol=1.5e-12 )
}

\arguments{

  \item{start}{string describing the initial solution for the optimization algorithm. It
  can have one of the follwing values: 
  
    \itemize{
      \item \dQuote{allub}: all variables are initialized at their upper bounds.

      \item \dQuote{alllb}: all variables are initialized at their lower bounds. 

      \item \dQuote{compubwithlb}: the solutions with all variables at their lower bounds and
      with all variables at their upper bounds are compared, and the one with the best (lower)
      criterion value is selected as the initial solution. 
      
      \item \code{warmstarts} the initial solution is set at the values given by argument \code{alpha0}
    }
  }
     
   \item{alpha0}{vector with initial values for the optimization algorithm. Only used
   wehn argument \code{start} is equal to \dQuote{warmstarts}.}  

  \item{epsilon}{stoping criterion for the optimization algorithm. When gradient values, or criteria increments,
    are below \code{epsilon} the algorithm stops, and the current solution is returned as the pressumed optimum.}

   \item{maxiter}{maximum number of iterations for the optimization algorithm.}

  \item{tol}{variable numerical tolerance. Variables closer to their lower or upper bounds than \code{tol}
    are assumed to be exactly at the corresponding bounds.}  

}

\value{a list with components \code{start}, \code{alpha0}, \code{epsilon}, \code{maxiter} and \code{tol} .}

\seealso{\code{\link{trainSVM}}}
  
