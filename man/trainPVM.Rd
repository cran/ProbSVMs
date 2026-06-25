\name{trainPVM}

\alias{trainPVM}

\title{Training of Probabilistic Vector Machines}

\description{
trainPVM creates objects of class 'PVM' (Probabilistic Vector Machines)
by training a sequence of kernel-based weighted Supported Vector Machines
with different class-weight specifications.
}

\usage{

trainPVM(x=NULL, y, scaled=TRUE, K=NULL, loss=c("WW","LLW"),
         withbias=(length(levels(y))==2),
           kernel=c("rbfdot","vanilladot","polydot"), kpar=list(sigma="d2median"),
         C=0.25, lambda=NULL, tunex=x, tuney=y, tuneK=K,
         grid=NULL, dpiinv=ceiling(sqrt(length(y))/0.2), keepdt=TRUE, \dots )
}

\arguments{
  \item{x}{a matrix or data frame containing the training data, with observations in rows and variables in columns.
     If x is equal to NULL, a kernel matrix should be provided directly in argument \code{K}.}

  \item{y}{a factor describing the response vector, with one label by observation.}

  \item{scaled}{A logical flag indicating wheter or not the variables should be scaled to unit variance.}

  \item{K}{A symmetrix kernel matrix. When equal to NULL, the training data should be given by argument \code{x} and the
    kernel matrix is computed from it, and from the values of arguments \code{kernel} and \code{kpar}. When K is not NULL,
    it is assumed to be equal to a pre-computed kernel matrix and the arguments \code{x}, \code{kernel}, and \code{kpar} are ignored.}

  \item{loss}{a string specifying the multiclass large margin loss to be employed. Currently the following alternatives are implemented:
    \dQuote{LLW} for the Lee, Lin and Wahba loss (see reference[3]), and \dQuote{WW} for the Weston and Watkins loss (see reference[4]).
    In two-class problems these two losses are equivalent, and also equivalent to the classical hinge loss.}

  \item{withbias}{a logical flag indicating if a bias term (intercept) should be included in the classification rule. For problems with
    more than two classes, currently only rules without bias terms are implemented. However, in two-class problems the default is to use rules
    with a bias term.}

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
        or \dQuote{"d2q01q09hmean"}.  In the three later cases this value is computed automatically, from the
        distribution of all euclidean distances between \code{dt1}  observation pairs. 
        See \code{\link{GetrbfdotSigPar}} for details.

      \item \code{degree, scale, offset}: for the Polynomial kernel \dQuote{polydot}
    }
  }

  \item{C}{cost of constraints violation in the base SVMs (default: 0.25). This is the "C"-constant of the regularization term in the Lagrange
    formulation. When set to the string \dQuote{tuneit} the value of "C" is tuned to the data provided
    in arguments \code{tuneK} (default), or \code{tunex} and \code{tuney} (if \code{tuneK} is set to NULL).
    The value of this argument is overriden if the argument \code{lambda} has a value different than NULL.}

  \item{lambda}{regularization parameter in the functional analysis formulation of SVMs. When set to a non-NULL value, overrides
    the value of argument \code{C}. Specifying the amount of regularization by arguments \code{C} or \code{lambda} leads to equivalent
    models if the corresponding constants satisfy the relation \code{lambda} = 1/(2 n \code{C}),
    where n is the training sample size. When set to the string \dQuote{tuneit} the value of \code{lambda} is tuned to the data provided
    in arguments \code{tuneK} (default), or \code{tunex} and \code{tuney} (if \code{tuneK} is set to NULL).}

  \item{tunex}{a matrix or data frame containing the data (observations in rows and variables in columns) used to tune the hyper-parameters
     \code{C} and \code{lambda}. It is ignored if none of these arguments is set to the string \dQuote{tuneit}, or if argument \code{tuneK}
     is not equal to NULL.}

  \item{tuney}{a factor describing the response vector (with one label by observation) used to tune the hyper-parameters
     \code{C} and \code{lambda}. It is ignored if none of these arguments is set to the string \dQuote{tuneit}.}

  \item{tuneK}{a kernel matrix based on the trainining data (columns), and the data used (rows) for tuning the hyper-parameters \code{C} and \code{lambda}. 
  It is ignored if none of these arguments is set to the string \dQuote{tuneit}. When set to a non-NULL value overrides the value of argument \code{tunex}.}

  \item{grid}{A matrix with the class-weight specifications used to train the PVM. The rows correspond to different specifications, and the columns
    to the corresponding class weights. All rows should have non-negative elements adding up to one. When set to NULL (default) the grid matrix
    is automatically created based on the value of argument \code{dpiinv}. }

  \item{dpiinv}{the distance between two consecutive weights in the (one-dimensional) specifications of class-weights used to construct a grid matrix.
    Note that weights are uniformily distributed only along one dimension, while the remaing components of the grid are set at random. However, all grid
    dimensions (classes) are  chosen in turn as the one forced to have uniform weights. See reference [2] for further details.
    This argument is ignored when argument \code{grid} is not NULL.}

  \item{keepdt}{a logical flag indicating if the orginal training data should be returned together with the trainded PVM}

  \item{\dots}{other arguments to be passed to trainSVM}

}

\details{

trainPVM trains the Probabilistic Vector Machines (PVMs) proposed by Duarte Silva in reference [2]. PVMs are distribution-free estimators of class
probabilities based on the predictions given by a sequence of kernel-based weighted Supported Vector Machines (SVMs) with different class-weight specifications.
A grid matrix with these specifications is usually created automatically, with its resolution level controled by the argument  \code{dpiinv}, but can also
be directly provided by argument \code{grid}.

Currently there two variants of multiclass PVMs implemented in trainPVM: (i) machines that use the SVM loss proposed by Lee, Lin and Wahba (see reference [3]),
and machine that use the loss proposed by Weston and Watkins (see reference [4]). The former variant has better asymptotic properties, but there are
some empirical evidence suggesting that the later often gives more reliable results in many applications (see, e.g., reference [1]). For two-class problems
these two variants are equivalent.

In problems with more than two classes currently trainPVM only implements PVMs based on classification rules derived from weighted kernel-based SVMs without
bias terms. In two-class problems rules with (default) or without bias terms can be used, with this choice controled by the argument \code{with bias}. In problems
with more than two classes dropping the bias terms has some computational advantages, and it has been argued that it should not have a noticable impact on
statistical performance (see references [1] and [2] for further discussion).

The amount of regularization provided by the underlying SVMs is controled by the value of the hyper-paramters \code{C} or \code{lambda}. While \code{C} (the
cost of constraint violations) is the traditional regularization hyper-parameter used in the majority of the SVM literature, there are also alternative formulations
in which SVM training is presented as a particular case of regularized function estimation in functional analysis. In this perspective,
the SVM constraint violations are viewed as unit-cost losses, and the margin component of the traditional SVM criterion is understood as a complexity penalty
that is weighted by the regularization hyper-paramter \code{lambda}. The two formulations are mathematically equivalent if \code{lambda} is set to 1/(2 n \code{C}),
or if \code{C} is set to 1/(2 n \code{lambda}), with n being the number of observations in the training data. Note that, as \code{C} and \code{lambda} are inversely
related, higher values of \code{lambda} and lower values of \code{C} correspond to higher regularization and smoother classification rules, while lower \code{lambda}
values and higher \code{C} values lead to more complex (and flexible) rules with better training sample performance, but possibly with worse generability potential.

In trainPVM the values of \code{C} and \code{lambda} can be set by the user, or tuned automatically to the training or additional data, provided by arguments \code{tunex},
\code{tuney}, or \code{tuneK}. In that case, \code{C} and \code{lambda} are found by optimizing the log-likelihood of the class probability estimates on the tunning data.
However, the tuning procedure is computionally intensive, and can be very time consuming. See reference [2] for further details.

}

\value{A object of class \dQuote{PVM} containing the trained Probability Vector Machine. Class \dQuote{PVM} has a predict method, that estimates class probabilites
from \dQuote{PVM} objects, and data frames or Kernel matrices of new (or training) data.

An object inheriting from class \dQuote{PVM} is a list containing at least the following components:

    \describe{
    
      \item{nclasses}{the number of different classess.}
                         
      \item{kernel}{the value of argument \code{kernel} in the call to \code{trainPVM}.}
      
      \item{kpar}{the value of argument \code{kpar} in the call to \code{trainPVM}.}
      
      \item{grid}{the grid matrix with the class-weight specifications used to train the PVM.}

      \item{svms}{a list with as many components as the number of different weighted SVMs trained.
       Each component of this list is an object of class \code{kernelSVM} (if argument \code{withbias} is FALSE) or class 
       \code{kernelksvms} (if argument \code{withbias} is TRUE) containing the trained SVM for a particular weight specification.} 
       
     }

  }

\seealso{\code{\link{predict.PVM}}, \code{\link{trainSVM}}, \code{\link{SetUpOptPar}}, \code{\link{GetrbfdotSigPar}}, \code{\link{plot.ClassProb}}}

\references{

[1] Dogan U.; Glasmachers T. and Igel, C. (2016) A unified view on multi-class support vector classification.
\emph{Journal of Machine Learning Research}, Vol. 17 (45), 1--32.

[2] Duarte Silva, A.P. (2025) Probabilistic Vector Machines. \emph{Computers and Operations
Research}, Vol. 183, 107203. <doi:10.1016/j.cor.2025.107203>

[3] Lee, Y.; Lin, Y. and Wahba, G. (2004) Multicategory support vector machines:
Theory and application to the classification of microarray data and satellite radiance data.
\emph{Journal of the American Statistical Association}, Vol. 99, 67--81. <doi:10.1198/016214504000000098>

[4] Weston, J. and Watkins, C. (1999) Support vector machines for multi-class pattern recognition.
In \emph{Proceedings of the 7th European Symposium on Artificial Neural Networks}, 219--224.

}

\examples{

\donttest{

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

# Estimate class probabilities based on the WW loss

WWpvm <- trainPVM(crabs[,10:14],crabs$grp)
WWprob <- predict(WWpvm,newdata=crabs[,10:14],trndt="newdata")


# Display the probabilites of the predicted classes

plot(WWprob, , trdata=crabs[,10:14], grouping=crabs$grp)
plot(WWprob, , trdata=crabs[,10:14], grouping=crabs$grp, type="stckdbarplt")
WWprob

# Repeat the analysis, using the first 45 observations in each class for training,
# and the last 5 for probability estimation

trind <- c(1:45,51:95,101:145,151:195)
evalind <- c(46:50,96:100,146:150,196:200)
trdata <- crabs[trind,10:14]
trgrp <- crabs$grp[trind]
evaldata <- crabs[evalind,10:14]
WWpvm1 <- trainPVM(trdata,trgrp)
WWprob1 <- predict(WWpvm1,newdata=evaldata,trndt=trdata)

plot(WWprob1, type="stckdbarplt")
WWprob1

}

}
