\name{trainSVM}

\alias{trainSVM}

\title{Efficient training of (multiclass) kernel-based Support Vector Machines}

\description{
trainSVM creates objects of class \code{kernelSVM} or class \code{kernelSVMs}
by training one or several kernel-based (possibly weighted) Support Vector Machines
for general  multiclass problems, using an 'one-in-all' global model approach.
}

\usage{

trainSVM(x=NULL, y, class.weights=rep(1.,length(levels(y))), scaled=TRUE, K=NULL,
         loss=c("WW","LLW"), kernel=c("rbfdot","vanilladot","polydot"),
             kpar=list(sigma="d2median"), C=1., lambda=NULL, keepdt=TRUE, retotpst=FALSE,
         OptCntrl=SetUpOptPar(), TunCntrl=SetUpTunPar() )
}

\arguments{
  \item{x}{a matrix or data frame containing the training data, with observations in rows and variables in columns.
     If x is equal to NULL, a kernel matrix should be provided directly in argument \code{K}.}

  \item{y}{a factor describing the response vector, with one label by observation.}

  \item{class.weights}{a vector or a matrix of class weights.
     When class.weights is a vector, its length shoul equal to number of classes, and the i-th component represents the weight given
     to the i-th class in the SVM loss function. For classical (non-weighted) SVMs all its components are equal to one (default).
     In weighted SVM variants, the class.weights components may be different and should add up to one.
     When \code{class.weights} is a matrix the number of columns should equal the number of classes, and each row corresponds to a
     different weight specification. In that case the resulting SVM is trained as many times as the number of rows.
     When \code{class.weights} is a vector or a single row matrix, \code{trainSVM} returns an object of class \code{kernelSVM} (one trained SVM),
     and when \code{class.weights} is a matrix with multiple rows, \code{trainSVM} returns an object of class \code{kernelSVMs} (several trained
     SVMs).}

  \item{scaled}{A logical flag indicating wheter or not the variables should be scaled to unit variance.}

  \item{K}{A symmetrix kernel matrix. When equal to NULL, the training data should be given by argument \code{x} and the
    kernel matrix is computed from it and from the values of arguments \code{kernel}, and \code{kpar}. When K is not NULL,
    it is assumed to be equal to a pre-computed kernel matrix and the aguments \code{x}, \code{kernel}, and \code{kpar} are ignored.}

  \item{loss}{a string specifying the multiclass large margin loss to be employed. Currently the following alternatives are implemented:
    \dQuote{WW} for the Weston and Watkins loss (see reference [6]), and \dQuote{LLW} for the Lee, Lin and Wahba loss (see reference [4]).
    In two-class problems these two losses are equivalent and also equivalent to the classical hinge loss. In problems
    with more than two classes these losses differ and the \dQuote{LLW} loss has better assymptotic properties, but the \dQuote{WW} loss
    was found empirically to usually give more reliable results in many pratical applications.}

  \item{kernel}{the kernel function used for training and prediction.
    Currentlty this argument can be set to one of the following strings:

    \itemize{
      \item \dQuote{rbfdot} Radial Basis kernel \dQuote{Gaussian}

      \item \dQuote{vanilladot} Linear kernel

      \item \dQuote{polydot} Polynomial kernel
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

  \item{C}{cost of constraints violation (default: 1). This is the "C"-constant of the regularization term in the Lagrange
    formulation. When set to the string \dQuote{tuneit} the value of "C" is tuned acccording
    to the procedure defined by the value of argument \code{TunCntrl}. The value of this argument is overriden if the argument
    \code{lambda} has a valued different than NULL.}

  \item{lambda}{regularization parameter in the functional analysis formulation of SVMs. When set to a non-NULL value, overrides
    the value of argument \code{C}. Specifying the amount of regularization by arguments \code{C} or \code{lambda} leads to equivalent
    models if the corresponding constants satisfy the relation \code{lambda} = 1/(2 n \code{C}),
    where n is the training sample size. When set to the string \dQuote{tuneit} the value of \code{lambda} is tuned acccording
    to the procedure defined by the value of argument \code{TunCntrl}.}

  \item{keepdt}{a logical flag indicating if the orginal training data should be returned together with the trainded SVM(s)}

  \item{retotpst}{a logical flag indicating if internal optimization results and statistics should be returned together with the trainded SVM(s)}

  \item{OptCntrl}{a list with several control arguments for the optimization algorithm used in the SVM(s) training. See \code{\link{SetUpOptPar}} for details. }

  \item{TunCntrl}{a list with several arguments controlling the hyperparameter (\code{C} or \code{lambda}) tuning. Only used when either \code{C} or \code{lambda}
  are set to the string \dQuote{tuneit}. See \code{\link{SetUpTunPar}} for details.}

}


\details{

\code{trainSVM} trains multiclass kernel-based (possibly weighted) Supported Vector Machines (SVMs) without bias terms.
It implements Duarte Silva's adaptation (see reference [3]) of Dogan, Glasmachers and Igel' s efficient training
algorithm (see reference [1]) for 'all-in-one' multiclass SVMs. It can be shown that in kernel-based SVMs dropping
bias terms has considerable computational advantages, and it has been argued that the resulting rules have essentially
the same statistcal properties (at least asymptotically for universal kernels) as the corresponding
rules with bias terms. See references [2], [3] and [5] for further discussion.


Currently there are two variants of multiclass SVMs implemented in trainSVM:
(i) machines that use the SVM loss proposed by Weston and Watkins (see reference [6]), and machines that use
the loss proposed by Lee, Lin and Wahba (see reference [4]). The later variant has better asymptotic properties,
but there are some empirical evidence suggesting that the former often gives more reliable results in many
applications (see, e.g. reference [2]). For two-class problems  these two variants are equivalent among themselves
and to the classical SVM hinge loss.


The amount of regularization of the resulting SVMs is controled by the value of the hyper-paramters \code{C} or \code{lambda}.
While \code{C} (the cost of constraint violations) is the traditional regularization hyper-parameter used in the majority of the SVM literature,
there are also alternative formulations  in which SVM training is presented as a particular case of regularized function estimation
in functional analysis. Under this perspective, the SVM constraint violations are viewed as unit-cost losses, and the margin component
of the traditional SVM criterion is understood as a complexity penalty that is weighted by the regularization hyper-paramter \code{lambda}.
The two formulations are mathematically equivalent if \code{lambda} is set to 1/(2 n \code{C}), or if \code{C} is set to 1/(2 n \code{lambda}),
with n being the number of observations in the training data. Note that, as \code{C} and \code{lambda} are inversely
related, higher values of \code{lambda} and lower values of \code{C} correspond to higher regularization and smoother classification rules,
while lower \code{lambda} values and higher \code{C} values lead to more complex (and flexible) rules, with better training sample performance,
but potentially worse generability ability.

In \code{trainSVM} the values of \code{C} and \code{lambda} can be set by the user, or tuned automatically to the training or tuning data acccording
to the procedure defined by the value of argument \code{TunCntrl}. See the documentation of \code{\link{SetUpTunPar}} for further details.


}

\value{An object inheriting from class \dQuote{kernelSVM} containing the trained Support Vector Machine(s). It could be
either an \dQuote{kernelSVM} object if argument \code{class.weights} is a vector or a single-row
matrix, or an  \dQuote{kernelSVMs} object otherwise. Class \dQuote{kernelSVM} has a predict method,
that gives class predictions from \dQuote{kernelSVM} objects, and data frames or kernel matrices of new (or training) data.

An object inheriting from class \dQuote{kernelSVM} is a list containing, at least, the following components:

    \describe{

      \item{alpha}{in \dQuote{kernelSVM} classes, a matrix (observations in rows, classes in columns) containing the 
        resulting support vectors. In \dQuote{kernelSVMs} classes, a three-dimensional array, in which each level of the third dimension 
        corresponds to a different trained SVM, and the first two dimensions contain the correspondig support vectors, with observations
        in rows and classes in columns.}

      \item{grplvls}{the levels of the factor representing the different classes.}

      \item{lambda}{the value of the \code{lambda} regularization hyperparameter used in the SVM(s) training.}

      \item{C}{the value of the \code{C} regularization hyperparameter used in the SVM(s) training.}

      \item{x}{when argument \code{keepdt} is TRUE, a matrix or data frame containing the training data, with 
              observations in rows and variables in columns. Otherwise, the value NULL.}

      \item{scale}{when argument \code{scaled} is TRUE, a vector with the variables standard deviations, that were used for
              data scaling. Otherwise, the value NULL.}   

      \item{kernel}{the value of argument \code{kernel} in the call to \code{trainSVM}.}

      \item{kpar}{the value of argument \code{kpar} in the call to \code{trainSVM}.}

      \item{optlist}{this component is only non NULL when argument \code{retotpst} is set to TRUE. In that case it is a list with the following components:

          \describe{
          
                \item{optvalue}{the returned value of the optimization model criterion used for SVM training .}
                
                \item{iterations}{the number of iterations used by the SVM training optimization algorithm.}
                
                \item{hitrate}{the classification hit rate (when known) of the trained SVM. Tipically this value is only available when arguments \code{C}
                or \code{lamdba} are set to the string \dQuote{tuneit}, and refers to the optimal hit rate obtained in the tuning data.
                When \code{C} or \code{lambda} is not tuned, this component is set to NULL.}
           }
        }
     } 
  }


\seealso{\code{\link{predict.kernelSVM}}, \code{\link{SetUpOptPar}}, \code{\link{SetUpTunPar}}, \code{\link{GetrbfdotSigPar}}.}

\references{

[1] Dogan, U.; Glasmachers, T. and Igel, C. (2011) Fast training of multi-class support vector machines.
\emph{Technical report. Department of Computer Sciences, University of Copenhagen}.

[2] Dogan U.; Glasmachers T. and Igel, C. (2016) A unified view on multi-class support vector classification.
\emph{Journal of Machine Learning Research}, Vol. 17 (45), 1--32.

[3] Duarte Silva, A.P. (2025) Probabilistic Vector Machines. \emph{Computers and Operations
Research}, Vol. 183, 107203. <doi:10.1016/j.cor.2025.107203>

[4] Lee, Y.; Lin, Y. and Wahba, G. (2004) Multicategory support vector machines:
Theory and application to the classification of microarray data and satellite radiance data.
\emph{Journal of the American Statistical Association}, Vol. 99, 67--81. >doi:10.1198/016214504000000098>

[5] Poggio, T.; Mukherjee, S.;  Rifkin, R.; Raklin, A. and Verri, A. (2002). B. In
\emph{Uncertainty in Geometric Computations}, Springer US.

[6] Weston, J. and Watkins, C. (1999) Support vector machines for multi-class pattern recognition.
In \emph{Proceedings of the 7th European Symposium on Artificial Neural Networks}, 219--224.

}


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


# train the WW and LLW SVMs, with their default setings, or the crabs data

WWsvm <- trainSVM(crabs[,10:14],crabs$grp,loss="WW")
LLWsvm <- trainSVM(crabs[,10:14],crabs$grp,loss="LLW")

# Get in-sample classification results

WWpred <- predict(WWsvm,newdata=crabs[,10:14])
WWpred
LLWpred <- predict(LLWsvm,newdata=crabs[,10:14])
LLWpred

# Compare classifications with true assignments

print(WWpred==crabs$grp)
print(LLWpred==crabs$grp)

\donttest{

# Repeat the analysis, by tuning
# the regularization hyper-paremeter to the training data

WWsvm1 <- trainSVM(crabs[,10:14],crabs$grp,C="tuneit")
WWpred1 <- predict(WWsvm1,newdata=crabs[,10:14])
WWpred1
print(WWpred1==crabs$grp)

LLWsvm1 <- trainSVM(crabs[,10:14],crabs$grp,C="tuneit",loss="LLW")
LLWpred1 <- predict(LLWsvm1,newdata=crabs[,10:14])
LLWpred1
print(LLWpred1==crabs$grp)


# Repeat the analysis, for the WW loss only, using the first 45 observations
# in each class for training, and the last 5 for prediction

trind <- c(1:45,51:95,101:145,151:195)
evalind <- c(46:50,96:100,146:150,196:200)
trdata <- crabs[trind,10:14]
trgrp <- crabs$grp[trind]
evaldata <- crabs[evalind,10:14]
WWsvm2 <- trainSVM(trdata,trgrp,C="tuneit")
WWpred2 <- predict(WWsvm2,newdata=evaldata)
WWpred2
print(WWpred2==crabs$grp[evalind])


# Now, tune C by a cross-validated estimate of the error rate in the training data

WWsvm3 <- trainSVM(trdata,trgrp,C="tuneit",TunCntrl=list(crossval=TRUE))
WWpred3 <- predict(WWsvm3,newdata=evaldata)
WWpred3
print(WWpred3==crabs$grp[evalind])


## Repeat the analysis by using a pre-computed kernel matrix given by the K argument
## (Note: this is recommended in multiple analysis based on the same data, in order
## to avoid uncessary repeated evalutions of the kernel funtion

# All data as training data

  GramMat <- makeKMat(scale(crabs[,10:14],center=FALSE),kernel="rbfdot",kpar=list(sigma="d2median"))
  # Gram matrix, with sigma automatically adjusted to the scaled training data
  
  WWsvm4 <- trainSVM(K=GramMat,y=crabs$grp,C="tuneit")
  WWpred4 <- predict(WWsvm4,KMat=GramMat)
  print(WWpred4==crabs$grp)
  
# Split 45 observations per group for training and 5 observations per group for prediction

  scltrdata <- scale(crabs[trind,10:14],center=FALSE) 
  sclevaldata <- scale(crabs[evalind,10:14],center=FALSE,scale=attr(scltrdata,"scaled:scale")) 
  
  GramMat <- makeKMat(scltrdata,kernel="rbfdot",kpar=list(sigma="d2median"))
  # Gram matrix, with sigma automatically adjusted to the scaled training data
  trnevalKernMat <- makeKMat(scltrdata,sclevaldata,kernel="rbfdot",kpar=list(sigma="d2median"))
  #Kernel matrix betwenn the (commonly scaled) training and evaluation data
    
  WWsvm5 <- trainSVM(K=GramMat,y=crabs$grp[trind],C="tuneit")
  WWpred5 <- predict(WWsvm5,KMat=trnevalKernMat)
  print(WWpred5==crabs$grp[evalind])


}

}
