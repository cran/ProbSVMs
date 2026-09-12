\name{predict.PVM}

\alias{predict.PVM}

\title{S3 predict method for class 'PVM'}

\description{predict class probabilities based on objects of class 'PVM'.}

\usage{

\method{predict}{PVM}(object, \dots, newdata, eta="adjtogrid", probepsilon="adjtogrid", 
trndt=NULL, retallprd=FALSE, newdataasKmat=FALSE)
}

\arguments{

  \item{object}{Object of class inheriting from \dQuote{PVM}.}

   \item{\dots}{further arguments passed to or from other methods.}

  \item{newdata}{either a matrix or data frame with the data
    (with observations in rows and variables in columns)
    for which class probabilites are to be computed, or a Kernel matrix
    with the values resulting from applying the kernel function to all pairs
    between newdata and training data observations. In the later case argument
    \code{newdataasKmat} should be set to TRUE, and the newdata rows are associated
    with the new observations, while  its columns correspond to training data
    observations.}

   \item{eta}{value of the eta hyper-parameter used in the conversion of SVM(s)
     predictions to class probabilities. It gives the relative importance of
     undesirable versus desirable \eqn{l_1} norm deviations in the goodness of fit
     criterion optimized to find these probabilities. When set to the string \dQuote{adjtogrid}
     (default) it ensures that the undesirable deviations dominate the desirable ones, which is
     necessary for the consistency properties of the resulting estimators. See reference [1] 
     for further details.}

   \item{probepsilon}{the minimum possible value for the class probabilities estimated 
     by \code{predict.PVM}. It could be either the string \dQuote{adjtogrid} (default),
     a small strictly positive constant, or zero. When equal to \dQuote{adjtogrid} 
     it is set to half the minimum weight provided by the grid matrix of weight 
     specifications. This choice is consistent with the resolution implied by this
     matrix. If set to zero, it admits some class probabilites to be estimated
     exactly by zero. This choice implies that the corresponding classes will be impossible
     for the correspondig attribute vector values. See reference [1] for furhter
     discussion.}
     
  \item{trndt}{either a NULL value,  a matrix or data frame with the training data (with
    observations in rows and variables in columns) or the string
    \dQuote{newdata}, if (and only if) class probabilites are to be computed
    for the same data as the one used for training. When equal to NULL the training data
    will be retrieved from the value of argument object (if argument newdataasKmat is set
    to FALSE), or assumed to be already used in the computation of the kernel matrix given 
    by argument newdata (if argument newdataasKmat is set to TRUE) }

  \item{retallprd}{a boolean flag stating if all weigthed SVM predctions
    should be returned together with the estimated class probabilites.}

  \item{newdataasKmat}{a boolean flag stating if the new data is given
    as a pre-computed Kernel matrix.}

}

\value{when \code{retallprd} is set to FALSE, a \code{ClassProb} object with the
  probability estimates for the new data. Class \code{ClassProb} are matrices
  or data frames of class probabilities (observations in rows, classes in columns),
  plus a class attribute, and have a plot method for relevant displays of these
  probabilities.  
  When \code{retallprd} is set to TRUE a list with two components named (i) Pprob:
  a \code{ClassProb} object; and (ii) classpred: a matrix with observations in
  rows and class-weight specifications in columns, containing the weighted SVM
  predictions.}

\seealso{\code{\link{trainPVM}}, \code{\link{plot.ClassProb}} }


\references{

[1] Duarte Silva, A.P. (2025) Probabilistic Vector Machines. \emph{Computers and Operations
Research}, Vol. 183, 107203. <doi:10.1016/j.cor.2025.107203>

}


\examples{

\donttest{

set.seed(1)  # Fix the seed in order to get reproducible results
             # (given the the building ot the grid matrix has a random component)  

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

WWpvm <- trainPVM(crabs[,10:14],crabs$grp)
WWprob <- predict(WWpvm,newdata=crabs[,10:14],trndt="newdata")

# Display the probabilities of the predicted classes

plot(WWprob, , trdata=crabs[,10:14], grouping=crabs$grp)
plot(WWprob, , trdata=crabs[,10:14], grouping=crabs$grp, type="stckdbarplt")
WWprob

# Repeat the analysis, using the first 45 observations in each class for training,
# and the last 5 for probability estimation

trdata <- crabs[c(1:45,51:95,101:145,151:195),10:14]
trgrp <- crabs$grp[c(1:45,51:95,101:145,151:195)]
evaldata <- crabs[c(46:50,96:100,146:150,196:200),10:14]
WWpvm1 <- trainPVM(trdata,trgrp)
WWprob1 <- predict(WWpvm1,newdata=evaldata)

plot(WWprob1, type="stckdbarplt")
WWprob1

}

}
