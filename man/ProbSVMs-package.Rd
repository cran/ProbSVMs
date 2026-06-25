\name{ProbSVMs-package}

\alias{ProbSVMs-package}
\alias{ProbSVMs}

\alias{buildgrid}
\alias{GenCrossVFolds}
\alias{evalPP}
\alias{l2norm2}
\alias{GausKern}
\alias{MedianDist1Set}
\alias{MedianDist2Sets}
\alias{LPdual}
\alias{is.wholenumber}
\alias{trainandevalPVM}
\alias{Csearch}
\alias{sigest1}
\alias{makeKMat0}

\docType{package}

\title{Probabilistic Support Vector Machines}

\description{Implements 'all-in-one' kernel-based multiclass Support Vector Machines for supervised classification and class probability estimation.}

\details{
\tabular{ll}{
Package: \tab ProbSVMs\cr
Type: \tab Package\cr
Version: \tab 0.1-0\cr
Date: \tab 2026-06-22\cr
License: \tab GPL-2\cr
LazyLoad: \tab yes\cr
}

Package \pkg{ProbSVMs} implements 'all-in-one' kernel-based multiclass Support Vector Machines for supervised classification and class probability estimation.
The most original novelty of \pkg{ProbSVMs} is the Probabilistic Vector Machines (PVMs) proposed by Duarte Silva in reference [3]. PVMs are distribution-free
estimators of class probabilities based on the predictions given by a sequence of kernel-based weighted Supported Vector Machines (SVMs).

Currently there two variants of multiclass PVMs implemented in \pkg{ProbSVMs}: (i) machines that use the SVM loss proposed by Lee, Lin and Wahba (see reference [5]),
and machine that use the loss proposed by Weston and Watkins (see reference [6]). The former variant has better asymptotic properties, but there are
some empirical evidence suggesting that the later often gives more reliable results in many applications (see, e.g., reference [2]). For two-class problems
these two variants are equivalent.

In problems with more than two classes currently \pkg{ProbSVMs} only implements PVMs based on classification rules derived from weighted kernel-based SVMs without
bias terms. In two-class problems both rules with (default) or without bias terms can be used. In problems
with more than two classes dropping the bias terms has some computational advantages, and it has been argued that it should not have a noticable impact on
statistical performance (see references [2] and [3] for further discussion).

In adition, \pkg{ProbSVMs} also implements Duarte Silva's adaptation (see reference [3]) of Dogan, Glasmachers and Igel' s efficient training
algorithm (see reference [1]) for 'all-in-one' multiclass kernel-based SVMs without bias terms. While an original implementation of this
algorithm is availalble on the machine learning Shark library (see reference [4]), unlike \pkg{ProbSVMs} that implementation does not seems to cover
weighted SVM versions, nor does it seem to include any interface to the R environment.

in \pkg{ProbSVMs}, PVMs can be trained by the function \code{\link{trainPVM}} which generates an object of class \code{PVM}, which has a predict method for computing
reliable class probability estimates. Similarly, kernel-based SVMs can be trained by the function \code{\link{trainSVM}} which generates an object of class
\code{kernelSVM} (single trained SVM)) or \code{kernelSVMs} (multiple trained SVMs) that have predict methods for finding class predictions.

}

\author{
Antonio Pedro Duarte Silva  <psilva@ucp.pt>

Maintainer: Antonio Pedro Duarte Silva  <psilva@ucp.pt>
}

\seealso{\code{\link{trainPVM}}, \code{\link{trainSVM}}, \code{\link{predict.PVM}}, \code{\link{predict.kernelSVM}}, \code{\link{plot.ClassProb}}}

\references{


[1] Dogan, U.; Glasmachers, T. and Igel, C. (2011) Fast training of multi-class support vector machines.
\emph{Technical report. Department of Computer Sciences, University of Copenhagen}.

[2] Dogan U.; Glasmachers T. and Igel, C. (2016) A unified view on multi-class support vector classification.
\emph{Journal of Machine Learning Research}, Vol. 17 (45), 1--32.

[3] Duarte Silva, A.P. (2025) Probabilistic Vector Machines. \emph{Computers and Operations
Research}, Vol. 183, 107203. <doi:10.1016/j.cor.2025.107203>

[4] Igel, C.; Heidrich-Meisner, V. and Glasmachers, T. (2008) Shark. \emph{Journal of Machine Learning Research},
Vol 9 (6), 993-996.

[5] Lee, Y.; Lin, Y. and Wahba, G. (2004) Multicategory support vector machines:
Theory and application to the classification of microarray data and satellite radiance data.
\emph{Journal of the American Statistical Association}, Vol. 99, 67--81. <doi:10.1198/016214504000000098>

[6] Weston, J. and Watkins, C. (1999) Support vector machines for multi-class pattern recognition.
In \emph{Proceedings of the 7th European Symposium on Artificial Neural Networks}, 219--224.

}


\keyword{ProbSVMs}

\examples{

# train the Weston and Watkins SVM on the Iris data


data(iris)
WWsvm <- trainSVM(iris[,1:4],iris$Species,loss="WW")

# Get in-sample classification results

WWpred <- predict(WWsvm,newdata=iris[,1:4])
WWpred

# Compare classifications with true assignments

cat("Original classes:\n")
print(iris$Species)
print(WWpred==iris$Species)

# Estimate class probabilities based on the WW loss

WWpvm <- trainPVM(iris[,1:4],iris$Species,loss="WW")
WWprob <- predict(WWpvm,newdata=iris[,1:4],trndt="newdata")

# Display the probabilities of the predicted classes
# on the space of the petal measurements

plot(WWprob, projecton="OrigDt", trdata=iris, axis=c(3,4))

}
