\name{plot.ClassProb}

\alias{plot.ClassProb}

\title{S3 plot method for \code{ClassProb} objects}

\description{Method for displaying the class probabilities contained in \code{ClassProb} objects.}

\usage{

\method{plot}{ClassProb}(x, \dots, type=c("scatterplt","stckdbarplt"),
 projecton=c("DiscFact","PCs","OrigDt"), trdata=NULL, grouping=NULL,
 newdata="trdata", ShownCl="PredCl", axis=c(1,2), obs=1:nrow(x),
 withxlabels=length(obs)<=30, title=NULL, pntsize = 2.5, threecolors=TRUE,
 clgrdparl = list(low="white",mid="blue",high="darkred",midpoint=0.7))

}


\arguments{

  \item{x}{an object of class inheriting from \dQuote{ClassProb}.}

 \item{\dots}{further arguments passed to or from other methods.}

 \item{type}{type of plot to be displayed.
     Currentlty this argument can be set to one of the following strings:

    \itemize{
      \item \dQuote{scatterplt} scatter plot of projections on a space defined by the value of argument \code{projecton}.
       Probabilities are represented by different colors according to the colour scale defined by the value of arguments
       \code{threecolors} and  \code{clgrdparl}.

      \item \dQuote{stckdbarplt} stacked bar plot. Probabilites are represented by the size of different colour groups.
    }
  }

  \item{projecton}{Type of projection used to display the class probabilites. The value of this argument is ignored when 
  argument \code{type} is equal to \dQuote{stckdbarplt}. Currentlty this argument can be set to one of the following strings:

    \itemize{

      \item \dQuote{DiscFact} Discriminant space formed by two (problems with three o more classes) or the unique
      (problems with two classes) canonical Linear Discriminant Factors (LDFs). With problems with more than three
      groups the set of LDFs chosen is determined by the value of argument \code{axis}. By default the first two
      LDFs are chosen.

      \item \dQuote{PCs} Principal space formed by two Principal Components. The set of PCs chosen is determined
      by the value of argument \code{axis}. By default the first two PCs are chosen.

      \item \dQuote{OrigDt} Space of the original variables corresponding to \code{trdata} columns determined by
      the value of argument \code{axis}.
    }
  }

  \item{trdata}{either a matrix or data frame with the training data (with observations in rows and variables in columns)
  based on which the scatter plot projections are to be found. The value of this argument is ignored when argument \code{type} is equal
  to \dQuote{stckdbarplt}.}

  \item{grouping}{a factor describing the response vector, with one label by observation, in the training data. This argument
  is only required for projections onto canonical Linear Discriminant spaces, and its value is ignored when argument argument \code{type}
  is equal to \dQuote{stckdbarplt}, or argument \code{projecton} is different from \dQuote{DiscFact}.}

  \item{newdata}{either a matrix or data frame with the data (with observations in rows and variables in columns)
    for which class probabilites are to be displayed.}

  \item{ShownCl}{a description of the class whose probabilities are to be displayed by a colour grid in the scatter plot projection
  displays. This could be either the string \dQuote{PredCl} (default), if for each observation the probabilities of the predicted
  clase are to be displayed, or the position or name of a particular class. The value of this argument is ignored when argument \code{type}
  is equal to \dQuote{stckdbarplt}.}

  \item{axis}{the set of LDFs or PCs (when argument \code{projecton} is equal to \dQuote{DiscFact} or \dQuote{PCs}) or column positions
  (when argument \code{projecton} is equal to \dQuote{OrigDt}) used to define the projections in scatter plot displays. The value of this
  argument is ignored when argument \code{type} is equal to \dQuote{stckdbarplt}.}

  \item{obs}{the set of \code{newdata} observations for which class probabilities are to be displayed. Its value can be either equal to
  the row indices or row names of the data provided by argument \code{newdat}.}

  \item{withxlabels}{a logical flag indicating if, in stacked bar plot displays, the observation names are to be shown on the x axis.
   The value of this argument is ignored when argument \code{type} is equal to \dQuote{scatterplt}.}

  \item{title}{an overall title for the plot.}

  \item{pntsize}{a numerical value giving the amount by which, in scatter plots, plotting symbols should be magnified relative
  to the \code{base} \code{plot} default. The value of this argument is ignored when argument \code{type} is equal to \dQuote{stckdbarplt}.}

  \item{threecolors}{logical flag indicating if, in scatter plot displays, the colour scheme used for representing class probabilities
  should be based on a  three-colour scale (TRUE), or a more conventional (but somewhow less discriminatory) gradient colour scale (FALSE).
  The value of this argument is ignored when argument \code{type} is equal to \dQuote{stckdbarplt}.}

  \item{clgrdparl}{a list of arguments used to define the colour scale used for representing class probabilities in scatter plot projections.
  The value of this argument is ignored when argument \code{type} is equal to \dQuote{scatterplt}. The components of this list are :

    \itemize{

      \item \code{low}: the colour associated with low class probabilities.

      \item \code{mid}: the colour associated with medium class probabilities. The value of this argument is ignored when argument \code{threecolors}
      is equal to FALSE.

      \item \code{high}: the colour associated with high class probabilities.

      \item \code{midpoint}: the medium class probability associated with the colour given by the \code{mid} component.
      The value of this argument is ignored when argument \code{threecolors} is equal to FALSE.

    }
  }
}

\value{No return value, called for side effects}

\seealso{ \code{\link{trainPVM}}, \code{\link{predict.PVM}}, \code{\link{as.ClassProb}} }

\examples{

# Train the Weston and Watkins PVM on the Iris data 
# and obtain the training sample class probabilities


data(iris)

WWpvm <- trainPVM(iris[,1:4],iris$Species)
WWprob <- predict(WWpvm,newdata=iris[,1:4],trndt="newdata")

# Display the probabilites of the predicted classes
# on the two-dimensional canonical discriminant space

plot(WWprob, trdata=iris[,1:4], grouping=iris$Species)

# Display the same probabilities on the space of the petal measurements

plot(WWprob, projecton="OrigDt", trdata=iris, axis=c(3,4))

# Show all class probabilities by a stacked bar plot

plot(WWprob, type="stckdbarplt")

\donttest{

# Display the scatter plots, focusing now on the virginica probabilities

plot(WWprob, trdata=iris[,1:4], grouping=iris$Species, ShownCl="virginica")
plot(WWprob, projecton="OrigDt", trdata=iris, axis=c(3,4), ShownCl="virginica")

# Repeat the previous analysis, using the first 40 observations in each class
# for training, and the last 10 for probability estimation

trdata <- iris[c(1:40,51:90,101:140),1:4]
trSpecies <- iris$Species[c(1:40,51:90,101:140)]
evaldata <- iris[c(41:50,91:100,141:150),1:4]
WWpvm1 <- trainPVM(trdata,trSpecies)
WWprob1 <- predict(WWpvm1,newdata=evaldata)

plot(WWprob1, projecton="DiscFact", trdata=trdata,  newdata=evaldata, grouping=trSpecies)
plot(WWprob1, projecton="OrigDt", newdata=evaldata)
plot(WWprob1, type="stckdbarplt", title = "Estimated class probabilities")
plot(WWprob1, projecton="DiscFact", trdata=trdata,  newdata=evaldata,
 grouping=trSpecies, ShownCl="virginica")
plot(WWprob1, projecton="OrigDt", newdata=evaldata, ShownCl="virginica")

}

}
