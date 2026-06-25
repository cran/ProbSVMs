require(ggplot2)

as.ClassProb <- function(x, ...)
{
  if ( !is.matrix(x) && !is.data.frame(x) ) stop("Wrong class type for argument x\n")
  if (!isTRUE(all.equal(rowSums(x),rep(1.,nrow(x))))) stop("Not all x rows sum up to one as they should\n")
  class(x) <- "ClassProb"
  
  x
}

plot.ClassProb <- function(x, ..., type=c("scatterplt","stckdbarplt"), projecton=c("DiscFact","PCs","OrigDt"),
                           trdata=NULL, grouping=NULL, newdata="trdata",
                           ShownCl="PredCl", axis=c(1,2), obs=1:nrow(x), withxlabels=length(obs)<=30,
#                         title=NULL, pntsize = 2.5, threecolors=TRUE, clgrdparl = list(low="white",mid="#fee08b",high="darkpurple",midpoint=0.7)
                           title=NULL, pntsize = 2.5, threecolors=TRUE, clgrdparl = list(low="white",mid="blue",high="darkred",midpoint=0.7)
)
{
  type <- match.arg(type)
  projecton <- match.arg(projecton)
  PredCl <- ClProb <- LDF <-  Observations <- Classes <- Probabilities <- NULL        # dummy command to avoid unecessary NOTES from R CMD check !!!

  if (is.null(clgrdparl$low)) clgrdparl$low <- "white"
  if (is.null(clgrdparl$mid)) clgrdparl$mid <- "blue"
  if (is.null(clgrdparl$high)) clgrdparl$high <- "darkred"
  if (is.null(clgrdparl$midpoint)) clgrdparl$midpoint <- 0.7

  if (type == "stckdbarplt") {

    x1 <- as.data.frame(cbind(Observations=row.names(x[obs,]),x[obs,]))
    x2 <- melt(x1,id.vars="Observations",measure.vars=2:ncol(x1),variable.name = "Classes",value.name="Probabilities")
    x2$Observations <- factor(x2$Observations,levels=unique(x2$Observations),ordered=TRUE)
    x2$Probabilities <- as.numeric(x2$Probabilities)
    p <- ggplot(x2,aes(x=Observations,y=Probabilities,fill=Classes)) + geom_bar(stat="identity") + scale_y_continuous(breaks=seq(0.1,0.9,0.1))
    if (!withxlabels) p <- p + theme(axis.text.x = element_blank())
    else p <- p + theme(axis.text.x=element_text(angle=90,hjust=1))

  } else if (type == "scatterplt") {

    DispTwoLDF <- TRUE

    if ( (projecton!="OrigDt" || !is.null(trdata)) && !is.data.frame(trdata) && !is.matrix(trdata) ) stop("Wrong value for argument trdata\n")
    if ( !is.null(dim(newdata)) || newdata!="trdata" ) {
      if ( !is.data.frame(newdata) && !is.matrix(newdata) ) stop("Wrong value for argument newdata\n")
      if ( !is.null(dim(trdata)) && ncol(newdata)!=ncol(trdata) ) stop("Dimensions of the data given by arguments trdata and newdata do not agree with each other\n")
    }

    if (projecton == "DiscFact") {
      if (is.null(grouping)) stop("Argument grouping is required for a projection into a factorial discriminant space")
      if (ncol(x)<3) {
         DispTwoLDF <- FALSE
         if (is.null(dim(newdata)) && newdata=="trdata") dispdata <- as.data.frame(predict(lda(trdata,grouping))$x[obs])
         else dispdata <- as.data.frame(predict(lda(trdata,grouping),newdata=newdata)$x[obs])
         colnames(dispdata) <- "LDF"
      } else {
        if (is.null(dim(newdata)) && newdata=="trdata") dispdata <- as.data.frame(predict(lda(trdata,grouping))$x[obs,axis])
        else dispdata <- as.data.frame(predict(lda(trdata,grouping),newdata=newdata)$x[obs,axis])
        colnames(dispdata) <- paste("LDF",axis,sep="")
      }
    } else if (projecton == "PCs") {
      if (is.null(dim(newdata)) && newdata=="trdata") dispdata <- as.data.frame(prcomp(trdata,center=TRUE,scale=TRUE)$x[obs,axis])
      else dispdata <- as.data.frame(predict(prcomp(trdata,center=TRUE,scale=TRUE),newdata=newdata)[obs,axis])
      colnames(dispdata) <- paste("PC",axis,sep="")

    }  else if (projecton == "OrigDt") {
      if ( (is.null(dim(newdata)) && newdata=="trdata") || is.null(newdata) ) {
        if ( !is.data.frame(trdata) && !is.matrix(trdata) ) stop("Wrong value for argument trdata\n")
        dispdata <- trdata[obs,axis]
      }
      else {
        if ( !is.data.frame(newdata) && !is.matrix(newdata) ) stop("Wrong value for argument newdata\n")
        dispdata <- newdata[obs,axis]
      }
    }
    dispdata$PredCl <- factor(apply(x[obs,],1,function(row) colnames(x)[which.max(row)]))
    if (ShownCl=="PredCl") {
      dispdata$ClProb <- apply(x[obs,],1,max)
    } else {
      if (!is.numeric(ShownCl) && !is.integer(ShownCl) && !is.element(ShownCl,colnames(x)))
        stop("Wrong value for argument ShownCl\n")
      dispdata$ClProb <- x[obs,ShownCl]
    }
    if (DispTwoLDF)  p <- ggplot(dispdata,aes(x=.data[[colnames(dispdata)[1]]],y=.data[[colnames(dispdata)[2]]],shape=PredCl,colour=ClProb))
    else p <- ggplot(dispdata,aes(x=LDF,y=ClProb,shape=PredCl,colour=ClProb))
    p <- p + geom_point(size=pntsize)

#    if (threecolors) p <- p + scale_color_gradient2(low=clgrdparl$low,mid=clgrdparl$mid,high=clgrdparl$high,midpoint=clgrdparl$midpoint,limits=c(0.,1.))
#    else p <- p + scale_color_gradient(low=clgrdparl$low,high=clgrdparl$high,limits=c(0.,1.))
    if (threecolors) p <- p + scale_color_gradient2(low=clgrdparl$low,mid=clgrdparl$mid,high=clgrdparl$high,midpoint=clgrdparl$midpoint)
    else p <- p + scale_color_gradient(low=clgrdparl$low,high=clgrdparl$high)
  }

  if (!is.null(title)) p <- p + ggtitle(title) + theme(plot.title=element_text(hjust=0.5))
  print(p)

}
