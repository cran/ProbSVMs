SetUpTunPar <- function(tunex=NULL, tuney=NULL, tuneK=NULL, Csrchpar=list(Cpowerbase=2.,Cgridinlev=7,Cnloops=2),
                        crossval=FALSE, crossvalpar=list(Strfolds=TRUE,kfold=10,CVrep=3))
{
  if ( !is.null(tunex) && !is.matrix(tunex) && !is.data.frame(tunex) ) stop("Wrong type for argument tunex\n")
  if ( !is.null(tuneK) && !is.matrix(tuneK) ) stop("Wrong type for argument tuneK\n")
  if ( !is.null(tuney) && !is.factor(tuney) ) stop("Argument tuney is not a factor as it should\n")

  list(tunex=tunex,tuney=tuney,tuneK=tuneK,Csrchpar=Csrchpar,crossval=crossval,crossvalpar=crossvalpar)
}

CompTunPar <- function(parlist)
{
  if (is.null(parlist$Csrchpar)) parlist$Csrchpar <- list(Cpowerbase=2.,Cgridinlev=7,Cnloops=2)
  if (!is.numeric(parlist$Csrchpar$Cpowerbase) || !is.numeric(parlist$Csrchpar$Cgridinlev) || !is.numeric(parlist$Csrchpar$Cnloops)
      || !is.wholenumber(parlist$Csrchpar$Cgridinlev) || !is.wholenumber(parlist$Csrchpar$Cnloops) )
    stop("Wrong value for the parlist$Csrchpar argument\n")
  if (is.null(parlist$crossval)) parlist$crossval <- FALSE
  if (is.null(parlist$crossvalpar)) parlist$crossvalpar <- list(Strfolds=TRUE,kfold=10,CVrep=3)
  if (is.null(parlist$crossvalpar$Strfolds)) parlist$crossvalpar$Strfolds <- TRUE
  if (is.null(parlist$crossvalpar$kfold)) parlist$crossvalpar$kfold <- 10
  if (is.null(parlist$crossvalpar$CVrep)) parlist$crossvalpar$CVrep <- 3
  if ( !is.logical(parlist$crossvalpar$Strfolds) || !is.numeric(parlist$crossvalpar$kfold) || !is.numeric(parlist$crossvalpar$CVrep) ||
        !is.wholenumber(parlist$crossvalpar$kfold) || !is.wholenumber(parlist$crossvalpar$CVrep)
  )
    stop("Wrong value for the parlist$crossvalpar argument\n")

  parlist
}
