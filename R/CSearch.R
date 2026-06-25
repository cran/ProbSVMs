Csearch <- function(Cgridinlev=7,nloops=2,powerbase=2,zereoexp=c("center","start"),algorithm,...)
{
  zereoexp <- match.arg(zereoexp)   

  validres <- matrix(nrow=Cgridinlev,ncol=nloops)
  Cgrid <- numeric(Cgridinlev)
  minb <- floor(Cgridinlev/2)
  meadind <- minb+1 
  bestb <- meadind
  if (zereoexp=="start") bestC <- powerbase^(meadind-1) 
  else bestC <- 1.
  bestres <- -Inf

  for (loop in 1:nloops) {
    if (loop==1 && zereoexp=="start") bind2 <- 1
    else bind2 <- -minb
    for (bind1 in 1:Cgridinlev) {
      if (loop==1) {
        if (zereoexp=="start") {
          if (nloops==1) Cgrid[bind1] <- powerbase^(bind2-1)
          else Cgrid[bind1] <- powerbase^(bind2*meadind-1)
        } else {
          if (nloops==1) Cgrid[bind1] <- powerbase^bind2
          else Cgrid[bind1] <- powerbase^(bind2*meadind)
        } 
      }
      else if (loop==2) Cgrid[bind1] <- powerbase^bind2*bestC
      else if (loop==3) Cgrid[bind1] <- powerbase^(bind2/meadind)*bestC 
      bind2 <- bind2+1
    }
    for (bind1 in 1:Cgridinlev) {
      if (bind1==1) b <- meadind
      else if (bind1<=meadind) b <- bind1-1
           else b <- bind1 
      if (b==meadind && loop>1) {
        validres[meadind,loop] <- validres[bestb,loop]
        bestb <- meadind
      } else {
        C <- Cgrid[b]      
        validres[b,loop] <- algorithm(C=C,...)
        if (validres[b,loop] > bestres) {
          bestb <- b
          bestC <- Cgrid[b]
          bestres <- validres[b,loop]
        }
      }
    }
  }
  return(bestC)
}

