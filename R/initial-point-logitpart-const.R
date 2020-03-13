#' Fitting Linear Models - Constat
#'
#'@export

constLogitpart     <- function(X1,X2,y,coef1,alphahat,from,to,by,tau,status=status, ylim=c(0,2000),xlim=c(0,100),plot=TRUE){

  X1  <- model.matrix(~ X1)
  X2  <- model.matrix(~ X2)
  coef2 <- initialpoint

  loglik <- function(coef2){

    mu1        <- (X1 %*% coef1)
    mu2        <- (X2 %*% coef2)
    zetai1     <- (2 / alphahat) * cosh((y - mu1) / 2) #Nao censurado
    zetai2     <- (2 / alphahat) * sinh((y - mu1) / 2) #Nao censurado
    zetaic2    <- (2 / alphahat) * sinh((tau - mu1) / 2) #censurado
    #lambda    <- (dnorm(zetaic2)/pnorm(zetaic2))

    result <- sum((  (1-status)*(log(1 + (exp(mu2) * (abs(pnorm(zetaic2) - 1))))- log(1 + (exp(mu2))))) + (status) * (-log(2) - (log(2 * pi) / 2) + (mu2) + log(zetai1) - ((1 / 2) * (zetai2 ^ 2))- log(1 + (exp(mu2)))))


    return(-result)
  }
  seq         <-  seq(from=from,to=to,by=by)
  vect        <-  Vectorize(loglik)(seq)

  if(plot==TRUE){
    plot(seq,vect,type="l",main="" ,xlab="r",ylab="y",lwd=2,
         xlim=xlim,ylim=ylim)
  }
  else{
    result <- list(vect = vect,seq  = seq)
    return(result)
  }
}
