#' Fitting Linear Models - CI
#'
#'@export
LbsIc     <- function(model){
  k1  <- dim(model$x1)[2]
  k2  <- dim(model$x2)[2]
  n   <- dim(model$x1)[1]

  loglikIC <- model$loglink


  AIC     <- - 2 * loglikIC + 2 * (k1 + k2 + 1)


  BIC     <- - 2 * loglikIC + (k1 + k2 + 1)  * log(n)

  result  <- list ( AIC = AIC,
                    BIC = BIC
  )
  return(result)
}

#' Fitting Linear Models - Profile
#'
#'
#'@export

LbsProfile     <- function(X1,X2,y,coef1,coef2,from,to,by,tau,
                           status=status, ylim=c(0,2000),xlim=c(0,100),plot=TRUE){

  X1  <- model.matrix(~ X1)
  X2  <- model.matrix(~ X2)

  loglik <- function(alphahat){

    mu1        <- (X1 %*% coef1)
    mu2        <- (X2 %*% coef2)
    zetai1     <- (2 / alphahat) * cosh((y - mu1) / 2) #Nao censurado
    zetai2     <- (2 / alphahat) * sinh((y - mu1) / 2) #Nao censurado
    zetaic2    <- (2 / alphahat) * sinh((log(tau) - mu1) / 2) #censurado
    #lambda    <- (dnorm(zetaic2)/pnorm(zetaic2))

    result <- sum((1-status)*((exp(mu2)*(pnorm(zetaic2)-1)-log(1+exp(mu2))) )+
                    ((status)*((mu2)+(-log(2)-log(2*pi/2))+log(abs(zetai1))
                                 -((zetai2^2)/2)-log(1+exp(mu2)))))

    return(-result)
  }
  seq         <-  seq(from=from,to=to,by=by)
  vect        <-  Vectorize(loglik)(seq)

  if(plot==TRUE){
    plot(seq,vect,type="l",main="" ,xlab="r",ylab="y",lwd=2,
         xlim=xlim,ylim=ylim)
  }
  else{
    result <- list(vect = vect,
                   seq  = seq)
    return(result)
  }
}
###################################################################
