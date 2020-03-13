#' Fitting Linear Models - CI
#'
#'@export
LbsIc     <- function(X1,X2,y,tau=0.1,status,coef1, coef2, alphahat){

  X1  <- model.matrix(~ X1)
  X2  <- model.matrix(~ X2)

  k1  <- dim(X1)[2]
  k2  <- dim(X2)[2]
  n   <- dim(X1)[1]

  loglikIC <- function(theta){

    theta1     <- theta[1:k1]
    theta2     <- theta[(k1+1):(k1+k2)]
    theta3     <- theta[((k1+k2)+1)]
    mu1        <- (X1 %*% theta1)
    mu2        <- (X2 %*% theta2)
    zetai1     <- (2 / theta3) * cosh((y - mu1) / 2) #Nao censurado
    zetai2     <- (2 / theta3) * sinh((y - mu1) / 2) #Nao censurado
    zetaic2    <- (2 / theta3) * sinh((log(tau) - mu1) / 2) #censurado


  #  result <- sum(((1-status)*(log(1 + (exp(mu2) * (abs(pnorm(zetaic2) - 1))))-
  #                               log(1 + (exp(mu2))))) + (status) * (-log(2) - (log(2 * pi) / 2) + (mu2) + log(zetai1) - ((1 / 2) * (zetai2 ^ 2))- log(1 + (exp(mu2)))))

    result <- sum((1-status)*(log(1 + ( exp(mu2) * pnorm(zetaic2) ))  - log(1 + exp(mu2) ) ) + status * (-log(2) - (log(2 * pi) / 2) + (mu2) +
                                                                                                            log(zetai1) - ((1 / 2) * (zetai2 ^ 2))- log(1 + (exp(mu2)))   )    )

    return(result)
  }


  AIC     <- - 2 * loglikIC(c(coef1,coef2,alphahat)) + 2 * (k1 + k2 + 1)


  BIC     <- - 2 * loglikIC(c(coef1,coef2,alphahat)) + (k1 + k2 + 1)  * log(n)

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
