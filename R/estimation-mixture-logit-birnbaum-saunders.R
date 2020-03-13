#' Monte Carlo simulations
#'
#' @description In this function, we evaluate the performance of the mixture model through MC simulations.
#'
#' @usage Lbs(X,y,status,tau=0,initialpoint,method="BFGS",iterations=10000)
#'
#' @param X.
#' @param y shape parameter.
#' @param status Monte Carlo replications.
#' @param tau If is TRUE print the output.
#' @param initialpoint
#' @param method
#' @iterations
#'
#'
#'@export

descriptiveSummary <- function(x){
  n            <- length(x)
  skewness     <- (1 / n) * sum(((x - mean(x)) / sd(x)) ^ 3)
  kurtosis     <- ((1 / n) * sum(((x - mean(x)) / sd(x)) ^ 4) - 3 )
  cVariation   <- (sd(x) / mean(x))
  statistics   <- list(minimum              = round(min(x), 3),
                       median               = round(median(x), 3),
                       mean                 = round(mean(x), 3),
                       maximum              = round(max(x), 3),
                       standardDeviation    = round(sd(x), 3),
                       coefficientVariation = round(cVariation * 100, 3),
                       coefficientSkewness  = round(skewness, 3),
                       coefficientKurtosis  = round(kurtosis, 3)
                       #range                = round(max(x) - min(x), 3),
  )
  return(statistics)
}


#' Monte Carlo simulations
#'
#' @description In this function, we evaluate the performance of the mixture model through MC simulations.
#'
#' @usage Lbs(X,y,status,tau=0,initialpoint,method="BFGS",iterations=10000)
#'
#' @param X.
#' @param y shape parameter.
#' @param status Monte Carlo replications.
#' @param tau If is TRUE print the output.
#' @param initialpoint
#' @param method
#' @iterations
#'
#'
#'@export
#'
Lbs<- function(X1,X2,y,status,tau=0,initialpoint,method="BFGS",hessian="TRUE"){



    k1  <- ncol(X1)
    k2  <- ncol(X2)


  LogLik <- function(theta){

    theta1     <- theta[1:k1]
    theta2     <- theta[(k1+1):(k1+k2)]
    theta3     <- theta[((k1+k2)+1)]
    mu1        <- (X1 %*% theta1)
    mu2        <- (X2 %*% theta2)
    zetai1     <- (2 / theta3) * cosh((y - mu1) / 2) #Nao censurado
    zetai2     <- (2 / theta3) * sinh((y - mu1) / 2) #Nao censurado
    zetaic2    <- (2 / theta3) * sinh((log(tau) - mu1) / 2) #censurado

    result1 <- sum((1-status)*(log(1 + ( exp(mu2) * pnorm(zetaic2) ))  - log(1 + exp(mu2) ) ) + status * (-log(2) - (log(2 * pi) / 2) + (mu2) +
   log(zetai1) - ((1 / 2) * (zetai2 ^ 2))- log(1 + (exp(mu2)))   )    )


    return(-result1)
  }

  score <- function(theta){

    theta1   <- theta[1:k1]
    theta2   <- theta[(k1+1):(k1+k2)]
    theta3   <- theta[((k1+k2)+1)]
    mu1      <- (X1 %*% theta1)
    mu2      <- (X2 %*% theta2)
    zetai1   <- (2 / theta3) * cosh((y  - mu1) / 2) #Nao censurado
    zetai2   <- (2 / theta3) * sinh((y  - mu1) / 2) #Nao censurado
    zetaic1  <- (2 / theta3) * cosh((log(tau) - mu1) / 2) #censurado
    zetaic2  <- (2 / theta3) * sinh((log(tau) - mu1) / 2) #censurado


    Ualpha      <-  sum(((1-status)*(-1/theta3)*(exp(mu2)*dnorm(zetaic2)*zetaic2)/(1+(exp(mu2)*(pnorm(zetaic2)))))+((status)*(1/theta3)*((zetai2^2)-1)))

    Utheta1     <-  (1-status)*((-1/2)*((dnorm(zetaic2)*zetaic1*exp(mu2))/log(1+(exp(mu2)*(pnorm(zetaic2))))))
    + ((status)*(1/2)*((zetai1*zetai2)-(zetai2/zetai1)))


    Utheta2     <-  (1-status)*((exp(mu2)*(pnorm(zetaic2)))/(1+(exp(mu2)*((pnorm(zetaic2)))))  - (exp(mu2)/(1+exp(mu2))))
    + (status)*(1-(exp(mu2)/(1+exp(mu2))))

    result2     <- c(t(X1) %*% Utheta1,t(X2) %*% Utheta2, Ualpha)

    result2
  }

 ## est <- optim(initialpoint, LogLik,score,method = method, hessian = hessian)
  est <- optim(initialpoint, LogLik ,method = method, hessian = hessian)

  if(est$conv != 0)
    warning("FUNCTION DID NOT CONVERGE!")

  hessian             <- est$hessian
  I                   <- -solve(hessian)

  coef1               <- (est$par)[1:k1]
  coef2               <- (est$par)[(k1+1):(k1+k2)]
  alphahat            <- est$par[((k1+k2)+1)]

  stderrorsb1         <- sqrt(diag(abs(I)))[1:k1]
  stderrorsb2         <- sqrt(diag(abs(I)))[(k1+1):(k1+k2)]
  stderroralpha       <- sqrt(diag(abs(I)))[((k1+k2)+1)]

  zstats1             <- coef1 / stderrorsb1
  pvalues1            <- 2 * (1 - pnorm(abs(coef1 / stderrorsb1)))

  zstats2             <- coef2 / stderrorsb2
  pvalues2            <- 2 * (1 - pnorm(abs(coef2 / stderrorsb2)))

  pvaluea             <- 2 * (1 - pnorm(abs(alphahat / stderroralpha)))
  conv  <- est$conv



  result3 <- list( #est      = est,
                   #value    = est$value,
                   #conv     = conv,
                   coef1    = coef1,
                   coef2    = coef2,
                   alphahat = alphahat,
                   #I        = I,
                   stderrorsb1 = stderrorsb1,
                   stderrorsb2 = stderrorsb2,
                   stderroralpha = stderroralpha,
                   #zstats1 = zstats1,
                   #zstats2 = zstats2,
                   pvalues1 = pvalues1,
                   pvalues2 = pvalues2,
                   pvaluea = pvaluea
  )
  return(result3)
}


#' Monte Carlo simulations
#'
#' @description In this function, we evaluate the performance of the mixture model through MC simulations.
#'
#' @usage Lbs(X,y,status,tau=0,initialpoint,method="BFGS",iterations=10000)
#'
#' @param X.
#' @param y shape parameter.
#' @param status Monte Carlo replications.
#' @param tau If is TRUE print the output.
#' @param initialpoint
#' @param method
#' @iterations
#'
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

#' Monte Carlo simulations
#'
#' @description In this function, we evaluate the performance of the mixture model through MC simulations.
#'
#' @usage Lbs(X,y,status,tau=0,initialpoint,method="BFGS",iterations=10000)
#'
#' @param X.
#' @param y shape parameter.
#' @param status Monte Carlo replications.
#' @param tau If is TRUE print the output.
#' @param initialpoint
#' @param method
#' @iterations
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
