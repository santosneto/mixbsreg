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
tobitLSVectorY <- function(t,tau){
  if(any(t<=0)){
stop("Logarithmic functions are not defined for values less than or equal to zero")
  }

  else{

      if(any(tau!=min(t)))
        {
  stop("The threshold must be the minimum value of the selected variable.")
       }
       {
      y     <- ifelse(t==tau,log(tau),log(t))
    }
    }
  return(y)
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

tobitLSVectorStatus <- function(t,tau){
  if(any(t<=0)){
stop("Logarithmic functions are not defined for values less than or equal to zero.")
  }

  else{

    if(any(tau!=min(t)))
    {
      stop("The threshold must be the minimum value of the selected variable.")
    }
    {
      status <- ifelse(t==tau,0,1)
    }
    }
  return(status)
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

initialValuesLS <- function(y,X,status){

  lm                 <- lm.fit(X,y)
  status             <- as.matrix(status) # 0 para censura e 1 sem censura
  coef               <- lm$coef
  k                  <- length(coef) ### Numero de parametros
  n                  <- length(y)
  n1                 <- sum(status) ###Quantidade de observacoes nao censuradas
  n2                 <- n-n1        ###Quantidade de observavoes censuradas
  mu                 <- X %*% coef
  phi                <- ((t(y-mu)%*%(y-mu))/n)^(1/2)
  psiStar            <- c(coef,phi) ###Ponto inicial que sera usado na otimizacao

  result    <-list( coef              = coef,
                     phi              = phi,
                     psiStar           = psiStar,
                     n                 = n,
                     n1                = n1,
                     n2                = n2,
                     k                 = k
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

initialValuesBS <- function(y,X,status){

  lm                 <- lm.fit(X,y)
  status             <- as.matrix(status) # 0 para censura e 1 sem censura
  coef               <- lm$coef
  k                  <- length(coef) ### Numero de parametros
  n                  <- length(y)
  n1                 <- sum(status) ###Quantidade de observacoes nao censuradas
  n2                 <- n-n1        ###Quantidade de observavoes censuradas
  mu                 <- X %*% coef
  alpha              <-  sqrt(4/n * sum ((sinh((y-mu)/2))^2))
  psiStar            <- c(coef,alpha) ###Ponto inicial que sera usado na otimizacao

  result    <-list( coef              = coef,
                    alpha             = alpha,
                    psiStar           = psiStar,
                    n                 = n,
                    n1                = n1,
                    n2                = n2,
                    k                 = k
  )

  return(result)
}

