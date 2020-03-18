#' Descriptive Measures
#'
#' @description Central and dispersion measures.
#'
#' @usage descriptiveSummary(x)
#'
#' @param x data.
#'
#'@importFrom moments kurtosis skewness
#'
#'@export

descriptiveSummary <- function(x){
  n            <- length(x)
  cVariation   <- (sd(x) / mean(x))
  statistics   <- list(minimum              = min(x),
                       median               = median(x),
                       mean                 = mean(x),
                       maximum              = max(x),
                       standardDeviation    = sd(x),
                       coefficientVariation = cVariation,
                       coefficientSkewness  = skewness(x),
                       coefficientKurtosis  = kurtosis(x)
                       #range                = round(max(x) - min(x), 3),
  )
  return(statistics)
}


#' Response tobit
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


#' Status Vector
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

#' Initial Values
#'
#' @description Initial values for the parameters to be optimizes over.
#'
#' @usage initialValuesLS(y,X,status)
#'
#' @param y The response variable.
#' @param X The model matrix.
#' @param status The status indicator, normally 0=alive, 1=dead.
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

#' Initial Values for Birnbaum-Saunders model
#'
#' @description Initial values for the parameters to be optimizes over.
#'
#' @usage initialValuesLS(y,X,status)
#'
#' @param y The response variable.
#' @param X The model matrix.
#' @param status The status indicator, normally 0=alive, 1=dead.
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

