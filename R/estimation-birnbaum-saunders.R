#' Fitting Linear Models - tobit bs
#'
#' @description Lbs is used to fit mixture model.
#'
#' @usage tobitbs(X,y,status,tau=0,k,initialpoint,method="BFGS",iterations=10000,hessian="TRUE",logtau="FALSE")
#'
#' @param X The model matrix.
#' @param y The response used.
#' @param status The status indicator, normally 0=alive, 1=dead.
#' @param tau Is a prefixed limiting values. Default for ZERO.
#' @param k value.
#' @param initialpoint Initial values for the parameters to be optimized over.
#' @param method The method to be used. The default is BFGS.
#' @param iterations The maximum number of iterations. Default to 10000.
#' @param hessian Logical. Should a numerically differentiated Hessian matrix be returned?
#' @param logtau Logical. Default for FALSE.
#'
#'
#'@export

tobitbs  <- function(X,y,status,tau=0,k,initialpoint,method="BFGS",iterations=10000,hessian="TRUE",logtau="FALSE"){

  #if(xi[1]<=0) stop("the extra parameter must be positive!!",call.=FALSE)
  #alpha <- xi[1]
 # v <- function(z) 4*sinh(z)*cosh(z)/(alpha^2*z) - tanh(z)/z
#  vp <- function(z) ((cosh(z)*z-sinh(z))/z^2)*(4*cosh(z)/alpha^2 - 1/cosh(z)) + (sinh(z)/z)*(4*sinh(z)/alpha^2 + sinh(z)/(cosh(z)^2))
#  dg <- 2 + 4/(alpha^2) - (sqrt(2*pi)/alpha)*(1-2*(pnorm(sqrt(2)/alpha,mean=0,sd=sqrt(2)/2)-0.5))*exp(2/(alpha^2))
#  dshn <- function(z) 2*cosh(z)*exp(-(2/alpha^2)*(sinh(z))^2)/(alpha*sqrt(2*pi))
#  fgf <- function(z) dshn(z)*(4*sinh(z)*cosh(z)/(alpha^2) - tanh(z))^2*z^2
#  fg <- 2*integrate(fgf,0,60)$value
#  deviance.mu <- function(z){if(alpha<=2) abs(4*(sinh(z))^2/alpha^2 - log((cosh(z))^2))
#    else{2*log(dshn(acosh(alpha/2))/dshn(z))}}
#  tau <- uniroot(function(x) 4*x*sinh(x)*cosh(x)/(alpha^2) - tanh(x)*x -1, lower=0, upper=50)$root
#  deviance.phi <- function(z){
#    a <- 2*log(cosh(tau)*exp(-(2/alpha^2)*(sinh(tau))^2)) + log(tau^2)
#    b <- 2*log(cosh(z)*exp(-(2/alpha^2)*(sinh(z))^2)) + log(z^2)
#    ifelse(a<b,0,a-b)}
  cdf <- function(z,alpha) pnorm(2*sinh(z)/alpha)
  pdf <- function(z,alpha) (2*cosh(sqrt(z^2))*exp(-2*sinh(sqrt(z^2))*sinh(sqrt(z^2))/alpha^2)/(sqrt(2*pi)*alpha))
#  fgf <- function(z) dshn(z)*z^2
#  xix <- 2*integrate(fgf,0,20)$value



  loglikbs <- function(theta){

    theta1     <- theta[1:k]     ##regressors
    theta2     <- theta[k + 1]
    alpha      <- theta2
    phi_es     <- 2
    mu_es      <- (X %*% theta1)


    if(logtau=="TRUE"){tau1 <- tau} else{
      tau1 <- log(tau)
    }

    z_cen      <- (tau1 - mu_es)/phi_es
    z_unc      <- (y -    mu_es)/phi_es

   result1 <- sum((1 - status)*log(cdf(z=z_cen,alpha)) +   status*(log(pdf(z=z_unc,alpha)) - log(phi_es)))

   return(-result1)
  }

  score <- function(theta){

    theta1    <-theta[1:k]
    theta2    <-theta[k+1]
    mu        <-(X %*% theta1)
    tau1      <-  log(tau)
    zetai1    <- (2 / theta2) * cosh((y - mu) / 2) ##N?o censurado
    zetai2    <- (2 / theta2) * sinh((y - mu) / 2) ##N?o censurado
    delta1    <- ((y-mu)/2)
    deltac1   <- ((tau1-mu)/2)
    const     <- 1/sqrt(2*pi)
    zetaic1   <- (2 / theta2) * cosh((tau1 - mu) / 2) ##N?o censurado
    zetaic2   <- (2 / theta2) * sinh((tau1 - mu) / 2) ##N?o censurado
    lambda    <-(dnorm(zetaic2)/pnorm(zetaic2))

    Ualpha     <-  sum(status*((zetai2-1)/theta2) -
                         (1-status)*((lambda*zetaic2)/theta2))



    Utheta <-  status*((sinh(2*delta1)/(theta2^2)) - (tanh(delta1)/2)) - (1-status)*
      ((cosh(deltac1)*lambda)/theta2)

    result2     <- c(t(X) %*% Utheta, Ualpha)
    return(result2)
  }

  est <- optim(initialpoint, loglikbs, score, method = method, hessian = hessian,
               control = list(fnscale = 1, maxit = iterations, reltol = 1e-30))

  conv  <- est$conv

  if(est$conv != 0)
    warning("FUNCTION DID NOT CONVERGE!")

  SHess   = solve(-est$hessian)
  SE      = sqrt(diag(SHess))
  se.coef = SE
  tval    = c((est$par)[1:k],est$par[k+1])/se.coef
  matcoef = cbind(c((est$par)[1:k],est$par[k+1]), se.coef, tval, 2*(1-pnorm(abs(tval))))

  AIC    = 2*est$value + 2*(k+1)
  BIC    = 2*est$value + (k+1)*log(length(y))

  result3 <- list( #est      = est,
                  coef     = (est$par)[1:k],
                  alphahat = est$par[k+1],
                  value    = (-est$value),
                  conv     = conv,
                  matcoef  = matcoef,
                  ICVs     = c(AIC, BIC)
  )

  return(result3)
  }





