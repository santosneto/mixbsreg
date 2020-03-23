#' Fitting Linear Models
#'
#' @description Lbs is used to fit mixture model.
#'
#' @usage Lbs0(X,y,status,tau=0,initialpoint,method="BFGS",iterations=10000)
#' @usage Lbs(X1,X2,y,status,tau=0,initialpoint,method="BFGS",hessian="TRUE")
#'
#' @param X1 The model matrix.
#' @param X2 The model matrix.
#' @param X The model matrix.
#' @param y The response used.
#' @param status The status indicator, normally 0=alive, 1=dead.
#' @param tau Is a prefixed limiting values. Default for ZERO.
#' @param initialpoint Initial values for the parameters to be optimized over.
#' @param method The method to be used. The default is BFGS.
#' @param iterations The maximum number of iterations. Default to 10000.
#' @param hessian Logical. Should a numerically differentiated Hessian matrix be returned?
#'
#'
#'@export


Lbs <- function(X1,X2,y,status,tau=0,initialpoint,method="BFGS",hessian="TRUE"){
  k1  <- ncol(X1)
  k2  <- ncol(X2)
   n  <- length(y)

  LogLik <- function(theta){

    theta1     <- theta[1:k1]
    theta2     <- theta[(k1+1):(k1+k2)]
    theta3     <- theta[((k1+k2)+1)]

    mu1        <- (X1 %*% theta1)
    mu2        <- (X2 %*% theta2)

    zetai1     <- (2 / theta3) * cosh((y - mu1) / 2) #Nao censurado
    zetai2     <- (2 / theta3) * sinh((y - mu1) / 2) #Nao censurado
    zetaic2    <- (2 / theta3) * sinh((log(tau) - mu1) / 2) #censurado

    result1 <- sum((1-status)*(log(1 + ( exp(mu2) * pnorm(zetaic2) ))  - log(1 + exp(mu2) ) ) + status * (-log(2) - (log(2*pi)/2) + (mu2) + log(zetai1) - ((1 / 2) * (zetai2 ^ 2))- log(1 + (exp(mu2)))   )    )


    return(result1)
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

    return(result2)
  }

  ## est <- optim(initialpoint, LogLik,score,method = method, hessian = hessian)
  est <- optim(initialpoint, fn=LogLik ,gr=score,method = method, hessian = hessian,control=list(fnscale=-1))

  if(est$conv != 0)
    warning("FUNCTION DID NOT CONVERGE!")

  hessian             <- -as.matrix(est$hessian)
  I                   <- solve(hessian)

  coef1               <- (est$par)[1:k1]
  p_cont                   <- length(coef1)
  coef2               <- (est$par)[(k1+1):(k1+k2)]
  p_log                  <- length(coef2)
  alphahat            <- est$par[((k1+k2)+1)]
  p_disp              <- length(alphahat)
  stderrorsb1         <- sqrt(diag(I))[1:k1]
  stderrorsb2         <- sqrt(diag(I))[(k1+1):(k1+k2)]
  stderroralpha       <- sqrt(diag(I))[((k1+k2)+1)]
  df                  <- n - p_log - p_cont - p_disp
  zstats1             <- coef1/stderrorsb1
  pvalues1            <- 2*pt(abs(zstats1), df, lower.tail = FALSE)

  zstats2             <- coef2/stderrorsb2
  pvalues2            <- 2*pt(abs(zstats2 ), df, lower.tail = FALSE)

  zstatsa             <- alphahat/stderroralpha
  pvaluea             <- 2*pt(abs(zstatsa), df, lower.tail = FALSE)

  conv  <- est$conv
  loglink_out <-  est$value


  result3 <- list( #est      = est,
    #value    = est$value,
    x1     = X1,
    x2     = X2,
    status = status,
    y =y,
    coef1    = coef1,
    coef2    = coef2,
    alphahat = alphahat,
    vcov        = I,
    stderrorsb1 = stderrorsb1,
    stderrorsb2 = stderrorsb2,
    stderroralpha = stderroralpha,
    #zstats1 = zstats1,
    #zstats2 = zstats2,
    pvalues1 = pvalues1,
    pvalues2 = pvalues2,
    pvaluea = pvaluea,
    loglink = loglink_out
  )
  return(result3)
}


#' Fitting Linear Models
#'
#' @description Lbs is used to fit mixture model.
#'
#' @usage Lbs0(X,y,status,tau=0,initialpoint,method="BFGS",iterations=10000)
#' @usage Lbs(X1,X2,y,status,tau=0,initialpoint,method="BFGS",hessian="TRUE")
#'
#' @param X1 The model matrix.
#' @param X2 The model matrix.
#' @param X  The model matrix.
#' @param y The response used.
#' @param status The status indicator, normally 0=alive, 1=dead.
#' @param tau Is a prefixed limiting values. Default for ZERO.
#' @param initialpoint Initial values for the parameters to be optimized over.
#' @param method The method to be used. The default is BFGS.
#' @param iterations The maximum number of iterations. Default to 10000.
#' @param hessian Logical. Should a numerically differentiated Hessian matrix be returned?
#'
#'
#'@export

Lbs0 <- function(X,y,status,tau=0,initialpoint,method="BFGS",
               iterations=10000){

    X1  <- X
    X2  <- X

    k1  <- dim(X1)[2]
    k2  <- dim(X2)[2]


  loglik <- function(theta){

    theta1     <- theta[1:k1]
    theta2     <- theta[(k1+1):(k1+k2)]
    theta3     <- theta[((k1+k2)+1)]
    mu1        <- (X1 %*% theta1)
    mu2        <- (X2 %*% theta2)
    zetai1     <- (2 / theta3) * cosh((y - mu1) / 2) #non censured
    zetai2     <- (2 / theta3) * sinh((y - mu1) / 2) #non censured
    zetaic2    <- (2 / theta3) * sinh((log(tau) - mu1) / 2) #censured


    result1 <- sum((status*(log(1 + (exp(mu2) * (abs(pnorm(zetaic2) - 1))))- log(1 + (exp(mu2))))) + (1 - status) * (-log(2) - (log(2 * pi) / 2) + (mu2) + log(zetai1) - ((1 / 2) * (zetai2 ^ 2))- log(1 + (exp(mu2)))))


    return(-result1)
  }

  score <- function(theta){

    theta1   <- theta[1:k1]
    theta2   <- theta[(k1+1):(k1+k2)]
    theta3   <- theta[((k1+k2)+1)]
    mu1      <- (X1 %*% theta1)
    mu2      <- (X2 %*% theta2)
    zetai1   <- (2 / theta3) * cosh((y  - mu1) / 2) #non censured
    zetai2   <- (2 / theta3) * sinh((y  - mu1) / 2) #non censured
    zetaic1  <- (2 / theta3) * cosh((log(tau) - mu1) / 2) #censured
    zetaic2  <- (2 / theta3) * sinh((log(tau) - mu1) / 2) #censured


    Ualpha      <-  sum((status*(-1/theta3)*(exp(mu2)*dnorm(zetaic2)*zetaic2)/(1+(exp(mu2)*(pnorm(zetaic2)-1))))+
                          ((1-status)*(1/theta3)*((zetai2^2)-1)))

    Utheta1     <-  status*((-1/2)*((dnorm(zetaic2)*zetaic1*exp(mu2))/log(1+(exp(mu2)*abs(pnorm(zetaic2)-1)))))
    + ((1-status)*(1/2)*((zetai1*zetai2)-(zetai2/zetai1)))


    Utheta2     <-  status*((exp(mu2)*(pnorm(zetaic2)-1))/(1+(exp(mu2)*abs((pnorm(zetaic2)-1))))  - (exp(mu2)/(1+exp(mu2))))
    + (1-status)*(1-(exp(mu2)/(1+exp(mu2))))

    result2     <- c(t(X1) %*% Utheta1,t(X2) %*% Utheta2, Ualpha)
    return(result2)
  }

  est <- optim(initialpoint, loglik, score ,method = method, hessian = FALSE,
               control = list(fnscale = 1, maxit = iterations, reltol = 1e-30))

  coef1               <- (est$par)[1:k1]
  coef2               <- (est$par)[(k1+1):(k1+k2)]
  alphahat            <- est$par[((k1+k2)+1)]

  conv  <- est$conv

  if(est$conv != 0)
    warning("FUNCTION DID NOT CONVERGE!")

  result3 <- list( est      = est,
                   value    = est$value,
                   conv     = conv,
                   coef1    = coef1,
                   coef2    = coef2,
                   alphahat = alphahat
  )
  return(result3)
}


#' Monte Carlo simulations
#'
#' @description In this function, we evaluate the performance of the mixture model through MC simulations.
#'
#' @usage sim_mix_bs(n=100,alpha=0.1,nrep=5000,print=FALSE)
#'
#' @param n sample size.
#' @param alpha shape parameter.
#' @param nrep Monte Carlo replications.
#' @param print If is TRUE print the output.
#'
#'@importFrom VGAM rbisa
#'
#'@export

sim_mix_bs <- function(n=100,beta1=c(0.2,0.5),beta2=c(1,2),alpha=0.1,nrep=5000,print=FALSE)
{
#MC simulation logit/BS model
NREP <- nrep
n    <- n
set.seed(c(2000,1900),kind="Marsaglia-Multicarry")
x <- runif(n)
X <- cbind(1, x)

#Setup for parameters
beta1 <- beta1
beta2 <- beta2
alpha <- alpha

beta10 <- beta1[1]
beta11 <- beta1[2]
beta20 <- beta2[1]
beta21 <- beta2[2]

mu1   <- as.vector(X%*%beta1)
mu2   <- as.vector(X%*%beta2)


k1    <- length(beta1)
k2    <- length(beta2)

#auxiliar matrixes and vectors
MCalpha  <- numeric(NREP)
MCbeta20 <- numeric(NREP)
MCbeta21 <- numeric(NREP)
MCbeta10 <- numeric(NREP)
MCbeta11 <- numeric(NREP)
t        <- numeric(n)

p        <- exp(mu2)/(1+exp(mu2))
tau      <- 1

cont  <- 0
contv <- 0

while(cont < NREP) {

    cont <- cont + 1

  #  prop <- (cont / NREP) * 100
   # if(prop==0 || prop==25 || prop==50 ||prop==75|| prop==100)
  #  cat(paste(prop,"%"),"\n")

    x <- runif(n)
    X <- cbind(1, x)

    mu1   <- as.vector(X%*%beta1)
    mu2   <- as.vector(X%*%beta2)

    for(j in 1:n)
    {
    t[j]     <- sample(c(1,exp(mu1[j]) * rbisa(n=1,alpha,mu1[j])),size=1,replace=TRUE,prob=c((1-p)[j],p[j]))
    }

  ygerado  <- log(t)
  status   <- ifelse(ygerado==log(tau),1,0)

  b10 <- rnorm(1,mean=beta10,sd=0.095)
  b11 <- rnorm(1,mean=beta11,sd=0.095)
  b20 <- rnorm(1,mean=beta20,sd=0.095)
  b21 <- rnorm(1,mean=beta21,sd=0.095)
  alp <- rnorm(1,mean=alpha,sd=0.018)

  start  <- c(b10,b11,b20,b21,alp) # initial values

  opt    <- suppressWarnings(Lbs0(X=X,y=ygerado,tau=tau,initialpoint = start,status=status))

  if(opt$conv!=0){

    cont  <- cont - 1
    contv <- contv + 1
  }
     else{

fitc1  <- opt$coef1
fitc2  <- opt$coef2
fital  <- opt$alphahat

MCbeta10[cont] <- fitc1[1]
MCbeta11[cont] <- fitc1[2]
MCbeta20[cont] <- fitc2[1]
MCbeta21[cont] <- fitc2[2]
MCalpha[cont]  <- fital

}
}
#MC estimators
MCbeta20hat <- mean(MCbeta20)
MCbeta21hat <- mean(MCbeta21)
MCbeta10hat <- mean(MCbeta10)
MCbeta11hat <- mean(MCbeta11)
MCalphahat  <- mean(MCalpha)

out_mean <- rbind(MCalphahat,MCbeta10hat,MCbeta11hat,MCbeta20hat,MCbeta21hat)
#bias
bias20       <- MCbeta20hat - beta20
bias21       <- MCbeta21hat - beta21
bias10       <- MCbeta10hat - beta10
bias11       <- MCbeta11hat - beta11
biasalpha    <- MCalphahat  - alpha

out_bias <- rbind(biasalpha,bias10,bias11,bias20,bias21)

#MSE

MSE20     <- mean((MCbeta20 - beta20) ^ 2)
MSE21     <- mean((MCbeta21 - beta21) ^ 2)
MSE10     <- mean((MCbeta10 - beta10) ^ 2)
MSE11     <- mean((MCbeta11 - beta11) ^ 2)
MSEalpha  <- mean((MCalpha  - alpha)  ^ 2)


out_mse <- rbind(MSEalpha,MSE10,MSE11,MSE20,MSE21)

out_end <- cbind(out_mean,out_bias,out_mse)
rownames(out_end) <- c('alpha','beta(1)0','beta(1)1','beta(2)0','beta(2)1')
colnames(out_end) <- c('Mean','Bias','MSE')

if(print == TRUE)
{
cat(sprintf('Proportion of censored data: %0.2f \n',sum(status)/length(status)))
cat(sprintf('MC Replications: %i \n',NREP))
cat(sprintf('Sample size: %i \n',n))
cat(sprintf('Shape: %0.1f \n',alpha))
cat(sprintf('Estimates: \n'))
print(round(out_end,5))
}else return(round(out_end,5))
}



