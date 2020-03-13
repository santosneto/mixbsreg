
require(VGAM)
require(rootSolve)

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


tobitResidualsNO <- function(model,nboot = 19,alpha=0.01,tau=0.1,intercept = "TRUE")
{

  n  <- summary(model)$n[1]
  y  <- as.numeric(model$y)[1:n]
  c <-  (1*(y>tau))
  muhat <- model$linear.predictors
  sigmahat <- model$scale
  deltahat <-(y-muhat)/sigmahat
  X <- model.matrix(model)
  var.explic <- X[,-1]

  S <- 1-pnorm(deltahat)
  rM <- c+log(S)


  rCS     <- -log(S)


  alpha1 <- ceiling(nboot*alpha)
  alpha2 <- ceiling(nboot*(1-alpha))
  e <- matrix(0,n,nboot)

  for(i in 1:nboot){

    ygerado  <- sigmahat*rnorm(n,0,1)+ muhat
    n1  <- summary(model)$n[2] #n?mero de obs. cens.
    pc       <- n1/n  #propor??o de obs. cens.
    tau1      <- sort(ygerado)[pc*n]
    yestrela <- ifelse(ygerado>tau1,ygerado,tau)
    status   <- (1*(ygerado>tau1))

    if(intercept == "FALSE"){form <- yestrela ~ var.explic - 1} else{form <- yestrela ~ var.explic}

    model1 <- tobit(form,left = tau)
    muhat1 <- model1$linear.predictors
    sigmahat1 <- model1$scale
    deltahat1 <- (yestrela - muhat1)/sigmahat1
    S1 <-1- pnorm(deltahat1)
    rCS1<- -log(S1) # Res?duo componente do desvio Martingal

    e[,i]    <- sort(rCS1)
  }
  e1<- numeric(n)
  e2<- numeric(n)

  for(j in 1:n){

    eo    <- sort(e[j,])
    e1[j] <- eo[alpha1]
    e2[j] <- eo[alpha2]
  }

  a  <-  qqplot(qexp(ppoints(500)),e1,plot.it=FALSE)$x
  a1 <-  qqplot(qexp(ppoints(500)),e1,plot.it=FALSE)$y
  b  <-  qqplot(qexp(ppoints(500)),e2,plot.it=FALSE)$x
  b1 <-  qqplot(qexp(ppoints(500)),e2,plot.it=FALSE)$y
  r  <-  qqplot(qexp(ppoints(500)),rCS,plot.it=FALSE)$x
  r1 <-  qqplot(qexp(ppoints(500)),rCS,plot.it=FALSE)$y

  xx <- c(a,rev(b))
  yy <- c(a1,rev(b1))
  med   <- apply(e,1,mean)
  faixa <- range(rCS,e1,e2,med)
  par(mar=c(4.0,4.0,0.1,0.1))
  plot(r,r1,type = "n", ylim=faixa,axes=FALSE,xlab="",ylab="")
  par(new=T)
  polygon(xx,yy,col="antiquewhite3",border=NA)
  par(new=T)
  qqplot(qexp(ppoints(500)),rCS,main="", ylim=faixa, ylab="R",xlab="Q", cex=0.7, pch=3)
  par(new=T)
  qqplot(qexp(ppoints(500)),e1,axes=F,type="l",main="",ylim=faixa,lty=1, xlab="", ylab="",col="gray")
  par(new=T)
  qqplot(qexp(ppoints(500)),e2,axes=F,type="l",main="",ylim=faixa,lty=1, xlab="", ylab="",col="gray")
  par(new=T)
  qqplot(qexp(ppoints(500)),med,axes=F,type="l",main="",ylim=faixa,lty=2, xlab="", ylab="",col="black")
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
#'

###################################################################
#                   Residuals - Birnbaum-Saunders
###################################################################

tobitResidualsBS <- function(y,alphahat,muhat,status){

  sfBS    <- 1-pnorm((2/alphahat)*sinh((y-muhat)/2))
  rCS     <- -log(sfBS)

  result  <- list( rCS = rCS
  )
  return(result)
}


tobitResidualsEnvelopeBS <- function(X,y,muhat,alphahat,n,k,status,nboot = 19,alpha=0.05)
{
  rcs1         <- tobitResidualsBS(y=y,alphahat=alphahat,muhat=muhat,status=status)$rCS
  alpha1       <- ceiling(nboot*alpha)
  alpha2       <- ceiling(nboot*(1-alpha))
  e            <- matrix(0,n,nboot)

  for(i in 1:nboot){
    ygerado     <-  muhat + 2*asinh((alphahat*rnorm(n))/2)
    n2          <-  initialValuesBS(X=X,y=y,status=status)$n2
    pc          <-  n2/n  #proporcao de obs. cens.
    tauenv      <-  sort(ygerado)[pc*n]
    yestrela    <-  ifelse(ygerado>tauenv,ygerado,0)
    statusenv   <-  (1*(ygerado>tauenv))
    psiStarenv  <-  initialValuesBS(X=X,y=yestrela,status=statusenv)$psiStar
    kenv        <-  initialValuesBS(X=X,y=yestrela,status=statusenv)$k
    model1      <-  tobitbs(X=X,y=yestrela,tau=tauenv,k=kenv,initialpoint=psiStarenv,status=statusenv,hessian="TRUE",method="BFGS",logtau="TRUE")
    muhat1      <- as.vector(X %*% model1$coef)
    alphahat1   <- model1$alphahat
    rCS11       <- tobitResidualsBS(y=yestrela,alphahat=alphahat1,muhat=muhat1,status=statusenv)$rCS
    e[,i]       <- sort(rCS11)
  }

  e1<- numeric(n)
  e2<- numeric(n)

  for(j in 1:n){

    eo    <- sort(e[j,])
    e1[j] <- eo[alpha1]
    e2[j] <- eo[alpha2]
  }

  a  <-  qqplot(qexp(ppoints(500)),e1,plot.it=FALSE)$x
  a1 <-  qqplot(qexp(ppoints(500)),e1,plot.it=FALSE)$y
  b  <-  qqplot(qexp(ppoints(500)),e2,plot.it=FALSE)$x
  b1 <-  qqplot(qexp(ppoints(500)),e2,plot.it=FALSE)$y
  r  <-  qqplot(qexp(ppoints(500)),rcs1,plot.it=FALSE)$x
  r1 <-  qqplot(qexp(ppoints(500)),rcs1,plot.it=FALSE)$y

  xx <- c(a,rev(b))
  yy <- c(a1,rev(b1))
  med   <- apply(e,1,mean)
  faixa <- range(rcs1,e1,e2,med)
  par(mar=c(4.0,4.0,0.1,0.1))
  plot(r,r1,type = "n", ylim=faixa,axes=FALSE,xlab="",ylab="")
  par(new=T)
  polygon(xx,yy,col="antiquewhite3",border=NA)
  par(new=T)
  qqplot(qexp(ppoints(500)),rcs1,main="", ylim=faixa, ylab="R",xlab="Q", cex=0.7, pch=3)
  par(new=T)
  qqplot(qexp(ppoints(500)),e1,axes=F,type="l",main="",ylim=faixa,lty=1, xlab="", ylab="",col="gray")
  par(new=T)
  qqplot(qexp(ppoints(500)),e2,axes=F,type="l",main="",ylim=faixa,lty=1, xlab="", ylab="",col="gray")
  par(new=T)
  qqplot(qexp(ppoints(500)),med,axes=F,type="l",main="",ylim=faixa,lty=2, xlab="", ylab="",col="black")
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

envelopet <- function(model,nboot = 19,alpha=0.01,tau=0.1,intercept = "TRUE")
{
  n <- summary(model)$n[1]
  y <- as.numeric(model$y)[1:n]
  nu <- model$parms
  muhat <- model$linear.predictors
  sigmahat <- model$scale
  deltahat <- (y - muhat)/sigmahat
  X <- model.matrix(model)
  var.explic <- X[, -1]
  nu <- model$parms
  S <- 1 - pt(deltahat, nu)
  rCS <- -log(S)


  alpha1 <- ceiling(nboot*alpha)
  alpha2 <- ceiling(nboot*(1-alpha))
  e <- matrix(0,n,nboot)

  for(i in 1:nboot){

    ygerado <- sigmahat * rt(n, nu) + muhat
    n1  <- summary(model)$n[2] #n?mero de obs. cens.
    pc       <- n1/n  #propor??o de obs. cens.
    tau1      <- sort(ygerado)[pc*n]
    yestrela <- ifelse(ygerado>tau1,ygerado,tau)
    status   <- (1*(ygerado>tau1))

    if(intercept == "FALSE"){form <- yestrela ~ var.explic - 1} else{form <- yestrela ~ var.explic}

    model1 <- tobit(form,dist="t",left=tau)
    muhat1 <- model1$linear.predictors
    sigmahat1 <- model1$scale
    deltahat1 <- (yestrela - muhat1)/sigmahat1
    nu1 <- model$parms
    S1 <- 1 - pt(deltahat1, nu1)
    rCS1<- -log(S1) # Res?duo componente do desvio Martingal

    e[,i]    <- sort(rCS1)
  }
  e1<- numeric(n)
  e2<- numeric(n)

  for(j in 1:n){

    eo    <- sort(e[j,])
    e1[j] <- eo[alpha1]
    e2[j] <- eo[alpha2]
  }

  a  <-  qqplot(qexp(ppoints(n)),e1,plot.it=FALSE)$x
  a1 <-  qqplot(qexp(ppoints(n)),e1,plot.it=FALSE)$y
  b  <-  qqplot(qexp(ppoints(n)),e2,plot.it=FALSE)$x
  b1 <-  qqplot(qexp(ppoints(n)),e2,plot.it=FALSE)$y
  r  <-  qqplot(qexp(ppoints(n)),rCS,plot.it=FALSE)$x
  r1 <-  qqplot(qexp(ppoints(n)),rCS,plot.it=FALSE)$y

  xx <- c(a,rev(b))
  yy <- c(a1,rev(b1))
  med   <- apply(e,1,mean)
  faixa <- range(rCS,e1,e2,med)
  par(mar=c(4.0,4.0,0.1,0.1))
  plot(r,r1,type = "n", ylim=faixa,axes=FALSE,xlab="",ylab="")
  par(new=T)
  polygon(xx,yy,col="antiquewhite3",border=NA)
  par(new=T)
  qqplot(qexp(ppoints(n)),rCS,main="", ylim=faixa, ylab="R",xlab="Q", cex=0.7, pch=3)
  par(new=T)
  qqplot(qexp(ppoints(n)),e1,axes=F,type="l",main="",ylim=faixa,lty=1, xlab="", ylab="",col="gray")
  par(new=T)
  qqplot(qexp(ppoints(n)),e2,axes=F,type="l",main="",ylim=faixa,lty=1, xlab="", ylab="",col="gray")
  par(new=T)
  qqplot(qexp(ppoints(n)),med,axes=F,type="l",main="",ylim=faixa,lty=2, xlab="", ylab="",col="black")
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

###################################################################
#                   Residuals - Birnbaum-Saunders/Bernoulli
###################################################################


CoxSnell_BSmixture <- function(x,X1,X2,theta,status,tau){

  k1 <- ncol(X1)
  k2 <- ncol(X2)

  theta1     <- theta[1:k1]
  theta2     <- theta[(k1+1):(k1+k2)]
  theta3     <- theta[((k1+k2)+1)]

    mu1        <- as.vector(X1 %*% theta1)
    mu2        <- as.vector(X2 %*% theta2)

  zeta1     <- (2 / theta3) * sinh((x - mu1) / 2)
  zeta1c    <- (2 / theta3) * sinh((log(tau) - mu1) / 2)
  pihat      <- 1 - (exp(mu2)/(1+exp(mu2)))

 CDF <- (pihat + (1 - pihat)*pnorm(zeta1c))*(1-status) +
        (pihat + (1 - pihat)*pnorm(zeta1c) +  (1 - pihat)* (pnorm(zeta1) - pnorm(zeta1c)))*status

 SF <- 1 - CDF

  return(-log(SF))
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


tobitResidualsEnvelopeBSBernoulli <- function(y,X1,X2,theta,status,tau,nboot = 19,alpha=0.05)
{
  n   <- length(y)
  rcs1         <- CoxSnell_BSmixture(x=y,X1=X1,X2=X2,theta=theta,status=status,tau=tau)
  alpha1       <- ceiling(nboot*alpha)
  alpha2       <- ceiling(nboot*(1-alpha))
  e            <- matrix(0,n,nboot)

  k1 <- ncol(X1)
  k2 <- ncol(X2)

  theta1     <- theta[1:k1]
  theta2     <- theta[(k1+1):(k1+k2)]
  theta3     <- theta[((k1+k2)+1)]

  mu1        <- as.vector(X1 %*% theta1)
  mu2        <- as.vector(X2 %*% theta2)
  pi         <- exp(mu2)/(1+exp(mu2))

  All <- cbind(y,X1,X2,status)



  for(i in 1:nboot){
     js <- sample(1:n, n,  replace = TRUE)

     Allger <- All[js,]
    #yger     <-  mu1 + 2*asinh((theta3*rnorm(n))/2)
    #ygerap   <- numeric()
    #for(j in 1:n)
    #{
    #ygerap[j]        <- sample(c(log(tau),yger[j]),size=1,replace=TRUE,prob=c((1-pi)[j],pi[j]))
    #}
    #ygerado  <- ygerap
    #status   <- ifelse(ygerado==log(tau),1,0)
     Allger
    psiStar    <- c(-0.910,  0.188,  0.074,  0.121, 0.650, 0.683,0.320,-0.274,1.165)
    est1      <-  Lbs(X1=Allger[,2:5],X2=Allger[,6:9],y=Allger[,1],status=Allger[,10],initialpoint = psiStar,tau=0.1,method="BFGS")


    coef1     <- est1$coef1
    coef2     <- est1$coef2
    alphahat  <- est1$alphahat
    thetahat = c(coef1,coef2,alphahat)

    rCS11       <-  CoxSnell_BSmixture(x=Allger[,1],X1=Allger[,2:5],X2=Allger[,6:9],theta=thetahat,status=Allger[,10],tau=tau)
    e[,i]       <- sort(rCS11)
  }

  e1<- numeric(n)
  e2<- numeric(n)

  for(j in 1:n){

    eo    <- sort(e[j,])
    e1[j] <- eo[alpha1]
    e2[j] <- eo[alpha2]
  }

  a  <-  qqplot(qexp(ppoints(500)),e1,plot.it=FALSE)$x
  a1 <-  qqplot(qexp(ppoints(500)),e1,plot.it=FALSE)$y
  b  <-  qqplot(qexp(ppoints(500)),e2,plot.it=FALSE)$x
  b1 <-  qqplot(qexp(ppoints(500)),e2,plot.it=FALSE)$y
  r  <-  qqplot(qexp(ppoints(500)),rcs1,plot.it=FALSE)$x
  r1 <-  qqplot(qexp(ppoints(500)),rcs1,plot.it=FALSE)$y

  xx <- c(a,rev(b))
  yy <- c(a1,rev(b1))
  med   <- apply(e,1,mean)
  faixa <- range(rcs1,e1,e2,med)
  par(mar=c(4.0,4.0,0.1,0.1))
  plot(r,r1,type = "n", ylim=faixa,axes=FALSE,xlab="",ylab="")
  par(new=T)
  polygon(xx,yy,col="antiquewhite3",border=NA)
  par(new=T)
  qqplot(qexp(ppoints(500)),rcs1,main="", ylim=faixa, ylab="R",xlab="Q", cex=0.7, pch=3)
  par(new=T)
  qqplot(qexp(ppoints(500)),e1,axes=F,type="l",main="",ylim=faixa,lty=1, xlab="", ylab="",col="gray")
  par(new=T)
  qqplot(qexp(ppoints(500)),e2,axes=F,type="l",main="",ylim=faixa,lty=1, xlab="", ylab="",col="gray")
  par(new=T)
  qqplot(qexp(ppoints(500)),med,axes=F,type="l",main="",ylim=faixa,lty=2, xlab="", ylab="",col="black")
}






