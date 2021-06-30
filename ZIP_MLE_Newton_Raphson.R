### zip.mle : compute MLE for ZIP using newton raphson method
# install.packages("Rfast")
library(Rfast)

##### zip.mle.NR : (my function) compute MLE for ZIP using newton raphson method
zip.mle.NR <- function (x, tol = 1e-09){
  n0 <- sum(x == 0)
  n <- length(x)
  prob <- n0/n
  n1 <- n - n0
  x1 <- x[x > 0]
  lambda <- mean(x)
  dlp <- (n0*(1-exp(-lambda))/(prob +(1-prob)*exp(-lambda))) - n1/(1-prob)
  dll <- (-n0*(1-prob)*exp(-lambda))/(prob +(1-prob)*exp(-lambda)) - n1 + sum(x1)/lambda
  dlpp <- (-n0*((1-exp(-lambda))^2))/((prob +(1-prob)*exp(-lambda))^2) - n1/((1-prob)^2)
  dlll <- (n0*prob*(1-prob)*exp(-lambda))/((prob +(1-prob)*exp(-lambda))^2) - sum(x1)/(lambda^2)
  dlpl <- (n0*exp(-lambda))/((prob+(1-prob)*exp(-lambda))^2)
  S <- c(dll, dlp) #score function
  H <-matrix(c(dlll, dlpl, dlpl, dlpp), 2,2,byrow=T)   #Hessian matrix
  Theta <- c(lambda, prob)
  Theta2 <- Theta - solve(H)%*%S
  i <- 0
  while ((abs(Theta2[1] - Theta[1])+abs(Theta2[2] - Theta[2])) > tol) {
    i <- i + 1
    Theta <- Theta2
    dlp <- (n0*(1-exp(-Theta[1]))/(Theta[2] +(1-Theta[2])*exp(-Theta[1]))) - n1/(1-Theta[2])
    dll <- (-n0*(1-Theta[2])*exp(-Theta[1]))/(Theta[2] +(1-Theta[2])*exp(-Theta[1])) - n1 + sum(x1)/Theta[1]
    dlpp <- (-n0*((1-exp(-Theta[1]))^2))/((Theta[2] +(1-Theta[2])*exp(-Theta[1]))^2) - n1/((1-Theta[2])^2)
    dlll <- (n0*Theta[2]*(1-Theta[2])*exp(-Theta[1]))/((Theta[2] +(1-Theta[2])*exp(-Theta[1]))^2) - sum(x1)/(Theta[1]^2)
    dlpl <- (n0*exp(-Theta[1]))/((Theta[2]+(1-Theta[2])*exp(-Theta[1]))^2)
    S <- c(dll, dlp) #score function
    H <-matrix(c(dlll, dlpl, dlpl, dlpp), 2,2,byrow=T)   #Hessian matrix
    Theta2 <- Theta - solve(H)%*%S
  }
  p = Theta[2]; lambda = Theta[1]
  loglik <- n0*log(p + (1-p)*exp(-lambda)) + n1*log(1-p) -n1*lambda + sum(x1)*log(lambda) - sum(log(factorial(x1)))
  param <- t(Theta)
  colnames(param) <- c("lambda", "pi")
  list(iters = i, loglik = loglik, param = param)
}



x <- rpois(1000, 10)
num = sample(1:1000,200)
x[num] <- 0
zip.mle(x)
zip.mle.NR(x)




### zip.mle code
zip.mle <- function (x, tol = 1e-09) {
  no <- sum(x == 0)
  n <- length(x)
  prop <- no/n
  n1 <- n - no
  x1 <- x[x > 0]
  sx <- sum(x1)
  m <- sx/n
  s <- (sum(x1^2) - m * sx)/(n - 1)
  l1 <- s/m + m - 1
  fx <- m - m * exp(-l1) - l1 + prop * l1
  der <- m * exp(-l1) - 1 + prop
  l2 <- l1 - fx/der
  i <- 2
  while (abs(l2 - l1) > tol) {
    i <- i + 1
    l1 <- l2
    fx <- m - m * exp(-l1) - l1 + prop * l1
    der <- m * exp(-l1) - 1 + prop
    l2 <- l1 - fx/der
  }
  p <- 1 - m/l2
  loglik <- no * log(p + (1 - p) * exp(-l2)) + n1 * log(1 - 
                                                          p) + sum(dpois(x1, l2, log = TRUE))
  param <- c(l2, p)
  names(param) <- c("lambda", "pi")
  list(iters = i, loglik = loglik, param = param)
}



############################################
zip.mle.NR.EM <- function(x,zhat,lambda,pis){
  dlp <- list()
  dll <- list()
  dlpp <- list()
  dlll <- list()
  dlpl <- list()
  Theta <- list()
  Theta2 <- list()
  
  for (l in 1:L){
    dlp[[l]] <- colSums(zhat[,l]*((Izero(x)*(1-exp(-lambda[[l]]))/(pis[[l]] +(1-pis[[l]])*exp(-lambda[[l]]))) - Inonzero(x)/(1-pis[[l]])))
    dll[[l]] <- colSums(zhat[,l]*((-Izero(x)*(1-pis[[l]])*exp(-lambda[[l]]))/(pis[[l]] +(1-pis[[l]])*exp(-lambda[[l]])) - Inonzero(x) + x/lambda[[l]]))
    dlpp[[l]] <- colSums(zhat[,l]*((-Izero(x)*((1-exp(-lambda[[l]]))^2))/((pis[[l]] +(1-pis[[l]])*exp(-lambda[[l]]))^2) - Inonzero(x)/((1-pis[[l]])^2)))
    dlll[[l]] <- colSums(zhat[,l]*((Izero(x)*pis[[l]]*(1-pis[[l]])*exp(-lambda[[l]]))/((pis[[l]] +(1-pis[[l]])*exp(-lambda[[l]]))^2) - x/(lambda[[l]]^2)))
    dlpl[[l]] <- colSums(zhat[,l]*((Izero(x)*exp(-lambda[[l]]))/((pis[[l]]+(1-pis[[l]])*exp(-lambda[[l]]))^2)))
    # S <- c(dll[[l]], dlp[[l]]) #score function
    Theta[[l]] <- cbind(lambda[[l]], pis[[l]])
    colnames(Theta[[l]]) <- c("lambda","pi")
    HS <- c()
    for (t in 1:d2){
      S <- c(dll[[l]][t], dlp[[l]][t])
      H <-matrix(c(dlll[[l]][t], dlpl[[l]][t], dlpl[[l]][t], dlpp[[l]][t]), 2,2,byrow=T)   #Hessian matrix
      HS <- cbind(HS,solve(H)%*%S)
    }
    Theta2[[l]] <- Theta[[l]] - t(HS)
  }
  
  i <- 0
  tol = 1e-05
  dif <- Inf
  while (sum(dif < tol) == 0) {
    i <- i + 1
    Theta <- Theta2
    for (l in 1:L){
      dlp[[l]] <- colSums(zhat[,l]*((Izero(x)*(1-exp(-lambda[[l]]))/(pis[[l]] +(1-pis[[l]])*exp(-lambda[[l]]))) - Inonzero(x)/(1-pis[[l]])))
      dll[[l]] <- colSums(zhat[,l]*((-Izero(x)*(1-pis[[l]])*exp(-lambda[[l]]))/(pis[[l]] +(1-pis[[l]])*exp(-lambda[[l]])) - Inonzero(x) + x/lambda[[l]]))
      dlpp[[l]] <- colSums(zhat[,l]*((-Izero(x)*((1-exp(-lambda[[l]]))^2))/((pis[[l]] +(1-pis[[l]])*exp(-lambda[[l]]))^2) - Inonzero(x)/((1-pis[[l]])^2)))
      dlll[[l]] <- colSums(zhat[,l]*((Izero(x)*pis[[l]]*(1-pis[[l]])*exp(-lambda[[l]]))/((pis[[l]] +(1-pis[[l]])*exp(-lambda[[l]]))^2) - x/(lambda[[l]]^2)))
      dlpl[[l]] <- colSums(zhat[,l]*((Izero(x)*exp(-lambda[[l]]))/((pis[[l]]+(1-pis[[l]])*exp(-lambda[[l]]))^2)))
      # S <- c(dll[[l]], dlp[[l]]) #score function
      Theta[[l]] <- cbind(lambda[[l]], pis[[l]])
      colnames(Theta[[l]]) <- c("lambda","pi")
      HS <- c()
      for (t in 1:d2){
        S <- c(dll[[l]][t], dlp[[l]][t])
        H <-matrix(c(dlll[[l]][t], dlpl[[l]][t], dlpl[[l]][t], dlpp[[l]][t]), 2,2,byrow=T)   #Hessian matrix
        HS <- cbind(HS,solve(H)%*%S)
      }
      Theta2[[l]] <- Theta[[l]] - t(HS)
    }
    dif = mapply(function(x,y) abs(x-y), x=Theta2, y=Theta)
    print(i)
  }
  
  
  p = Theta[2]; lambda = Theta[1]
  loglik <- n0*log(p + (1-p)*exp(-lambda)) + n1*log(1-p) -n1*lambda + sum(x1)*log(lambda) - sum(log(factorial(x1)))
  param <- t(Theta)
  colnames(param) <- c("lambda", "pi")
  list(iters = i, loglik = loglik, param = param)
  
  
  
  lambda.new <- lapply(Theta2,function(x) t(x[,1]))
  pis.new <- lapply(Theta2,function(x) t(x[,2]))
  return(list(lambda.new = lambda.new, pis.new = pis.new))
}































