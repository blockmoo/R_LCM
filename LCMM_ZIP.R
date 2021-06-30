library(dplyr)
library(tidyverse)
library(stringr)
library(Rfast)


######## 0 step : make dataset
############################
### 00.1.cate variable part
############################
# class 3
y_prob = list(c1p1 = c(0.6,0.3,0.1),
              c2p1 = c(0.2,0.4,0.4),
              c3p1 = c(0.2,0.5,0.3),
              c1p2 = c(0.8,0.2),
              c2p2 = c(0.3,0.7),
              c3p2 = c(0.5,0.5),
              c1p3 = c(0.3,0.65,0.05),
              c2p3 = c(0.6,0.3,0.1),
              c3p3 = c(0.2,0.66,0.14),
              c1p4 = c(0.25,0.25,0.25,0.25),
              c2p4 = c(0.1,0.1,0.6,0.2),
              c3p4 = c(0.2,0.4,0.3,0.1))


## # of categories each observed variable has
mnum = c()
for (i in 1:length(y_prob)) {
  mnum <- c(mnum, length(y_prob[[i]]))
}
mnum

## class prob
r = t(c(0.3,0.6,0.1)) #class 3
# r = t(c(0.4,0.1,0.5)) #class 3
# r = t(c(0.1,0.3,0.5,0.1)) #class 4

L=length(r)
M=length(y_prob) / L
n=3000


#initial_probs
initial_prob <- list()
for(m in 1:M){
  initial_prob[[m]] <- y_prob[((m-1)*L+1) : (m*L)] %>% as.data.frame() %>%  t()
}


## each class dataset
y1 <- vector(mode="list", length = length(y_prob))
for (l in 1:L){
  for (m in 1:M){
    y1[[m + M*(l-1)]] <- sample(1:mnum[m + M*(l-1)], n, replace = TRUE, prob = y_prob[[m + M*(l-1)]])
    print(m + M*(l-1))
  }
}


## combine class data for categorical variables (ex : cbind 4 variables in each class)
y1.dat = vector(mode="list", length = L)
for (l in 1:L){
  ydat = c()
  for (m in 1:M){
    ydat <- cbind(ydat, y1[[L*(m-1)+l]])
  }
  y1.dat[[l]] <- ydat
  colnames(y1.dat[[l]]) <-  paste0('cate',1:M)
}
head(y1.dat[[1]])



#################################
### 00.2.numerical variable part
#################################
## ZIP
poi.x <- function(n,lambda, zero){
  dat <- rep(0,n)
  num = sample(1:n,zero)
  dat[num] <- 0
  dat[-num] <- rpois(n-zero, lambda)
  return(dat)
}

y2.dat= vector(mode="list", length = L)

lamb = list(c1l1 = c(5,1,2,4),
              c2l1 = c(1.5,2,4,7),
              c3l1 = c(8,5,2,1))
zero = 500
for (l in 1:L){
  x <- matrix(poi.x(n,lamb[[l]][1],zero),ncol=1)
  for(i in 1:3){
    x <- cbind(x,poi.x(n,lamb[[l]][i+1],zero))
    colnames(x) <- NULL
  }
  colnames(x) <- paste0('num',1:4)
  y2.dat[[l]] <- x
}

## lambda & pis
apply(y2.dat[[1]],2,zip.mle)-> temp
sapply(temp, function(x) x$param)
# num1      num2      num3      num4
# lambda 4.9451107 1.0387527 2.0443318 3.9941394
# pi     0.2486518 0.2779802 0.2613723 0.2481484

apply(y2.dat[[2]],2,zip.mle) -> temp2
sapply(temp2, function(x) x$param)
# num1      num2      num3      num4
# lambda 1.5387865 2.0612362 4.0498668 7.0076622
# pi     0.2640305 0.2511775 0.2514569 0.2498211

apply(y2.dat[[3]],2,zip.mle) -> temp3
sapply(temp3, function(x) x$param)
# num1      num2      num3      num4
# lambda 8.0073326 4.879461 1.9726067 1.0052794
# pi     0.2497502 0.248790 0.2444008 0.2618967


initial_pi <- initial_lambda <- list()
for(l in 1:L){
  initial_pi[[l]] <- sapply(apply(y2.dat[[l]],2,zip.mle), function(x) x$param)[2,]
  initial_lambda[[l]] <- sapply(apply(y2.dat[[l]],2,zip.mle), function(x) x$param)[1,]
}


##########################################
### 00.3. combine cate & numeric variables
##########################################
y.dat <- list()
for (l in 1:L) {
  y.dat[[l]]  <- cbind(y1.dat[[l]], y2.dat[[l]])
}

## z : class selection matrix 
z <- t(rmultinom(n,size = 1, prob = r))


## new dataset using z & each class's y data 
# If z = (1,0), then get data from class1, elseif (0,1), then from class2
dat <- matrix(0,n,ncol(y.dat[[1]]))
colnames(dat) <- colnames(y.dat[[1]])

for (i in 1:n) {
  dat[i,] <- y.dat[[which(z[i,]==1)]][i,]
}
head(dat)

## trans num to factor (because LCA is used when all var are factor.)
dat <- as.data.frame(dat)
for(m in 1:M) {
  dat[,m] <- as.factor(dat[,m])
}
str(dat)
dim(dat) #1000 * 8


######################
#####  1 step : LCMM_ZIP
formula = cbind(cate1, cate2, cate3, cate4) ~1  # class3
# formula = cbind(cate1, cate2, cate3) ~1         # class4

cate_data = dat[,1:M] 
num_data = dat[,-c(1:M)]

# initial_prob=NULL;initial_lambda=NULL ;initial_pi=NULL
maxiter = 2000; tol = 1e-5; na.rm=TRUEnclass=3; 

LCMM_ZIP_model = function(formula, cate_data, num_data, initial_prob=NULL,initial_lambda=NULL ,nclass, 
                      initial_pi = NULL, maxiter = 2000, tol = 1e-5, na.rm=TRUE){
  mframe <- model.frame(formula, cate_data, na.action = NULL)
  data1 <- cate_data[rowSums(is.na(model.matrix(formula, mframe))) == 0, ]
  if (na.rm) {
    mframe <- model.frame(formula, data1)
    y <- model.response(mframe)
  }else {
    mframe <- model.frame(formula, data1, na.action = NULL)
    y <- model.response(mframe)
    y[is.na(y)] <- 0
  }
  if (any(sapply(lapply(as.data.frame(y), table), length) == 1)) {
    y <- y[, !(sapply(apply(y, 2, table), length) == 1)]
    cat("\n ALERT: at least one manifest variable contained only one\n outcome category, and has been r
        emoved from the analysis. \n\n")
  }
  x <- model.matrix(formula, mframe)
  num_data <- as.matrix(num_data)
  
  
  # #indicate func
  # indicate <- function(x,y) ifelse(x==y,1,0)
  
  diff = 100000
  iter = 1
  Q = 100000
  #llik <- matrix(NA, nrow = maxiter, ncol = 1)
  llik = c()
  llik[1] <- -Inf
  Q.list = c()
  
  ###### initial value #######
  L = nclass
  M = ncol(y)
  d2 = ncol(num_data)
  r.m <- t(matrix(apply(y, 2, max)))
  r = rep(1/nclass,nclass) # class prob
  n = nrow(y)
  
  ## initial prob 
  if (is.null(initial_prob)) {
    probs = list()
    for (m in 1:M) {
      probs[[m]] <- matrix(runif(L * r.m[m]),
                           nrow = L, ncol = r.m[m])
      probs[[m]] <- probs[[m]]/rowSums(probs[[m]])
    }
    probs.new <- probs
  }else {
    probs.new <- probs <- initial_prob
  }
  
  ## initial lambda
  get_lambda <- function(x){
    x1 <- x[x>0]
    sx <- sum(x1)
    m <- sx/n
    s <- (sum(x1^2) - m*sx)/(n-1)
    result <- s/m + m -1
    return(result)
  }
  
  # sam.num <- apply(num_data,2,mean) %>% max() %>%  round(0)
  if (is.null(initial_lambda)) {
    lambda = list()
    for (l in 1:L) {
      ##new2
      # lambda[[l]] <- apply(num_data,2,get_lambda)
      ##new
      lambda[[l]] <- apply(num_data,2,mean)  # sam.num (can change num)
      ##origin
      # lambda[[l]] <- sample(1:sam.num, d2, replace = T)  # sam.num (can change num)
    }
    lambda.new <- lambda
  }else {
    lambda.new <- lambda <- initial_lambda
  }
  
  
  ## initial pi
  if(is.null(initial_pi)){
    pis = list()
    for (l in 1:L) {
      pis[[l]] <- apply(num_data,2,FUN = function(x) sum(x == 0)/n)
    }
    pis.new <- pis
  }else {
    pis.new <- pis <- initial_pi
  }
  
  
  ### change y into dummy
  # make dummy fun
  indi = function(x,y) {
    mat <- rep(0,y)
    mat[x] <- 1
    return(mat)
  }
  
  # y_dummy
  yi <- list()
  flag = 0
  for (k in r.m) {
    flag = flag + 1
    yi[[flag]]<- sapply(y[,flag],function(x) indi(x,k)) %>%  t
  }
  
  #indicator function
  Izero <- function(x){ifelse(x==0,1,0)}
  Inonzero <- function(x){ifelse(x!=0,1,0)}
  
  
  # make z0 matrix
  find.z0 <- function(r, probs, lambda, pis){
    
    ## 01. cate part [rho^I(w=m)] 
    ### calculate sum(I(yi)*probs)
    yi.prob <- list()
    # tprobs <- sapply(probs,t)
    tprobs <- lapply(probs,t)
    
    # mapply() : calculate fun to each element in two list
    yi.prob = mapply(FUN=function(x,y) x%*%y, x = yi, y = tprobs, SIMPLIFY = FALSE) 
    xy <- Reduce("*", yi.prob) # Reduce() : Calculate values corresponding to each element in list
    
    ## 02. numarical part [ZIP의 loglik]
    temp <- matrix(1,n,L) # z0 : numerator of zhat
    
    # #indicator function
    # Izero <- function(x){ifelse(x==0,1,0)}
    # Inonzero <- function(x){ifelse(x!=0,1,0)}
    
    ## origin
    for(l in 1:L) {
      zero <- apply(num_data,1, Izero) %>%  as.data.frame() %>% t()
      zero_term <- log(pis[[l]] +(1-pis[[l]])*exp(-lambda[[l]]))
      temp1 = sweep(zero, MARGIN = 2, zero_term,`*`) # sweep : multiply matrix by vector
      nonzero <- apply(num_data,1, Inonzero) %>%  as.data.frame() %>% t()
      nonzero_term <- apply(num_data,1,FUN=function(x)log(1-pis[[l]]) + x*log(lambda[[l]]) - lambda[[l]] - log(factorial(x))) %>% as.data.frame() %>%  t()
      temp2 <- nonzero * nonzero_term
      temp[,l] = exp(rowSums(temp1 + temp2))
    }
    
    
    ### calculate z0
    # z00 <- matrix(rep(r, nrow(xy)),ncol=L, byrow = T) * xy * temp
    z00 <- sweep(xy * temp, MARGIN = 2, r, "*")
    return(z00)
  }
  
  z0 <- find.z0(r, probs, lambda, pis)
  zhat <- z0 * (1/rowSums(z0))
  
  ### MLE of ZIP using newton raphson in M-step
  # zip.mle.NR.EM <- function(x,zhat,lambda,pis){
  #   dlp <- list()
  #   dll <- list()
  #   dlpp <- list()
  #   dlll <- list()
  #   dlpl <- list()
  #   Theta <- list()
  #   Theta2 <- list()
  #   
  #   for (l in 1:L){
  #     dlp[[l]] <- colSums(zhat[,l]*((Izero(x)*(1-exp(-lambda[[l]]))/(pis[[l]] +(1-pis[[l]])*exp(-lambda[[l]]))) - Inonzero(x)/(1-pis[[l]])))
  #     dll[[l]] <- colSums(zhat[,l]*((-Izero(x)*(1-pis[[l]])*exp(-lambda[[l]]))/(pis[[l]] +(1-pis[[l]])*exp(-lambda[[l]])) - Inonzero(x) + x/lambda[[l]]))
  #     dlpp[[l]] <- colSums(zhat[,l]*((-Izero(x)*((1-exp(-lambda[[l]]))^2))/((pis[[l]] +(1-pis[[l]])*exp(-lambda[[l]]))^2) - Inonzero(x)/((1-pis[[l]])^2)))
  #     dlll[[l]] <- colSums(zhat[,l]*((Izero(x)*pis[[l]]*(1-pis[[l]])*exp(-lambda[[l]]))/((pis[[l]] +(1-pis[[l]])*exp(-lambda[[l]]))^2) - x/(lambda[[l]]^2)))
  #     dlpl[[l]] <- colSums(zhat[,l]*((Izero(x)*exp(-lambda[[l]]))/((pis[[l]]+(1-pis[[l]])*exp(-lambda[[l]]))^2)))
  #     # S <- c(dll[[l]], dlp[[l]]) #score function
  #     Theta[[l]] <- cbind(lambda[[l]], pis[[l]])
  #     colnames(Theta[[l]]) <- c("lambda","pi")
  #     HS <- c()
  #     for (t in 1:d2){
  #       S <- c(dll[[l]][t], dlp[[l]][t])
  #       H <-matrix(c(dlll[[l]][t], dlpl[[l]][t], dlpl[[l]][t], dlpp[[l]][t]), 2,2,byrow=T)   #Hessian matrix
  #       HS <- cbind(HS,solve(H)%*%S)
  #     }
  #     Theta2[[l]] <- Theta[[l]] - t(HS)
  #   }
  #   lambda.new <- lapply(Theta2,function(x) t(x[,1]))
  #   pis.new <- lapply(Theta2,function(x) t(x[,2]))
  #   return(list(lambda.new = lambda.new, pis.new = pis.new))
  # }
  
  
  ### Define constraints of constrOptim function
  # pi_num = unlist(pis) %>% length()
  # lamb_num = unlist(lambda) %>% length()
  # ui <- matrix(0,nrow=(2*pi_num +lamb_num), ncol = (pi_num +lamb_num))
  # for(i in 1:pi_num){
  #   ui[(2*i-1),i] <- 1
  #   ui[(2*i),i] <- -1
  # }
  # 
  # for(j in (2*pi_num+1):(2*pi_num +lamb_num)){
  #   ui[j,j-pi_num] <- 1
  # }
  # ci <- c(rep(c(0,-1),pi_num), rep(0,lamb_num))
  
  
  
  ### EM algorithm
  while (diff > tol & iter < maxiter) {
    
    ### E- step: compute 
    
    # Q 정의
    Q.new = sum(zhat * log(z0))
    
    
    
    ### M = step : update r, p11, p21, p12, p22
    
    #r.new
    r.new = colSums(zhat) / sum(zhat)
    
    
    ## prob.new (new version2)
    # probs.new <- sapply(yi, function(x) (t(zhat) %*% x) * matrix(rep(1/apply(zhat,2, sum),ncol(x)),ncol=ncol(x)))
    probs.new <- lapply(yi, function(x) (t(zhat) %*% x) * matrix(rep(1/apply(zhat,2, sum),ncol(x)),ncol=ncol(x)))
    
    
    # lambda.new <- zip.mle.NR.EM(num_data, zhat, lambda, pis)$lambda.new %>% lapply(as.numeric) 
    # pis.new <- zip.mle.NR.EM(num_data, zhat, lambda, pis)$pis.new %>% lapply(as.numeric) 
    # # pis.new <- lapply(pis.new, function(x) ifelse(x<0,abs(x),x)) # pis.new가 음수일 때 처리
    
    #### estimate (lambda, pis) using "optim"
    # zipzip <- function(pis_lambda_list, probs, r){
    #   # pis = para[1]
    #   # lambda = para[2]
    #   pis = list()
    #   lambda = list()
    #   num = length(pis_lambda_list) / (2*L)
    #   for (l in 1:L) {
    #     pis[[l]] <- pis_lambda_list[((l-1)*num+1):(l*num)]
    #     # print(pis_list[((l-1)*num+1):(l*num)])
    #   }
    #   for (l in (L+1):(2*L)) {
    #     lambda[[l-L]] <- pis_lambda_list[((l-1)*num+1):(l*num)]
    #     # print(pis_list[((l-1)*num+1):(l*num)])
    #   }
    #   # if(sum(unlist(lapply(pis, function(x) x<0))) != 0){
    #   #   return(Inf)
    #   # }
    #   # if(sum(unlist(lapply(pis, function(x) x>1))) != 0){
    #   #   return(Inf)
    #   # }
    #   z0 <- find.z0(r, probs, lambda, pis)
    #   zhat <- z0 * (1/rowSums(z0))
    #   mloglik = -sum(zhat*log(z0))
    #   return(mloglik)
    # } 
    # 
    # ## 1) optim
    # # optim의 default: 최소값 찾기
    # # 최댓값을 가지게 하는 para 찾으려면? control = list(fnscale = -1)
    # # optim(c(0.1, 10), function(y){zipzip(y[1],y[2],1000,x)},control = list(fnscale = -1)) 
    # low = c(rep(0.0001,unlist(pis) %>% length()),rep(0,unlist(lambda) %>% length()))
    # up = c(rep(0.9999,unlist(pis) %>% length()), rep(Inf,unlist(lambda) %>% length()))
    # optim_para = optim(par=unlist(c(pis,lambda)), fn=zipzip, probs = probs, r = r,
    #                    method="L-BFGS-B", lower=low,upper = up)
    # 
    # > system.time(optim(par=unlist(c(pis,lambda)), fn=zipzip, probs = probs, r = r, control=list(fnscale=-1)))
    # 사용자  시스템 elapsed 
    # 77.82    1.68   79.92 
    
    ## 2) constrOptim  -> 너무 오래걸리고 못 찾는듯...ㅎ
    # conoptim <- constrOptim(theta=unlist(c(pis,lambda)), f=zipzip, ui=ui, ci=ci, grad=NULL, probs = probs, r = r)
    
    
    #### optim 2차원(lambda, pi)에서만 돌아가도록 코드 수정 in 06/30
    zip_optim <- function( lambda_pis, zhat, num_data){
      lambda = lambda_pis[1]
      pis = lambda_pis[2]
      zero <- Izero(num_data)
      zero_term <- log(pis +(1-pis)*exp(-lambda))
      temp1 <- zero * zero_term
      nonzero <- Inonzero(num_data)
      nonzero_term <- log(1-pis) + num_data*log(lambda) - lambda
      temp2 <- nonzero * nonzero_term
      temp = sum(zhat * (temp1 + temp2))
      return(-temp)
    }
    
    
    lambda.new <- lambda
    pis.new <- pis
    for (i in 1:d2){
      for (l in 1:L){
        optim_para = optim(par=c(lambda[[l]][i],pis[[l]][i]), fn=zip_optim, zhat = zhat[,l], num_data = num_data[,i],
                           method="L-BFGS-B", lower=c(0,0.0001),upper = c(Inf,0.9999))
        lambda.new[[l]][i] <- optim_para$par[1]
        pis.new[[l]][i] <- optim_para$par[2]
      }
    }
    
    # 
    # pis.new = list()
    # lambda.new = list()
    # pis_lambda_list = unlist(c(pis,lambda))
    # num = length(pis_lambda_list) / (2*L)
    # for (l in 1:L) {
    #   pis.new[[l]] <- optim_para$par[((l-1)*num+1):(l*num)]
    #   # print(pis_list[((l-1)*num+1):(l*num)])
    # }
    # for (l in (L+1):(2*L)) {
    #   lambda.new[[l-L]] <- optim_para$par[((l-1)*num+1):(l*num)]
    #   # print(pis_list[((l-1)*num+1):(l*num)])
    # }
    
    
    z0.new <- find.z0(r.new, probs.new, lambda.new, pis.new)
    zhat.new <- z0.new * (1/rowSums(z0.new))
    
    
    # calculate observed loglikelihood to check convergence
    llik[iter+1] = sum(log(rowSums(z0.new)))
    diff = llik[iter+1] - llik[iter]
    if(diff < 0){
      break
    }
    
    
    # update para
    r = r.new
    probs = probs.new
    lambda = lambda.new
    pis = pis.new
    zhat = zhat.new
    z0 = z0.new
    Q = Q.new
    
    
    # cat("Iter", iter, ": r = ",r.new, ", llik = ",llik[iter+1],", diff = ", diff, '\n')
    
    ### new version (06/30)
    cat("Iter", iter, ": r = ",r , ", llik = ",llik[iter+1],", diff = ", diff, '\n')
    
    Q.list = c(Q.list, Q)
    iter=iter+1
    
  }
  ret = list()
  ret$nclass <- L
  ret$npara <- (L * sum(r.m - 1)) + (L - 1)
  ret$aic <- (-2 * tail(llik,1)) + (2 * ret$npar)
  ret$bic <- (-2 * tail(llik,1)) + (log(n) * ret$npar)
  ret$caic <- (-2 * tail(llik,1)) + ret$npar * (log(n) + 1)
  
  # prob
  for(m in 1:M){
    colnames(probs[[m]])<- 1:r.m[m]
    rownames(probs[[m]]) <- paste('class',1:L,":", sep = ' ')
  }
  names(probs.new) <- colnames(y)
  
  # lambda
  for(l in 1:L){
    names(lambda[[l]])<- 1:d2
    names(lambda)[l] <- paste('class',l)
  }
  
  # pi
  for(l in 1:L){
    names(pis[[l]])<- 1:d2
    names(pis)[l] <- paste('class',l)
  }
  
  # return(list(r=r.new, ollik = llik, probs = probs.new,lambda = lambda.new, pi = pis.new, Q = Q.list, z= zhat,
  #             ret = ret, diff=diff))
  
  ### new version (06/30)
  return(list(r=r, ollik = llik[1:iter], probs = probs, lambda = lambda, pi = pis, Q = Q.list,
              z= zhat, ret = ret))
}


lcmm_zip <- LCMM_ZIP_model(formula, cate_data, num_data, initial_prob=initial_prob,initial_lambda=initial_lambda,
                   initial_pi =initial_pi,nclass = 3, maxiter = 2000, tol = 1e-5, na.rm=TRUE)


lcmm_zip$ollik %>%  tail(1)
lcmm_zip$r
lcmm_zip$probs
lcmm_zip$lambda
lcmm_zip$pi


system.time(
  lcmm_zip <- LCMM_ZIP_model(formula, cate_data, num_data, initial_prob=initial_prob,initial_lambda=initial_lambda,
                             initial_pi =initial_pi,nclass = 3, maxiter = 2000, tol = 1e-5, na.rm=TRUE)
)
# 사용자  시스템 elapsed 
# 20.22    1.07   21.66 

save(num_data, cate_data, dat, LCMM_ZIP_model, lcmm_zip, file="LCMM_ZIP_0617.r")
save(y_prob, initial_lambda, initial_pi, optim_1000, optim_2000, file="LCMM_ZIP_0618.r")

load("LCMM_ZIP_0615.r")

###### optim 2차원으로 돌린 코드 (06/30)------------------------------




#---------------------------------------------------------------------
##### 병렬 처리 : simul 10번
library(foreach)
library(parallel)
library(doParallel)
cl <- makeCluster(2) #not to overload your computer
registerDoParallel(cl) # Ready to parallel
simul_times = 100
simulation3 <- list()

simul_LCMM_zip0630 <- foreach(simul = 1:simul_times) %dopar% {
  library(dplyr)
  library(tidyverse)
  library(stringr)
  library(Rfast)
  
  simulation3[[simul]] <- LCMM_ZIP_model(formula, cate_data, num_data, initial_prob=initial_prob,initial_lambda=initial_lambda,
                                         initial_pi =initial_pi,nclass = 3, maxiter = 2000, tol = 1e-5, na.rm=TRUE)
  
  simulation3
  
  
}
stopCluster(cl)

save(simul_LCMM_zip0630 , num_data, cate_data, LCMM_ZIP_model, y_prob, initial_lambda, 
     initial_pi, initial_prob, file="LCMM_ZIP_simul10.R")
load( file="LCMM_ZIP_simul10.R")

save(simul_low_add_decomp0627 , num_data, cate_data, LCMM_ZIP_model, y_prob, initial_lambda, 
     initial_pi, initial_prob, file="LCMM_ZIP_simul0627.R")
load( file="LCMM_ZIP_simul0627.R")

### simul mean values
# simul_r <- sapply(simul_low_add_decomp, FUN= function(y) sapply(y,FUN= function(x) sort(x$r)) %>%  unlist %>% matrix(ncol=3, byrow = T) %>% colMeans()) %>% rowMeans()
# r순서대로 prob, lambda, pi값 가져올 수 있도록,
# order_r <- sapply(simul_low_add_decomp[[10]],FUN= function(x) if(is.null(x) == F)order(x$r)) 
# simul_low_add_decomp[[10]][[10]]$probs %>% sapply(function(x) x[aaaa,])
# order_prob <- mapply(y=simul_low_add_decomp, z=order_r, FUN= function(y,z) sapply(y,FUN= function(x) {x$probs %>% sapply(function(w) w[z,])})) #%>%  unlist %>% matrix(ncol=3, byrow = T) %>% colMeans()) %>% rowMeans()
# sapply(order_prob, function(x) sapply(Reduce("+",x)))

result <- list()

for (i in 1:10) {
  # result[[i]] <- simul_low_add_decomp[[i]][[i]]
  result[[i]] <- simul_low_add_decomp0627[[i]][[i]]
} 
#### simul_r
simul_r <- sapply(result,FUN= function(x) sort(x$r)) %>% rowMeans()
order_r <- sapply(result,FUN= function(x) order(x$r))
col_list <- list()
for (i in 1:10) {
  col_list[[i]] <- order_r[,i]
}

#### simul_prob
result[[1]]$probs[]
sapply(result[[2]]$probs, function(x) x[col_list[[2]],]) -> a
sapply(result[[1]]$probs, function(x) x[col_list[[1]],]) -> b
temp <- list()
for (j in 1:10){
  temp[[j]] <- sapply(result[[j]]$probs, function(x) x[col_list[[j]],])
}
simul_prob <- temp[[1]]
for (i in 2:10){
  simul_prob[[1]] <- simul_prob[[1]] + temp[[i]][[1]]
  simul_prob[[2]] <- simul_prob[[2]] + temp[[i]][[2]]
  simul_prob[[3]] <- simul_prob[[3]] + temp[[i]][[3]]
  simul_prob[[4]] <- simul_prob[[4]] + temp[[i]][[4]]
}

for (i in 1:4){
  simul_prob[[i]] <- simul_prob[[i]]/10
}
simul_prob

#### simul_lambda
temp <- list()
for (j in 1:10){
  temp[[j]] <- result[[j]]$lambda[col_list[[j]]]
}

simul_lambda <- temp[[1]]
for (i in 2:10){
  simul_lambda[[1]] <- simul_lambda[[1]] + temp[[i]][[1]]
  simul_lambda[[2]] <- simul_lambda[[2]] + temp[[i]][[2]]
  simul_lambda[[3]] <- simul_lambda[[3]] + temp[[i]][[3]]
}

for (i in 1:3){
  simul_lambda[[i]] <- simul_lambda[[i]]/10
}
simul_lambda

#### simul_pis
temp <- list()
for (j in 1:10){
  temp[[j]] <- result[[j]]$pi[col_list[[j]]]
}
simul_pi <- temp[[1]]
for (i in 2:10){
  simul_pi[[1]] <- simul_pi[[1]] + temp[[i]][[1]]
  simul_pi[[2]] <- simul_pi[[2]] + temp[[i]][[2]]
  simul_pi[[3]] <- simul_pi[[3]] + temp[[i]][[3]]
}

for (i in 1:3){
  simul_pi[[i]] <- simul_pi[[i]]/10
}
simul_pi

#### simul_loglik
simul_loglik <- matrix(0,nrow=10)
for (i in 1:10){
  simul_loglik[i,1] <- result[[i]]$ollik %>% tail(1)
}
simul_loglik


##initial values
initial_pi
initial_lambda
initial_prob




















#################################################
### zip.mle : compute MLE for ZIP using newton raphson method
# install.packages("Rfast")
library(Rfast)

x <- rpois(100, 10)
num = sample(1:100,30)
x[num] <- 0
zip.mle(x)


##### NR in M-step (ZIP)
zip.mle.NR.EM <- function (x, zhat, lambda, pi, tol = 1e-09){
  #indicator function
  Izero <- function(x){ifelse(x==0,1,0)}
  Inonzero <- function(x){ifelse(x!=0,1,0)}
  
  #Newton Raphson
  dlp <- zhat*((Izero(x)*(1-exp(-lambda))/(pi +(1-pi)*exp(-lambda))) - Inonzero(x)/(1-pi))
  dll <- zhat*((-Izero(x)*(1-pi)*exp(-lambda))/(pi +(1-pi)*exp(-lambda)) - Inonzero(x) + x/lambda)
  dlpp <- zhat*((-Izero(x)*((1-exp(-lambda))^2))/((pi +(1-pi)*exp(-lambda))^2) - Inonzero(x)/((1-pi)^2))
  dlll <- zhat*((Izero(x)*pi*(1-pi)*exp(-lambda))/((pi +(1-pi)*exp(-lambda))^2) - x/(lambda^2))
  dlpl <- zhat*((Izero(x)*exp(-lambda))/((pi+(1-pi)*exp(-lambda))^2))
  S <- c(dll, dlp) #score function
  H <-matrix(c(dlll, dlpl, dlpl, dlpp), 2,2,byrow=T)   #Hessian matrix
  Theta <- c(lambda, pi)
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
num = sample(1:1000,600)
x[num] <- 0
zip.mle(x)
zip.mle.NR(x)

x <- rpois(100, 10)
num = sample(1:100,20)
x[num] <- 0



#######################
#### optim function ###
###
zipzip <- function(para,n,data){
  p = para[1]
  lambda = para[2]
  y <- data
  n0 <- sum(y == 0)
  n1 <- n - n0
  if(p < 0 | p > 1){
    return(Inf)
  }
  loglik = n0*log(p+(1-p)*exp(-lambda)) + n1*log(1-p) -n1*lambda + n*mean(y)*log(lambda) - sum(log(factorial(y)))
  return(loglik)
} 

# optim의 default: 최소값 찾기
# 최댓값을 가지게 하는 para 찾으려면? control = list(fnscale = -1)
# optim(c(0.1, 10), function(y){zipzip(y[1],y[2],1000,x)},control = list(fnscale = -1)) 
optim(par=c(0.1,2), fn=zipzip, n=1000, data=num_data[,1], control=list(fnscale=-1))
zip.mle(num_data[,1])



########### LCMM 수정(06/14 ~)
# estimate (lambda, pis) using "optim"
zipzip <- function(pis_list, lambda, probs, zhat){
  # pis = para[1]
  # lambda = para[2]
  pis = list()
  num = length(pis_list) / L
  for (l in 1:L) {
    pis[[l]] <- pis_list[((l-1)*num+1):(l*num)]
    # print(pis_list[((l-1)*num+1):(l*num)])
  }
  z0 <- find.z0(r, probs, lambda, pis)
  if(sum(unlist(lapply(pis, function(x) x<0))) != 0){
    return(Inf)
  }
  if(sum(unlist(lapply(pis, function(x) x>1))) != 0){
    return(Inf)
  }
  loglik = sum(zhat*log(z0))
  return(loglik)
} 

# optim의 default: 최소값 찾기
# 최댓값을 가지게 하는 para 찾으려면? control = list(fnscale = -1)
# optim(c(0.1, 10), function(y){zipzip(y[1],y[2],1000,x)},control = list(fnscale = -1)) 
optim(par=unlist(pis), fn=zipzip, lambda = lambda, probs = probs, zhat = zhat, control=list(fnscale=-1))


###############
# estimate (lambda, pis) using "optim"
zipzip <- function(pis_lambda_list, probs, zhat){
  # pis = para[1]
  # lambda = para[2]
  pis = list()
  lambda = list()
  num = length(pis_lambda_list) / (2*L)
  for (l in 1:L) {
    pis[[l]] <- pis_lambda_list[((l-1)*num+1):(l*num)]
    # print(pis_list[((l-1)*num+1):(l*num)])
  }
  for (l in (L+1):(2*L)) {
    lambda[[l-L]] <- pis_lambda_list[((l-1)*num+1):(l*num)]
    # print(pis_list[((l-1)*num+1):(l*num)])
  }
  z0 <- find.z0(r, probs, lambda, pis)
  if(sum(unlist(lapply(pis, function(x) x<0))) != 0){
    return(Inf)
  }
  if(sum(unlist(lapply(pis, function(x) x>1))) != 0){
    return(Inf)
  }
  loglik = sum(zhat*log(z0))
  return(loglik)
} 

# optim의 default: 최소값 찾기
# 최댓값을 가지게 하는 para 찾으려면? control = list(fnscale = -1)
# optim(c(0.1, 10), function(y){zipzip(y[1],y[2],1000,x)},control = list(fnscale = -1)) 
optim_para = optim(par=unlist(c(pis,lambda)), fn=zipzip, probs = probs, zhat = zhat, control=list(fnscale=-1))
pis.new = list()
lambda.new = list()
num = length(pis_lambda_list) / (2*L)
for (l in 1:L) {
  pis.new[[l]] <- optim_para$par[((l-1)*num+1):(l*num)]
  # print(pis_list[((l-1)*num+1):(l*num)])
}
for (l in (L+1):(2*L)) {
  lambda.new[[l-L]] <- optim_para$par[((l-1)*num+1):(l*num)]
  # print(pis_list[((l-1)*num+1):(l*num)])
}