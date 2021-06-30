## LCA function
## L : # of class  = 2 (l = 1,2)
## M : # of variable   -> multi : 4 with 3 values
## r.m : # of categories each observed variable has 
## n : obs 
# feedback : Don't use round func, when we estimate zhat & r


library(dplyr)
library(poLCA)


### 0.create Dataset 
#### 0.1 cateforical variables  (560th line : numerical variables)
#### True value
# class = 2, # of y : 5 with 3 categories
y_prob = list(c1p1 = c(0.6,0.3,0.1),
              c2p1 = c(0.2,0.4,0.4),
              c1p2 = c(0.2,0.5,0.3),
              c2p2 = c(0.8,0.2),
              c1p3 = c(0.3,0.7),
              c2p3 = c(0.5,0.5),
              c1p4 = c(0.3,0.65,0.05),
              c2p4 = c(0.6,0.3,0.1),
              c1p5 = c(0.2,0.66,0.14),
              c2p5 = c(0.25,0.25,0.25,0.25),
              c1p6 = c(0.1,0.1,0.6,0.2),
              c2p6 = c(0.2,0.4,0.3,0.1))


# # of categories each observed variable has
mnum = c()
for (i in 1:length(y_prob)) {
  mnum <- c(mnum, length(y_prob[[i]]))
}
mnum

# class prob
r = t(c(0.2,0.7,0.1))
#r = t(c(0.4,0.6))
L=length(r)
M=length(y_prob) / L
n=1000

# each class dataset
yy <- vector(mode="list", length = length(y_prob))
for (l in 1:L){
  for (m in 1:M){
    yy[[m + M*(l-1)]] <- sample(1:mnum[m + M*(l-1)], n, replace = TRUE, prob = y_prob[[m + M*(l-1)]])
    print(m + M*(l-1))
  }
}


# combine class data (ex : cbind 4 variables in each class)
cl = vector(mode="list", length = L)
for (l in 1:L){
  ydat = c()
  for (m in 1:M){
    ydat <- cbind(ydat, yy[[L*(m-1)+l]])
  }
  cl[[l]] <- ydat
}
head(cl[[1]])


# z : class selection matrix 
z <- t(rmultinom(n,size = 1, prob = r))


# new dataset using z & each class's y data 
# If z = (1,0), then get data from class1, elseif (0,1), then from class2
dat <- matrix(0,n,ncol(cl[[1]]))
for (i in 1:n) {
  dat[i,] <- cl[[which(z[i,]==1)]][i,] 
}
head(dat)


## trans num to factor (because LCA is used when all var are factor.)
data <- as.data.frame(dat)

for (i in 1:ncol(data)){
  data[,i] <- as.factor(data[,i])
}
head(data)



### poLCA(04/09)
names(data) # var name in cbind
f <- cbind(V1,V2,V3,V4) ~ 1 
#f <- cbind(V1,V2,V3,V4,V5,V6) ~ 1 
# lc <- poLCA(f, data, nclass=2, maxiter=2000, tol=1e-5)


###############################################

formula = cbind(V1,V2,V3,V4) ~ 1 
na.rm=TRUE


############################################################
########## 01-2. my function new version
new_LCM = function(formula, data, initial_prob=NULL ,nclass, 
                     maxiter = 2000, tol = 1e-5, na.rm=TRUE){
  mframe <- model.frame(formula, data, na.action = NULL)
  data1 <- data[rowSums(is.na(model.matrix(formula, mframe))) == 0, ]
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
  
  #indicate func
  indicate <- function(x,y) ifelse(x==y,1,0)
  
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
  }
  else {
    probs.new <- probs <- initial_prob
  }
  
  ## initial zhat = p(zi=l |Y,r,rhomk)
  z0 <- matrix(1,n,L) # z0 : numerator of zhat
  
  for (l in 1:L){
    sum.probs = rep(0,dim(y)[1])
    for (i in 1:n){
      for (m in 1:M){
        sum.probs[i] = sum.probs[i] + log(probs[[m]][l,y[i,m]])
      }
      z0[i,l] = exp(log(r[l]) + sum.probs[i])
    }
  }
  zhat <- z0 * (1/rowSums(z0))

  
  ### EM algorithm
  while (diff > tol & iter < maxiter) {
    
    ### E- step: compute 
    
    # Q 정의
    Q.new = sum(zhat * log(z0))
    
    
    
    ### M = step : update r, p11, p21, p12, p22
    
    #r.new
    r.new = colSums(zhat) / sum(zhat)
    
   
    #prob.new (new version)
    for (l in 1:L) {
      for (m in 1:M) {
        for(k in 1:r.m[m]) {
          probs.new[[m]][l,k] = sum(zhat[,l]*indicate(y[,m],k)) / sum(zhat[,l])
        }
      }
    }
    
    
    # new zhat
    z0.new <- matrix(1,n,L) # z0.new : numerator of zhat.new
    
    for (l in 1:L){
      sum.probs = rep(0,dim(y)[1])
      for (i in 1:n){
        for (m in 1:M){
          sum.probs[i] = sum.probs[i] + log(probs.new[[m]][l,y[i,m]])
        }
        z0.new[i,l] = exp(log(r.new[l]) + sum.probs[i])
      }
    }
    zhat.new <- z0.new * (1/rowSums(z0.new))
    
    
    # calculate observed loglikelihood to check convergence
    llik[iter+1] = sum(log(rowSums(z0.new)))
    diff = llik[iter+1] - llik[iter]
    
    
    # update para
    r = r.new
    probs = probs.new
    zhat = zhat.new
    z0 = z0.new
    Q = Q.new
    
    
    cat("Iter", iter, ": r = ",r.new, ", llik = ",llik[iter+1],", diff = ", diff, '\n')
    
    Q.list = c(Q.list, Q)
    iter=iter+1
    
  }
  ret = list()
  ret$nclass <- L
  ret$npara <- (L * sum(r.m - 1)) + (L - 1)
  ret$aic <- (-2 * tail(llik,1)) + (2 * ret$npar)
  ret$bic <- (-2 * tail(llik,1)) + (log(n) * ret$npar)
  ret$caic <- (-2 * tail(llik,1)) + ret$npar * (log(n) + 1)
  
  
  for(m in 1:M){
    colnames(probs[[m]])<- 1:r.m[m]
    rownames(probs[[m]]) <- paste('class',1:L,":", sep = ' ')
  }
  names(probs.new) <- colnames(y)

  return(list(r=r.new, ollik = llik, probs = probs, Q = Q.list, z= zhat,
              ret = ret))
}


############################################################
########## 01-3. my function new version (until 04/29)
new_LCM2 = function(formula, data, initial_prob=NULL ,nclass, 
                   maxiter = 2000, tol = 1e-5, na.rm=TRUE){
  mframe <- model.frame(formula, data, na.action = NULL)
  data1 <- data[rowSums(is.na(model.matrix(formula, mframe))) == 0, ]
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
  
  #indicate func
  indicate <- function(x,y) ifelse(x==y,1,0)
  
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
  }
  else {
    probs.new <- probs <- initial_prob
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
  
  
  # make z0 matrix
  find.z0 <- function(r, probs){
    
    ### calculate sum(I(yi)*probs)
    yi.prob <- list()
    # tprobs <- sapply(probs,t)
    tprobs <- lapply(probs,t)
    
    
    # mapply() : calculate fun to each element in two list
    yi.prob = mapply(FUN=function(x,y) x%*%y, x = yi, y = tprobs, SIMPLIFY = FALSE) 
    # => same result in for version
    # for (m in 1:M){
    #   yi.prob[[m]] <- yi[[m]] %*% tprobs[[m]]
    # }
    # xy <- apply(simplify2array(yi.prob), c(1,2), prod)
    xy <- Reduce("*", yi.prob) # Reduce() : Calculate values corresponding to each element in list
    
    ### calculate z0
    z00 <- matrix(rep(r, nrow(xy)),ncol=L, byrow = T) * xy
    return(z00)
  }
  
  z0 <- find.z0(r, probs)
  zhat <- z0 * (1/rowSums(z0))
  
  
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
    
    
    #(new version1)
    # for (m in 1:M) {
    #   probs.new[[m]] <- (t(zhat) %*% yi[[m]]) *matrix(rep(1/apply(zhat,2, sum),ncol(yi[[m]])),ncol=ncol(yi[[m]]))
    # }
    
    #(old version)
    # for (l in 1:L) {
    #   for (m in 1:M) {
    #     for(k in 1:r.m[m]) {
    #       probs.new[[m]][l,k] = sum(zhat[,l]*indicate(y[,m],k)) / sum(zhat[,l])
    #     }
    #   }
    # }
    
    
    z0.new <- find.z0(r.new, probs.new)
    zhat.new <- z0.new * (1/rowSums(z0.new))
    
    
    # calculate observed loglikelihood to check convergence
    llik[iter+1] = sum(log(rowSums(z0.new)))
    diff = llik[iter+1] - llik[iter]
    
    
    # update para
    r = r.new
    probs = probs.new
    zhat = zhat.new
    z0 = z0.new
    Q = Q.new
    
    
    cat("Iter", iter, ": r = ",r.new, ", llik = ",llik[iter+1],", diff = ", diff, '\n')
    
    Q.list = c(Q.list, Q)
    iter=iter+1
    
  }
  ret = list()
  ret$nclass <- L
  ret$npara <- (L * sum(r.m - 1)) + (L - 1)
  ret$aic <- (-2 * tail(llik,1)) + (2 * ret$npar)
  ret$bic <- (-2 * tail(llik,1)) + (log(n) * ret$npar)
  ret$caic <- (-2 * tail(llik,1)) + ret$npar * (log(n) + 1)
  
  
  for(m in 1:M){
    colnames(probs[[m]])<- 1:r.m[m]
    rownames(probs[[m]]) <- paste('class',1:L,":", sep = ' ')
  }
  names(probs.new) <- colnames(y)
  
  return(list(r=r.new, ollik = llik, probs = probs, Q = Q.list, z= zhat,
              ret = ret))
}



######################################  
########### 02. compare poLCA & my LCA
### 1) my LCA result
f <- cbind(V1,V2,V3,V4) ~ 1 
new <- LCM_model(f, data, initial_prob=NULL ,nclass = 3, 
                 maxiter = 2000, tol = 1e-5, na.rm=TRUE)
system.time(new <- LCM_model(f, data, initial_prob=NULL ,nclass = 3, 
                      maxiter = 2000, tol = 1e-5, na.rm=TRUE))
tail(new$ollik,1)
new$r
new$probs

set.seed(1234)
# version2 (04/22)
system.time(newnew <- new_LCM(f, data, initial_prob=NULL ,nclass = 3, 
                              maxiter = 2000, tol = 1e-5, na.rm=TRUE))

tail(newnew$ollik,1)
newnew$r
newnew$probs


# version3 (04/29)
system.time(new0426 <- new_LCM2(f, data, initial_prob=NULL ,nclass = 3, 
                              maxiter = 2000, tol = 1e-5, na.rm=TRUE))

tail(new0426$ollik,1)
new0426$r
new0426$probs


### 2) poLCA result
lc <- poLCA(f, data, nclass=3, maxiter=2000, tol=1e-5)
system.time(lc <- poLCA(f, data, nclass=3, maxiter=2000, tol=1e-5))

lc$llik
lc$P
lc$probs

### 3) version1 insert result of poLCA into my LCA
insert_po_my <- LCM_model(f, data, initial_prob=lc$probs ,nclass = 3, 
                 maxiter = 2000, tol = 1e-5, na.rm=TRUE)
system.time(insert_po_my <- LCM_model(f, data, initial_prob=lc$probs ,nclass = 3, 
                                 maxiter = 2000, tol = 1e-5, na.rm=TRUE))
insert_po_my$ret

tail(insert_po_my$ollik,1)
insert_po_my$r
insert_po_my$probs

# version2
system.time(insert_po_my2 <- new_LCM(f, data, initial_prob=lc$probs ,nclass = 3, 
                      maxiter = 2000, tol = 1e-5, na.rm=TRUE))


tail(insert_po_my2$ollik,1)
insert_po_my2$r
insert_po_my2$probs

################# compare time
system.time(lc <- poLCA(f, data, nclass=3, maxiter=2000, tol=1e-5))
system.time(newnew <- new_LCM(f, data, initial_prob=NULL ,nclass = 3, 
                              maxiter = 2000, tol = 1e-5, na.rm=TRUE))
system.time(new0426 <- new_LCM2(f, data, initial_prob=NULL ,nclass = 3, 
                                maxiter = 2000, tol = 1e-5, na.rm=TRUE))

system.time(new <- LCM_model(f, data, initial_prob=NULL ,nclass = 3, 
                             maxiter = 2000, tol = 1e-5, na.rm=TRUE))


#########################################
############ 04/26 update

### change y into dummy
# make dummy fun
indi = function(x,y) {
  mat <- rep(0,y)
  mat[x] <- 1
  return(mat)
}


# make z0 matrix
find.z0 <- function(y, r, probs){
  
  # y_dummy
  yi <- list()
  tprobs <- sapply(probs,t)
  
  iter = 0
  for (k in r.m) {
    iter = iter + 1
    yi[[iter]]<- sapply(y[,iter],function(x) indi(x,k)) %>%  t
  }
  
  ### caculate sum(I(yi) *probs)
  yi.prob <- list()
  
  
  for (m in 1:M){
    yi.prob[[m]] <- yi[[m]] %*% tprobs[[m]]
    # yi[[m]] <- yi[[m]] %*% tprobs[[m]]
  }
  xy <- apply(simplify2array(yi.prob), c(1,2), prod)
  # xy <- apply(simplify2array(yi), c(1,2), prod)

  ### calculate z0
  z00 <- matrix(rep(r, nrow(xy)),ncol=L, byrow = T) * xy
  return(z0)
}



a <- function(x){
  z0 <- matrix(1,n,L) # z0 : numerator of zhat
  
  for (l in 1:L){
    sum.probs = rep(0,dim(y)[1])
    for (i in 1:n){
      for (m in 1:M){
        sum.probs[i] = sum.probs[i] + log(probs[[m]][l,y[i,m]])
      }
      z0[i,l] = exp(log(r[l]) + sum.probs[i])
    }
  }
  return(z0)
}

system.time(find.z0(y,probs))
system.time(a(x))


############################################
#### 03. numerical variables (update 05/03)
# x1 ~ poi(5), x2 ~ poi(10), x3 ~ poi(3)
lambda = list(c1l1 = c(5,10,3),
              c2l1 = c(1.5,4,2),
              c3l1 = c(8,2,4))
lambda = list(c1l1 = c(5),
              c2l1 = c(10))
lambda = list(c1l1 = c(5,10,30,20),
              c2l1 = c(1.5,4,20,20),
              c3l1 = c(8,2,40,20))

## class prob
r = t(c(0.3,0.6,0.1))
r = t(c(0.25, 0.75))

L=length(r)
d2=length(lambda[[1]])
n=1000

## each class dataset
y2 <- vector(mode="list", length = d2)
temp <- matrix(0, nrow=n, ncol=d2)
for (l in 1:L){
  for (t in 1:d2){
    temp[,t] <- rpois(n,lambda[[l]][t])
  }
  y2[[l]] <- temp
}

# -> check up
# > y2[[1]] %>% colSums() / n
# [1] 4.967 9.827 2.989
# > y2[[2]] %>% colSums() / n
# [1] 1.564 4.089 1.990
# > y2[[3]] %>% colSums() / n
# [1] 8.092 1.934 3.984


## z : class selection matrix 
z <- t(rmultinom(n,size = 1, prob = r)) 
#-> check up
# > colSums(z) / n
# [1] 0.294 0.609 0.097

## new dataset using z & each class's y data 
# If z = (1,0), then get data from class1, elseif (0,1), then from class2
dat <- matrix(0,n,ncol=d2)
for (i in 1:n) {
  dat[i,] <- y2[[which(z[i,]==1)]][i,] 
}
head(dat)


#############################
#### 04. my function for numerical latent model (until 05/14)
LNM_model = function(data, initial_lambda=NULL ,nclass,
                     maxiter = 2000, tol = 1e-5, na.rm=TRUE){
  
  diff = 100000
  iter = 1
  Q = 100000
  llik = c()
  llik[1] <- -Inf
  Q.list = c()
  
  ###### initial value #######
  L = nclass
  d2 = ncol(data)
  r = rep(1/nclass,nclass) # class prob
  n = nrow(data)
  
  ## initial lambda
  if (is.null(initial_lambda)) {
    lambda = list()
    for (l in 1:L) {
      lambda[[l]] <- sample(1:200, d2, replace = T)
    }
    lambda.new <- lambda
  }
  else {
    lambda.new <- lambda <- initial_lambda
  }
  
  ## initial zhat = p(zi=l |Y,r,rhomk)
  z0 <- matrix(1,n,L) # z0 : numerator of zhat
  
  for(i in 1:n){
    for(l in 1:L){
      for(t in 1:d2) {
        z0[i,l] = z0[i,l]*(exp(-lambda[[l]][t])*(lambda[[l]][t]^data[i,t])/factorial(data[i,t]))
      }
      # z0[i,l] = r[l] * z0[i,l]
    }
  }
  zhat <- z0 * (1/rowSums(z0))
  
  
  ### EM algorithm
  while (diff > tol & iter < maxiter) {
    
    ### E- step: compute
    
    # Q 정의
    Q.new = sum(zhat * log(z0))
    
    ### M = step : update r, p11, p21, p12, p22
    
    #r.new
    r.new = colSums(zhat) / sum(zhat)
    # r.new = colSums(zhat) / sum(zhat)
    
    
    for (l in 1:L) {
      for (t in 1:d2) {
        lambda.new[[l]][t] = sum(zhat[,l]*data[,t]) / sum(zhat[,l])
        
      }
    }
    
    
    # new zhat
    z0.new <- matrix(1,n,L) # z0 : numerator of zhat
    
    for(i in 1:n){
      for(l in 1:L){
        for(t in 1:d2) {
          z0.new[i,l] = z0.new[i,l]*(exp(-lambda.new[[l]][t])*((lambda.new[[l]][t]^data[i,t])/factorial(data[i,t])))
        }
        # z0.new[i,l] = r.new[l] * z0.new[i,l]
      }
    }
    zhat.new <- z0.new * (1/rowSums(z0.new))
    
    
    # calculate observed loglikelihood to check convergence
    llik[iter+1] = sum(log(rowSums(z0.new)))
    diff = llik[iter+1] - llik[iter]
    
    
    # update para
    r = r.new
    lambda = lambda.new
    zhat = zhat.new
    z0 = z0.new
    Q = Q.new
    
    
    cat("Iter", iter, ": r = ",r.new, ", llik = ",llik[iter+1],", diff = ", diff, '\n')
    
    Q.list = c(Q.list, Q)
    iter=iter+1
    
  }
  ret = list()
  ret$nclass <- L
  ret$npara <- (L * d2) + (L - 1)
  ret$aic <- (-2 * tail(llik,1)) + (2 * ret$npar)
  ret$bic <- (-2 * tail(llik,1)) + (log(n) * ret$npar)
  ret$caic <- (-2 * tail(llik,1)) + ret$npar * (log(n) + 1)
  
  # lambda
  for(l in 1:L){
    names(lambda[[l]])<- 1:d2
    names(lambda)[l] <- paste('class',l)
  }
  
  return(list(r=r.new, ollik = llik, lambda = lambda, Q = Q.list, z= zhat,
              ret = ret))
}



lnm <- LNM_model(data = dat, initial_lambda=NULL ,nclass = 3,
                 maxiter = 2000, tol = 1e-5, na.rm=TRUE)

lnm$ollik %>%  tail(1)
lnm$r
lnm$lambda



############################################
#### 05. mixed variables (update 05/03)









############################################
######## historical code

########01. my LCM function
# LCM_model = function(formula, data, initial_prob=NULL ,nclass, 
#                      maxiter = 2000, tol = 1e-5, na.rm=TRUE){
#   mframe <- model.frame(formula, data, na.action = NULL)
#   data1 <- data[rowSums(is.na(model.matrix(formula, mframe))) == 0, ]
#   if (na.rm) {
#     mframe <- model.frame(formula, data1)
#     y <- model.response(mframe)
#   }else {
#     mframe <- model.frame(formula, data1, na.action = NULL)
#     y <- model.response(mframe)
#     y[is.na(y)] <- 0
#   }
#   if (any(sapply(lapply(as.data.frame(y), table), length) == 1)) {
#     y <- y[, !(sapply(apply(y, 2, table), length) == 1)]
#     cat("\n ALERT: at least one manifest variable contained only one\n outcome category, and has been r
# emoved from the analysis. \n\n")
#   }
#   x <- model.matrix(formula, mframe)
#   
#   #indicate func
#   indicate <- function(x,y) ifelse(x==y,1,0)
#   
#   diff = 100000
#   iter = 1
#   Q = 100000
#   llik = c()
#   llik[1] <- -Inf
#   Q.list = c()
# 
#   ###### initial value #######
#   L = nclass
#   M = ncol(y)
#   r.m <- t(matrix(apply(y, 2, max)))
#   r = rep(1/nclass,nclass) # class prob
#   n = nrow(y)
#   
#   ## initial prob 
#   if (is.null(initial_prob)) {
#     probs = list()
#     for (m in 1:M) {
#       probs[[m]] <- matrix(runif(L * r.m[m]),
#                            nrow = L, ncol = r.m[m])
#       probs[[m]] <- probs[[m]]/rowSums(probs[[m]])
#     }
#     probs.new <- probs
#   }
#   else {
#     probs.new <- probs <- initial_prob
#   }
# 
#   ## initial zhat = p(zi=l |Y,r,rhomk)
#   z0 <- matrix(1,n,L) # z0 : numerator of zhat
#   
#   for(i in 1:n){
#     for(l in 1:L){
#       for(m in 1:M) {
#         for(k in 1:r.m[m]) {
#           z0[i,l] = z0[i,l]*(probs[[m]][l,k]^indicate(y[i,m],k))
#         } 
#       }
#       z0[i,l] = r[l] * z0[i,l]
#     }
#   }
#   zhat <- z0 * (1/rowSums(z0))
# 
#   ### EM algorithm
#   while (diff > tol & iter < maxiter) {
#     
#     ### E- step: compute 
#     
#     # Q 정의
#     Q.new = sum(zhat * log(z0))
#     
#     ### M = step : update r, p11, p21, p12, p22
#     
#     #r.new
#     r.new = colSums(zhat) / sum(zhat)
#     # r.new = colSums(zhat) / sum(zhat)
#     
#     
#     # #prob.new(old version)
#     # for (l in 1:L) {
#     #   for (m in 1:M) {
#     #     for(k in 1:r.m[m]) {
#     #       num = 0
#     #       for(i in 1:n){num = num + zhat[i,l]*indicate(y[i,m],k)}
#     #       probs.new[[m]][l,k] = round(num / sum(zhat[,l]),4)
#     #       
#     #     }
#     #     colnames(probs.new[[m]])<- 1:r.m[m]
#     #     rownames(probs.new[[m]]) <- paste('class',1:L,":", sep = ' ')
#     #   }
#     # }
#     # names(probs.new) <- colnames(y)
#   
#     #prob.new (new version)
#     for (l in 1:L) {
#       for (m in 1:M) {
#         for(k in 1:r.m[m]) {
#           probs.new[[m]][l,k] = sum(zhat[,l]*indicate(y[,m],k)) / sum(zhat[,l])
#         }
#       }
#     }
#     
#     
#   
#     # new zhat
#     z0.new <- matrix(1,n,L) # z0.new : numerator of zhat.new
#     for(i in 1:n){
#       for(l in 1:L){
#         for(m in 1:M) {
#           for(k in 1:r.m[m]) {
#             z0.new[i,l] = z0.new[i,l]*(probs.new[[m]][l,k]^indicate(y[i,m],k))
#           } 
#         }
#         z0.new[i,l] = r.new[l] * z0.new[i,l]
#       }
#     }
#     zhat.new <- z0.new * (1/rowSums(z0.new))
#     
#     
#     # calculate observed loglikelihood to check convergence
#     llik[iter+1] = sum(log(rowSums(z0.new)))
#     diff = llik[iter+1] - llik[iter]
#     
#     
#     # update para
#     r = r.new
#     probs = probs.new
#     zhat = zhat.new
#     z0 = z0.new
#     Q = Q.new
#     
#     
#     cat("Iter", iter, ": r = ",r.new, ", llik = ",llik[iter+1],", diff = ", diff, '\n')
# 
#     Q.list = c(Q.list, Q)
#     iter=iter+1
# 
#   }
#   ret = list()
#   ret$nclass <- L
#   ret$npara <- (L * sum(r.m - 1)) + (L - 1)
#   ret$aic <- (-2 * tail(llik,1)) + (2 * ret$npar)
#   ret$bic <- (-2 * tail(llik,1)) + (log(n) * ret$npar)
#   ret$caic <- (-2 * tail(llik,1)) + ret$npar * (log(n) + 1)
#   
#   # prob
#   for(m in 1:M){
#     colnames(probs[[m]])<- 1:r.m[m]
#     rownames(probs[[m]]) <- paste('class',1:L,":", sep = ' ')
#   }
#   names(probs.new) <- colnames(y)
#   
#   return(list(r=r.new, ollik = llik, probs = probs, Q = Q.list, z= zhat,
#               ret = ret))
# }



