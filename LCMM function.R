library(dplyr)
library(ggplot2)

###################################################################
#  00. LCMM dataset
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
# class 4
y_prob = list(c1p1 = c(0.6,0.3,0.1),
              c2p1 = c(0.2,0.4,0.4),
              c3p1 = c(0.2,0.5,0.3),
              c4p1 = c(0.6,0.2,0.2),
              c1p2 = c(0.3,0.7),
              c2p2 = c(0.5,0.5),
              c3p2 = c(0.3,0.7),
              c4p2 = c(0.6,0.4),
              c1p3 = c(0.2,0.6,0.1,0.1),
              c2p3 = c(0.25,0.25,0.25,0.25),
              c3p3 = c(0.1,0.1,0.6,0.2),
              c4p3 = c(0.2,0.4,0.3,0.1))


## # of categories each observed variable has
mnum = c()
for (i in 1:length(y_prob)) {
  mnum <- c(mnum, length(y_prob[[i]]))
}
mnum

## class prob
r = t(c(0.3,0.6,0.1)) #class 3
r = t(c(0.4,0.1,0.5)) #class 3
r = t(c(0.1,0.3,0.5,0.1)) #class 4

L=length(r)
M=length(y_prob) / L
n=1000

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
# class 3
lambda = list(c1l1 = c(5,1,20),
              c2l1 = c(1.5,4,2),
              c3l1 = c(8,20,2))

# class 3
lambda = list(c1l1 = c(5,1,0,0,20,0,0,1),
              c2l1 = c(1.5,0,0,0,2,8,0,2),
              c3l1 = c(0,0,0,5,0,0,20,2))

# class 3
lambda = list(c1l1 = c(5,1,3,20,3),
              c2l1 = c(1.5,4,10,2,10),
              c3l1 = c(8,20,16,2,2))

# class 4
lambda = list(c1l1 = c(5,1,20),
              c2l1 = c(15,4,1),
              c3l1 = c(8,20,16),
              c4l1 = c(2,2,2))

# class 4
lambda = list(c1l1 = c(5,1,3,20,3),
              c2l1 = c(1.5,4,10,2,10),
              c3l1 = c(8,20,16,2,2),
              c4l1 = c(2,2,2,2,2))

d2=length(lambda[[1]])


## each class dataset for numerical variables
y2.dat <- vector(mode="list", length = d2)
temp <- matrix(0, nrow=n, ncol=d2)
for (l in 1:L){
  for (t in 1:d2){
    temp[,t] <- rpois(n,lambda[[l]][t])
  }
  y2.dat[[l]] <- temp
  colnames(y2.dat[[l]]) <-  paste0('num',1:d2)
}

#check-up
# > y2.dat[[1]] %>% colSums() / n
# num1   num2   num3   num4   num5 
# 4.986  0.960  3.004 20.068  2.978 
# > y2.dat[[2]] %>% colSums() / n
# num1  num2  num3  num4  num5 
# 1.449 3.961 9.864 2.106 9.989 
# > y2.dat[[3]] %>% colSums() / n
# num1   num2   num3   num4   num5 
# 8.033 20.129 16.080  2.056  2.003 

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
dim(dat) #1000 * 9


######################################
### 01. my function : LCMM
formula = cbind(cate1, cate2, cate3, cate4) ~1  # class3
formula = cbind(cate1, cate2, cate3) ~1         # class4

cate_data = dat[,1:M] 
num_data = dat[,-c(1:M)]



initial_prob=NULL;initial_lambda=NULL ;nclass=3; 
maxiter = 2000; tol = 1e-5; na.rm=TRUE

LCMM_model = function(formula, cate_data, num_data, initial_prob=NULL,initial_lambda=NULL ,nclass, 
                      maxiter = 2000, tol = 1e-5, na.rm=TRUE){
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
  sam.num <- apply(num_data,2,mean) %>% max() %>%  round(0)
  if (is.null(initial_lambda)) {
    lambda = list()
    for (l in 1:L) {
     ##origin
     lambda[[l]] <- sample(1:sam.num, d2, replace = T)  # sam.num (can change num)
    }
    lambda.new <- lambda
  }else {
    lambda.new <- lambda <- initial_lambda
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
  find.z0 <- function(r, probs, lambda){
    
    ## 01. cate part
    ### calculate sum(I(yi)*probs)
    yi.prob <- list()
    # tprobs <- sapply(probs,t)
    tprobs <- lapply(probs,t)
    
    # mapply() : calculate fun to each element in two list
    yi.prob = mapply(FUN=function(x,y) x%*%y, x = yi, y = tprobs, SIMPLIFY = FALSE) 
    xy <- Reduce("*", yi.prob) # Reduce() : Calculate values corresponding to each element in list
    
    ## 02. numarical part
    temp <- matrix(1,n,L) # z0 : numerator of zhat
    
    # for(i in 1:n){
    #   for(l in 1:L){
    #     for(t in 1:d2) {
    #       temp[i,l] = temp[i,l]*(exp(-lambda[[l]][t])*(lambda[[l]][t]^num_data[i,t])/factorial(num_data[i,t]))
    #     }
    #   }
    # }
    ## origin
    # for(l in 1:L) {
    #   temp1 = sweep(num_data, MARGIN = 2, log(lambda[[l]]),`*`) # sweep : multiply matrix by vector
    #   temp2 <- sweep(temp1, MARGIN = 2, lambda[[l]], `-`) - log(apply(num_data,2, factorial))
    #   temp[,l] = exp(rowSums(temp2))
    # }
    
    ## new (06/09)
    for(l in 1:L){
      lamb = matrix(rep(lambda[[l]], nrow(num_data)), nrow=nrow(num_data), byrow = T)
      temp1 = lamb ^ num_data
      temp2 = sweep(temp1, MARGIN = 2,exp(-lambda[[l]]), "*" )
      fac_inv = 1/factorial(num_data)
      temp3 = temp2  * fac_inv
      temp[,l] = apply(temp3,1,prod)
    }
    
    ### calculate z0
    # z00 <- matrix(rep(r, nrow(xy)),ncol=L, byrow = T) * xy * temp
    z00 <- sweep(xy * temp, MARGIN = 2, r, "*")
    return(z00)
  }
  
  z0 <- find.z0(r, probs, lambda)
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

    
    ## lambda.new
    # for (l in 1:L) {
    #   for (t in 1:d2) {
    #     lambda.new[[l]][t] = sum(zhat[,l]*num_data[,t]) / sum(zhat[,l])   
    #   }
    # }
    
    lambda.new <- t(zhat) %*% num_data * matrix(rep(1/apply(zhat,2, sum),ncol(num_data)),ncol=ncol(num_data))
    lambda.new <- as.list(as.data.frame(t(lambda.new)))
    
    z0.new <- find.z0(r.new, probs.new, lambda.new)
    zhat.new <- z0.new * (1/rowSums(z0.new))
    
    
    # calculate observed loglikelihood to check convergence
    llik[iter+1] = sum(log(rowSums(z0.new)))
    diff = llik[iter+1] - llik[iter]
    
    
    # update para
    r = r.new
    probs = probs.new
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
  
  return(list(r=r.new, ollik = llik, probs = probs,lambda = lambda, Q = Q.list, z= zhat,
              ret = ret))
}

lcmm <- LCMM_model(formula, cate_data, num_data, initial_prob=NULL,initial_lambda=NULL,
                   nclass = 3, maxiter = 2000, tol = 1e-5, na.rm=TRUE)
lcmm <- LCMM_model(formula, cate_data, num_data, initial_prob=NULL,initial_lambda=NULL,
                   nclass = 4, maxiter = 2000, tol = 1e-5, na.rm=TRUE)

lcmm$ollik %>%  tail(1)
lcmm$r
lcmm$probs
lcmm$lambda

##### iterate simulation for culcurate average , MSE 
###### total code(average probs)
lcmm.list <- list()
probs.list <- list()
lambda.list <- list()
r.list <- matrix(0,100,3)

for (i in 1:100){
  lcmm.list[[i]] <-  LCMM_model(formula, cate_data, num_data, initial_prob=NULL,initial_lambda=NULL,
                                nclass = 3, maxiter = 2000, tol = 1e-5, na.rm=TRUE)

  r.order <- order(lcmm.list[[i]]$r)
  #r.order
  lcmm.list[[i]]$r[r.order]
  
  probs.list[[i]]<- lapply(lcmm.list[[i]]$probs, FUN = function(x) x[r.order,])
  r.list[i,]<- lcmm.list[[i]]$r[r.order]
  lambda.list[[i]]<- lcmm.list[[i]]$lambda[r.order]
}

### bowplot for probs, r, lambda
lapply(probs.list, FUN = function(x) x[[1]][1,1]) %>%  unlist

box.probs <- matrix(0,nrow=100,ncol=12)

for (L in 1:3){
  for (i in 1:2) {
    box.probs[,i*L] <- lapply(probs.list, FUN = function(x) x[[4]][L,i]) %>%  unlist
  }
}

boxplot(box.probs)

#### total average probs
# average 100 times probs
sum.probs <- probs.list[[1]]
for (j in 2:100) {
  sum.probs <- mapply(FUN=function(x,y) x+y, x = sum.probs, y = probs.list[[j]], SIMPLIFY = FALSE)
}


divide.list <- sum.probs
divide.list[[1]] <- matrix(1/100, nrow=3, ncol=3)
divide.list[[2]] <- matrix(1/100, nrow=3, ncol=2)
divide.list[[3]] <- matrix(1/100, nrow=3, ncol=3)
divide.list[[4]] <- matrix(1/100, nrow=3, ncol=4)
ave.probs <- mapply(FUN=function(x,y) x*y, x = sum.probs, y = divide.list, SIMPLIFY = FALSE)
ave.probs


#### MSE
## 1) true probs
y.probs <- unlist(y_prob)
initial.probs <- list()
initial.probs[[1]] <- matrix(y.probs[1:9], nrow = 3, ncol =3, byrow = T)
initial.probs[[2]] <- matrix(y.probs[10:15], nrow = 3, ncol =2, byrow = T)
initial.probs[[3]] <- matrix(y.probs[16:24], nrow = 3, ncol =3, byrow = T)
initial.probs[[4]] <- matrix(y.probs[25:36], nrow = 3, ncol =4, byrow = T)
init.probs <- lapply(initial.probs, FUN=function(x) x[c(3,1,2),]) 
init.probs


## 2) error probs
p.error2.list <- list()
for (k in 1:100){
  p.error2.list[[k]] <- mapply(FUN=function(x,y) x-y, x=init.probs, y=probs.list[[k]])
  p.error2.list[[k]] <- mapply(FUN=function(x,y) x*y, x=p.error2.list[[k]], y=p.error2.list[[k]])
}


p.sse.list <- p.error2.list[[1]]
for (h in 2:100){
  p.sse.list <- mapply(FUN=function(x,y) x+y, x = p.sse.list, y = p.error2.list[[h]], SIMPLIFY = FALSE)
}
## 3) mse probs
mse.probs <- mapply(FUN=function(x,y) x*y, x = p.sse.list, y = divide.list, SIMPLIFY = FALSE)

hundred.list <- divide.list
hundred.list[[1]] <- matrix(100, nrow=3, ncol=3)
hundred.list[[2]] <- matrix(100, nrow=3, ncol=2)
hundred.list[[3]] <- matrix(100, nrow=3, ncol=3)
hundred.list[[4]] <- matrix(100, nrow=3, ncol=4)

mse100.probs <- mapply(FUN=function(x,y) x*y, x = mse.probs, y = hundred.list, SIMPLIFY = FALSE)
mse100.probs

###################################
## total average r
###### total code
r.list

r[order(r)]

# average r
apply(r.list, 2,mean)

# mse r
error.r <- sweep(r.list, MARGIN = 2, r[order(r)],`-`)
sse.r <- error.r * error.r

mse100.r <-apply(sse.r, 2, mean) * 100
mse100.r


###################################
## total average lambda
lambda.list
# average 20 times lambda
sum.lambda <- lambda.list[[1]]
for (j in 2:100) {
  sum.lambda <- mapply(FUN=function(x,y) x+y, x = sum.lambda, y = lambda.list[[j]], SIMPLIFY = FALSE)
}


divide.list <- sum.lambda
divide.list[[1]] <- rep(1/100,3)
divide.list[[2]] <- rep(1/100,3)
divide.list[[3]] <- rep(1/100,3)
ave.lambda <- mapply(FUN=function(x,y) x*y, x = sum.lambda, y = divide.list, SIMPLIFY = FALSE)
ave.lambda


#### MSE
## 1) true lambda
init.lambda <-lambda[order(r)]

## 2) error lambda
p.error2.list <- list()
for (k in 1:100){
  p.error2.list[[k]] <- unlist(init.lambda)- unlist(lambda.list[[k]])
  p.error2.list[[k]] <- p.error2.list[[k]] * p.error2.list[[k]]
}

p.sse.list <- Reduce("+", p.error2.list)
# p.sse.list <- p.error2.list[[1]]
# for (h in 2:20){
#   p.sse.list <- p.sse.list + p.error2.list[[h]]
# }
## 3) mse lambda
mse.lambda <- p.sse.list / 100

mse100.lambda <- mse.lambda *100
mse100.lambda 

matrix(mse100.lambda, nrow=3)
matrix(unlist(init.lambda), nrow=3)



#########################################
####### etc
#########################################
library(tidyverse)
library(dplyr)
survey <- read.csv("https://raw.githubusercontent.com/whipson/tidytuesday/master/young_people.csv") %>% 
  dplyr::select(History:Pets)
head(survey)
names(survey)

# install.packages("tidyLPA")
library(tidyLPA)
pisaUSA15[1:100, ] %>%
  dplyr::select(broad_interest, enjoyment, self_efficacy) %>%
  single_imputation() %>%
  estimate_profiles(3)

head(pisaUSA15)


pisaUSA15[1:100, ] %>%
  dplyr::select(broad_interest, enjoyment, self_efficacy) %>%
  single_imputation() %>%
  scale() %>%
  estimate_profiles(3) %>% 
  plot_profiles()

iris_subset <- sample(nrow(iris), 20) # so examples execute quickly
results <- iris %>%
  subset(select = c("Sepal.Length", "Sepal.Width",
                    "Petal.Length", "Petal.Width")) %>%
  estimate_profiles(1:3) %>%
  compare_solutions()


iris %>%
  subset(select = c("Sepal.Length", "Sepal.Width",
                    "Petal.Length")) %>%
  estimate_profiles(n_profiles = 1:4)


iris_sample <- iris[c(1:4, 51:54, 101:104), ] 

iris_sample %>%
  subset(select = c("Sepal.Length", "Sepal.Width","Petal.Length")) %>%
  estimate_profiles(n_profiles = 1:4, variances = c("equal", "varying"),
                    covariances = c("zero", "zero"))

if(interactive()){
  results <- iris %>%
    dplyr::select(Sepal.Length, Sepal.Width, Petal.Length, Petal.Width) %>%
    estimate_profiles(3)
  get_estimates(results)
  get_estimates(results[[1]])
}

dim(iris)

data(id_edu)
head(id_edu)
dim(id_edu)
aa <-  id_edu[,1:15] %>%  estimate_profiles(5)
get_estimates(aa)
head(id_edu)

aa <- dat %>% estimate_profiles(3)
get_estimates(aa) %>%  View
names(get_estimates(aa))

head(dat)


save(LCMM_model, "LCMM_model.r") # 객체저장
save.image(file="LCMM_model.RData") #작업공간 전체저장
load(file="LCMM_model.RData")


###########################
## Newton Raphson Method in R
###########################

newton <- function(fun, tol = 1e-7, x0=1, N=300) {
  h <- 1e-7
  i <- 1
  x1 <- x0
  p <- numeric(N)
  while(i <= N){
    df.dx <- (fun(x0 + h) - fun(x0))/h
    x1 <- x0 - (fun(x0)/df.dx)
    p[i] <- x1
    i = i+1
    if(abs(x1-x0) < tol) break
    x0 = x1
  }
  return(p[1:(i-1)])
}

fun <- function(x){
  x^3 + 2
}

newton(fun, x0=2)

## visualization























































