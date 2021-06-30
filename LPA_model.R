LPA_model = function(num_data, initial_lambda=NULL ,nclass, 
                      maxiter = 2000, tol = 1e-5, na.rm=TRUE){
  # mframe <- model.frame(formula, cate_data, na.action = NULL)
  # data1 <- cate_data[rowSums(is.na(model.matrix(formula, mframe))) == 0, ]
  # if (na.rm) {
  #   mframe <- model.frame(formula, data1)
  #   y <- model.response(mframe)
  # }else {
  #   mframe <- model.frame(formula, data1, na.action = NULL)
  #   y <- model.response(mframe)
  #   y[is.na(y)] <- 0
  # }
  # if (any(sapply(lapply(as.data.frame(y), table), length) == 1)) {
  #   y <- y[, !(sapply(apply(y, 2, table), length) == 1)]
  #   cat("\n ALERT: at least one manifest variable contained only one\n outcome category, and has been r
  #       emoved from the analysis. \n\n")
  # }
  # x <- model.matrix(formula, mframe)
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
  # M = ncol(y)
  d2 = ncol(num_data)
  # r.m <- t(matrix(apply(y, 2, max)))
  r = rep(1/nclass,nclass) # class prob
  n = nrow(num_data)
  
  # ## initial prob 
  # if (is.null(initial_prob)) {
  #   probs = list()
  #   for (m in 1:M) {
  #     probs[[m]] <- matrix(runif(L * r.m[m]),
  #                          nrow = L, ncol = r.m[m])
  #     probs[[m]] <- probs[[m]]/rowSums(probs[[m]])
  #   }
  #   probs.new <- probs
  # }
  # else {
  #   probs.new <- probs <- initial_prob
  # }
  
  ## initial lambda
  sam.num <- apply(num_data,2,mean) %>% max() %>%  round(0)
  if (is.null(initial_lambda)) {
    lambda = list()
    for (l in 1:L) {
      lambda[[l]] <- sample(1:sam.num, d2, replace = T)  # sam.num (can change num)
    }
    lambda.new <- lambda
  }
  else {
    lambda.new <- lambda <- initial_lambda
  }
  
  ### change y into dummy
  # # make dummy fun
  # indi = function(x,y) {
  #   mat <- rep(0,y)
  #   mat[x] <- 1
  #   return(mat)
  # }
  
  # # y_dummy
  # yi <- list()
  # flag = 0
  # for (k in r.m) {
  #   flag = flag + 1
  #   yi[[flag]]<- sapply(y[,flag],function(x) indi(x,k)) %>%  t
  # }
  
  
  # make z0 matrix
  find.z0 <- function(r,lambda){
    
    ## 01. cate part
    # ### calculate sum(I(yi)*probs)
    # yi.prob <- list()
    # # tprobs <- sapply(probs,t)
    # tprobs <- lapply(probs,t)
    # 
    # # mapply() : calculate fun to each element in two list
    # yi.prob = mapply(FUN=function(x,y) x%*%y, x = yi, y = tprobs, SIMPLIFY = FALSE) 
    # xy <- Reduce("*", yi.prob) # Reduce() : Calculate values corresponding to each element in list
    # 
    ## 02. numarical part
    temp <- matrix(1,n,L) # z0 : numerator of zhat
    
    # for(i in 1:n){
    #   for(l in 1:L){
    #     for(t in 1:d2) {
    #       temp[i,l] = temp[i,l]*(exp(-lambda[[l]][t])*(lambda[[l]][t]^num_data[i,t])/factorial(num_data[i,t]))
    #     }
    #   }
    # }
    for(l in 1:L) {
      temp1 = sweep(num_data, MARGIN = 2, log(lambda[[l]]),`*`) # sweep : multiply matrix by vector 
      temp2 <- sweep(temp1, MARGIN = 2, lambda[[l]], `-`) - log(apply(num_data,2, factorial))
      temp[,l] = exp(rowSums(temp2))
    }
    
    ### calculate z0
    # z00 <- matrix(rep(r, nrow(xy)),ncol=L, byrow = T) * xy * temp
    z00 <- sweep(temp, MARGIN = 2, r, "*")
    return(z00)
  }
  
  z0 <- find.z0(r, lambda)
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
    # probs.new <- lapply(yi, function(x) (t(zhat) %*% x) * matrix(rep(1/apply(zhat,2, sum),ncol(x)),ncol=ncol(x)))
    
    
    ## lambda.new
    # for (l in 1:L) {
    #   for (t in 1:d2) {
    #     lambda.new[[l]][t] = sum(zhat[,l]*num_data[,t]) / sum(zhat[,l])   
    #   }
    # }
    
    lambda.new <- t(zhat) %*% num_data * matrix(rep(1/apply(zhat,2, sum),ncol(num_data)),ncol=ncol(num_data))
    lambda.new <- as.list(as.data.frame(t(lambda.new)))
    
    z0.new <- find.z0(r.new, lambda.new)
    zhat.new <- z0.new * (1/rowSums(z0.new))
    
    
    # calculate observed loglikelihood to check convergence
    llik[iter+1] = sum(log(rowSums(z0.new)))
    diff = llik[iter+1] - llik[iter]
    
    
    # update para
    r = r.new
    # probs = probs.new
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
  # ret$npara <- (L * sum(r.m - 1)) + (L - 1)
  ret$aic <- (-2 * tail(llik,1)) + (2 * ret$npar)
  ret$bic <- (-2 * tail(llik,1)) + (log(n) * ret$npar)
  ret$caic <- (-2 * tail(llik,1)) + ret$npar * (log(n) + 1)
  
  # # prob
  # for(m in 1:M){
  #   colnames(probs[[m]])<- 1:r.m[m]
  #   rownames(probs[[m]]) <- paste('class',1:L,":", sep = ' ')
  # }
  # names(probs.new) <- colnames(y)
  # 
  # lambda
  for(l in 1:L){
    names(lambda[[l]])<- 1:d2
    names(lambda)[l] <- paste('class',l)
  }
  
  return(list(r=r.new, ollik = llik, lambda = lambda, Q = Q.list, z= zhat,
              ret = ret))
}

new <- LPA_model(num_data, nclass=3, initial_lambda=NULL ,
                 maxiter = 2000, tol = 1e-5, na.rm=TRUE)
apply(new$z,1,max)
colnames(new$z)[apply(new$z, 1, which.max)]

head(iris)
