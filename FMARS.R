library(dplyr)
library(pracma)
library(frechet)
library(Rcpp)
library(parallel)
library(spherepc)


mars_train_c <- function(y,x,degree=2,mx,sample_num,type)
{
  if(sample_num>0){
    fwd_out = forward_c(y, x,degree,mx,sample_num,type)
    #bwd_out = backward(fwd_out,y,x)
    result <- fwd_out$splits
  }
  else{
    fwd_out = forward_zero(y, x,degree,mx,type)
    result <- fwd_out$splits
  }
  
  return(result)
}

forward_zero <- function(y,x,degree,mx,type)
{
  n <- nrow(x)  #样本个数
  p <- ncol(x)  #样本维数
  B <- as.matrix(rep(1,n),nrow=n)
  splits = list(data.frame(
    m = 0,
    v = 0
  ))
  M <- 1
  mx <- mx%/%2+1
  
  while (M < mx) {
    lof_best = Inf
    split_best = c(m = 0,v = 0)
    for (m in 1:M) {   #B中以有的元素，乘以新的分裂点
      #print(m)
      #保证阶数
      splits_loc <- splits[[m]]
      if (dim(splits_loc)[1]>degree)
      {
        next
      }
      for (v in c(1:p)) {      #对于x的每一列分别寻找
        # Select variable to split on
        #print(tt)
        Bnew =data.frame(B, Btem1 = B[, m]* x[,v])
        B_test <- as.matrix(Bnew[-1])
        
        #检测矩阵是否可逆
        fit <- try(solve(var(B_test)),silent = TRUE)
        if ('try-error' %in% class(fit)){
          next}
        
        if(det(var(B_test))<0.00000000000000000000000000000000000000001){
          next
        }
        lof1 = loss(y,Bnew,type)  #把全是1的那一项删除
        lof <- mean(lof1)
        
        if (lof < lof_best) {
          lof_best = lof
          split_best = c(m = m,v = v)
        } 
      } 
    } 
    #rint(M)
    m = split_best[['m']]
    v = split_best[['v']]
    B = cbind(B, B[, m] * x[,v])
    right = rbind(splits[[m]], c(m, v))   #选择和第m个相乘，就在第m个基础上加上这个的信息
    splits = c(splits,  list(right))   #更新splits
    #print(splits)
    M = M + 1
    #print(splits)
  } 
  #colnames(B) = paste0('B', seq(0, control$Mmax))
  splits[[1]] <- NULL   #删除全为1的项
  B <- B[,-1]
  #print(B)#看一下输出结果
  #print(det(var(B)))
  
  #print(111111111)
  #print(det(var(B_test))
  return(list(B = B, splits = splits))
  
}


predict_loss_distribution <- function(apply_data,y,B,Pmat,Amat,b0)
{
  p <- length(apply_data)
  p1 <- p-20
  p2 <- p-21
  y_true <- apply_data[p1:p]
  xout <- apply_data[1:p2]
  n <- dim(y)[1]  #样本个数
  #dime <- dim(y)[2]  #矩阵维数
  B <- as.matrix(B)
  xbar <- colMeans(B)
  #Sigma <- cov(xin) * (n - 1)/n
  Sigma <- cov(B)
  invSigma <- solve(Sigma)
  #trapz需要的参数
  nqSup <- 21
  qSup <- seq(0,1,length.out = nqSup)
  xp <- c(qSup, qSup[nqSup:1])
  n_sup <- 2 * nqSup
  s <- 1 + t(t(B) - xbar) %*% invSigma %*% (xout - xbar)
  s <- as.vector(s)
  gx <- colMeans(y * s)
  res <- do.call(osqp::solve_osqp, list(P = Pmat, q = -gx, A = Amat, l = b0, pars = osqp::osqpSettings(verbose = FALSE)))
  y_hat <- res$x
  dist <- pracma::trapz(qSup, (y_true - y_hat)^2)
  return(dist)
}

predict_loss_sphere <- function(index,y_trans,B,y,function1)
{
  n = nrow(B)
  p = ncol(B)
  xout <- B[index,]
  y_true <- y[index,]
  invVa = solve(var(B))
  mx = apply(B, 2, mean)
  s <- 1 + t(t(B) - mx) %*% invVa %*% (xout - mx)
  s <- as.vector(s)
  y_hat <- spherepc::ExtrinsicMean(data=y_trans,weights=s)
  y_pred <- spherepc::Trans.Euclid(y_hat)
  cross <- crossprod(y_pred,y_true)[1,1]
  dist <- acos(cross)^2
  
  return(dist)
}

predict_loss_matrix <- function(index,y,B,L,D)
{
  n = nrow(B)
  p = ncol(B)
  xout <- B[index,]
  y_true <- y[index,,]
  invVa = solve(var(B))
  mx = apply(B, 2, mean)
  ss = 0
  U = 0
  E = 0
  s <- 1 + t(t(B) - mx) %*% invVa %*% (xout - mx)
  s <- as.vector(s)
  for (i in 1:n) {
    ss = ss + s[i]
    U = U + s[i] * L[[i]]
    E = E + s[i] * log(D[[i]])}
  SS = U/ss + diag(exp(E/ss))
  y_hat = t(SS) %*% SS
  
  #如果不是对称正定的，转换成对称正定的
  fit <- try(chol(y_hat),silent = TRUE)
  if ('try-error' %in% class(fit)){
    dist <- Inf
  }
  else{
    p <- dim(y_hat)[1]
    M.a <- chol(y_hat ,pivot = TRUE)
    D.a <- diag(M.a)
    L.a <- M.a - diag(D.a)
    M.b <- chol(y_true ,pivot = TRUE)
    D.b <- diag(M.b)
    L.b <- M.b - diag(D.b)
    dist <- sqrt(sum((L.a - L.b)^2) + sum((log(D.a) - log(D.b))^2))
    
  }
  
  return(dist)
}


loss <- function(y,B1,type)
{
  B <- B1[,-1]
  fit <- try(solve(var(B)),silent = TRUE)
  if ('try-error' %in% class(fit)){
    #print(B_test)
    lof=c(100000000000000)
  }
  else
  {
    if(type == 'distribution'){
      nqSup <- 21
      qSup <- seq(0,1,length.out = nqSup)
      A <- cbind(diag(nqSup), rep(0, nqSup)) + cbind(rep(0, nqSup), -diag(nqSup))
      A <- A[, -c(1, ncol(A))]
      b0 <- rep(0, nqSup - 1)
      Pmat <- as(diag(nqSup), "sparseMatrix")
      Amat <- as(t(A), "sparseMatrix")
      
      apply_data <- cbind(B,y)
      lof <- parApply(cl=cl,apply_data,1,FUN=predict_loss_distribution,y=y,B=B,Pmat=Pmat,Amat=Amat,b0=b0)
    }
    else if(type == 'matrix'){
      n <- dim(y)[1]  
      dime <- dim(y)[2]  
      lof <- 0
      B <- as.matrix(B)
      #print(B)
      index_list <- matrix(c(1:n),nrow=n)
      MM = list()
      for (i in 1:n) {
        MM[[i]] = y[i, ,]
      }
      LL = lapply(MM, chol)   
      L = lapply(LL, function(X) X - diag(diag(X)))
      D = lapply(LL, function(X) diag(X))
      
      lof <- parApply(cl=cl,index_list,1,FUN=predict_loss_matrix,y=y,B=B,L=L,D=D)
    }
    else if(type =='sphere'){
      n <- dim(y)[1] 
      dime <- dim(y)[2]  
      lof <- 0
      B <- as.matrix(B)
      index_list <- matrix(c(1:n),nrow=n)
      S2_obs <- matrix(0,n,2)
      for(i in c(1:n)){
        S2_obs[i,] = Trans.sph(y[i,])
      }
      lof <- parApply(cl=cl,index_list,1,FUN=predict_loss_sphere,y_trans=S2_obs,B=B,y=y)
    }
  }
  
  return (lof)
}

backward = function(fwd_out,y,x) 
{
  B = fwd_out$B
  splits = fwd_out$splits
  splits[[1]] <- NULL   
  
  lof <- gcv_c(y,B,splits)
  B <- B[,2:ncol(B)]
  mx <- ncol(B)

  lof_list <- c()  
  B_list <- list(B)
  splits_list <- list(splits)
  
  B_loc <- B
  splits_loc <- splits
  
  for (M in c(mx:2)){   
    lof_best=Inf
    for (m in c(1:M)){
      B_m <- B_loc[,-m]
      splits_m <- splits_loc[-m]
      lof_m <- gcv_c(y,B_loc,splits_m)
      
      if (lof<lof_best){
        lof_best <- lof_m
        B_best <- B_m
        splits_best <- splits_m
      }
      
    }
    
    lof_list <-  c(lof_list, lof_best) 
    B_list <-  c(B_list, list(B_best))  
    splits_list <-  c(splits_list, list(splits_best))
    
    B_loc <- B_best 
    splits_loc <- splits_best
    
  }
  best_num <- which.min(lof_list)
  best_splits <- splits_list[[best_num]]
  
  return(best_splits)
}

gcv_c <- function(y,B,splits)
{
  lof <- loss_c(y,B)
  lof1 <- mean(lof)
  c <- 3
  r <- ncol(B)-1       
  k <- count_knots_c(splits) 
  n <- nrow(y)
  
  par_num <- r+c*k
  result <- lof1/(1-par_num/n)^2
  
  return(result)
  
}

count_knots_c <- function(splits)
{
  vt <- data.frame(m = 0,v = 0,s = NA,t = NA)
  for (i in c(1:length(splits))){
    vt <- cbind(vt,splits[[i]])
    
  }
  vt1 <- t(vt)
  vt2 <- as.data.frame(vt1)
  vt3 <- distinct(vt2,3,4)
  num <- dim(vt3)[1]
  
  return(num)
}

mars_predict_list <- function(xout,splits,M,xin,type)
{
  k <- nrow(xin)
  n <- nrow(xout)
  p <- ncol(xout)
  m <- ncol(M)
  
  
  if(type =='sphere' ){
    result <-  array(0,c(n,m))
    S2_obs <- matrix(0,k,2)
    for(i in c(1:k)){
      S2_obs[i,] = Trans.sph(M[i,])
    }
    
    for(i in c(1:n)){
      x <- matrix(xout[i,],ncol=p)
      predict_result <- mars_predict(x,mars.fit,S2_obs,xin,type)
      result[i,] <- predict_result
    }
  }
  else if(type == 'distribution'){
    result <-  array(0,c(n,m))
    for(i in c(1:n)){
      x <- matrix(xout[i,],ncol=p)
      predict_result <- mars_predict(x,mars.fit,M,xin,type)
      result[i,] <- predict_result}
  }
  else if(type == 'matrix'){
    result <-  array(0,c(n,m,m))
    for(i in c(1:n)){
      x <- matrix(xout[i,],ncol=p)
      predict_result <- mars_predict(x,mars.fit,M,xin,type)
      result[i,,] <- predict_result}
  }
  
  return(result)
}

mars_predict <- function(xout,splits,y,xin,type)
{
  if(sample_num>0){
    Bin <- trans_x2B(xin,splits)
    #print(Bin)
    Bout <- trans_x2B(xout,splits)
  }
  else{
    Bin <- trans_x2B_zero(xin,splits)
    #print(Bin)
    Bout <- trans_x2B_zero(xout,splits)
  }
  xin <- Bin
  xout <- Bout
  
  if(type == 'distribution'){
    qin <- y
    k <- nrow(xout)
    n <- nrow(xin)
    m <- ncol(qin)
    xbar <- colMeans(xin)
    Sigma <- cov(xin) * (n - 1)/n
    invSigma <- solve(Sigma)
    A <- cbind(diag(m), rep(0, m)) + cbind(rep(0, m), -diag(m))
    A <- A[, -c(1, ncol(A))]
    b0 <- rep(0, m - 1)
    Pmat <- as(diag(m), "sparseMatrix")
    Amat <- as(t(A), "sparseMatrix")
    s <- 1 + t(t(xin) - xbar) %*% invSigma %*% (xout - xbar)
    s <- as.vector(s)
    gx <- colMeans(qin * s)
    res <- do.call(osqp::solve_osqp, list(P = Pmat, q = -gx, A = Amat, l = b0, pars = osqp::osqpSettings(verbose = FALSE)))
    y_hat <- res$x
  }
  
  else if(type == 'matrix'){
    n = nrow(Bin)
    p = ncol(Bin)
    xout <- Bout
    invVa = solve(var(Bin))
    mx = apply(Bin, 2, mean)
    MM = list()
    for (i in 1:n) {
      MM[[i]] = y[i, ,]
    }
    M = lapply(MM, function(X) (X + t(X))/2)
    #print(M)
    LL = lapply(MM, chol)   
    L = lapply(LL, function(X) X - diag(diag(X)))
    D = lapply(LL, function(X) diag(X))
    
    ss = 0
    U = 0
    E = 0
    s <- 1 + t(t(Bin) - mx) %*% invVa %*% (xout - mx)
    s <- as.vector(s)
    for (i in 1:n) {
      ss = ss + s[i]
      U = U + s[i] * L[[i]]
      E = E + s[i] * log(D[[i]])}
    SS = U/ss + diag(exp(E/ss))
    y_hat = t(SS) %*% SS
    
    fit <- try(chol(y_hat),silent = TRUE)
    if ('try-error' %in% class(fit)){
      P = eigen(y_hat)$vectors
      Lambd_alpha = diag(pmax(0, eigen(y_hat)$values))
      y_hat = P %*% Lambd_alpha %*% t(P)
      y_hat = as.matrix(Matrix::forceSymmetric(y_hat))
    }
  }
  else if(type =='sphere'){
    
    n = nrow(Bin)
    p = ncol(Bin)
    xout <- Bout
    invVa = solve(var(Bin))
    mx = apply(Bin, 2, mean)
    
    s <- 1 + t(t(Bin) - mx) %*% invVa %*% (xout - mx)
    s <- as.vector(s)
    y_pred <- ExtrinsicMean(data=y,weights=s)
    y_hat <- Trans.Euclid(y_pred)
  }
  return(y_hat)
}

trans_x2B <- function(x,splits)
{
  n <- nrow(x)
  B <- matrix(1,nrow=n)
  for (i in c(1:length(splits))){
    vt <- splits[[i]]
    B_loc <- matrix(1,nrow=n)
    for (j in c(2:ncol(vt))){
      v <- vt[2,j]
      s <- vt[3,j]
      t <- vt[4,j]
      B_loc <- B_loc*pmax(0, s * (x[,v] - t))
    }
    B <- cbind(B,B_loc)
    
  }
  B <- B[,-1]
  return(B)
}

trans_x2B_zero <- function(x,splits)
{
  n <- nrow(x)
  B <- matrix(1,nrow=n)
  for (i in c(1:length(splits))){
    vt <- splits[[i]]    
    vt <- vt[-1,]
    B_loc <- matrix(1,nrow=n)
    for (j in c(1:dim(vt)[1])){
      v <- vt[j,'v']
      B_loc <- B_loc*x[,v]
    }
    B <- cbind(B,B_loc)
    
  }
  B <- B[,-1]
  return(B)
}

test_loss <- function(a,b,type)
{
  lof <- 0
  n <- dim(a)[1] 
  
  if(type == 'distribution'){
    nqSup <- 21     #分位数个数
    qSup <- seq(0,1,length.out = nqSup)
    for (j in c(1:n)){
      density1 <- a[j,]
      density2 <- b[j,]
      lof_loc <- pracma::trapz(qSup, (density1 - density2)^2)
      lof=lof+lof_loc
    }
  }
  else if(type == 'matrix'){
    for (j in c(1:n)){
      matrix1 <- a[j,,]
      matrix2 <- b[j,,]
      #print(eigen(matrix2)$values)
      lof_loc <- (DistCholesky(matrix1,matrix2))**2
      lof=lof+lof_loc
    }
  }
  else if(type == 'sphere'){
    for (j in c(1:n)){
      cross <- crossprod(a[j,],b[j,])[1,1]
      lof_loc <- acos(cross)^2
      lof=lof+lof_loc
    }
  }
  
  lof <- lof/n
  return(lof)
}

cppFunction('List forward_c (NumericVector y, NumericMatrix x, int degree, int mx,int sample_num,CharacterVector type){
    int n = x.nrow();
    int p = x.ncol();
    NumericMatrix B(n,1);
    for (int i =0; i<n ; i++){
      B(i,0)=1;
    }
    NumericVector a=NumericVector::create(0,0,NA_REAL,NA_REAL);
    
    DataFrame splits_first=DataFrame::create(Named("a") = a);
    List splits(mx+1);
    splits[0] = splits_first;
    int M = 1;
    Environment myEnv = Environment::global_env();
    Function loss = myEnv["loss"];
    Function scale_x = myEnv["scale_x"];
    
    NumericMatrix x_new = scale_x(Named("X_obs",x),Named("sample_num",sample_num));
    
    
    while(M <= mx){
      int best_m = 0;
      int best_v = 0;
      double best_t = 0;
      Rcout<<M;
      double lof_best = 1000000;
      for(int m=0; m<M ; m=m+1){
        DataFrame splits_loc = splits[m];
        if (splits_loc.ncol() > degree){
          continue ;
        }
        for (int v=0; v<p; v=v+1){
          NumericVector tt = x_new(_,v);
            for (int it =0; it<tt.size(); it++){
              double t = tt(it);
              NumericVector h1 = B(_,m) * pmax(0, 1 * (x(_,v) - t));
              NumericVector h2 = B(_,m) * pmax(0, (-1) * (x(_,v) - t));
              DataFrame Bnew = cbind(B, h1, h2);
              
              NumericMatrix B_test = internal::convert_using_rfunction(Bnew,"as.matrix");
        
              NumericVector lof = loss(Named("y", y), Named("B", B_test), Named("type",type));
              double lof1 = mean(lof);


              if (lof1 <lof_best){
                lof_best = lof1;
                best_m = m;
                best_v = v;
                best_t = t;

              }
            }
        }
      }
                          

    NumericVector h1 = B(_,best_m) * pmax(0, 1 * (x(_,best_v) - best_t));
    NumericVector h2 = B(_,best_m) * pmax(0, (-1) * (x(_,best_v) - best_t));
    
    B = cbind(B, h1, h2);
    
    NumericVector split_left=NumericVector::create(best_m+1,best_v+1,1,best_t);
    NumericVector split_right=NumericVector::create(best_m+1,best_v+1,-1,best_t);

    DataFrame splits_right = splits[best_m]; 
    splits_right.push_back(split_right);  
    DataFrame splits_left = splits[best_m]; 
    splits_left.push_back(split_left);

    splits[M] = splits_left  ; 
    splits[M+1] = splits_right ;
    M = M + 2  ;

    }
    
    List result = List::create(Named("B")=B, Named("splits")=splits);
    return result;
}')



matrix.exp <- function(A)
{
  eig <- eigen(A)
  EA=eig$vectors%*%diag(exp(eig$values))%*%t(eig$vectors)
  return((EA+t(EA))/2)
}

DistCholesky=function (A = NULL, B = NULL) 
{
  
  p <- dim(A)[1]
  M.a <- chol(A ,pivot = TRUE)
  D.a <- diag(M.a)
  L.a <- M.a - diag(D.a)
  M.b <- chol(B ,pivot = TRUE)
  D.b <- diag(M.b)
  L.b <- M.b - diag(D.b)
  dist <- sqrt(sum((L.a - L.b)^2) + sum((log(D.a) - log(D.b))^2))
  
  return(dist)
}


scale_x <- function(X_obs,sample_num)
{
  p <- ncol(X_obs)
  prob <- c(1:sample_num)
  prob_loc <- 1/sample_num*prob
  prob_loc <- prob_loc[-sample_num]
  
  X <- matrix(0,sample_num-1,p)
  
  for (i in c(1:p)){
    X_i <- X_obs[,i]
    X_qua <- quantile(X_i,probs=prob_loc)
    for (j in c(1:sample_num-1)){
      X[j,i] <- X_qua[j]
    }
  }
  #print(X)
  return(X)
}











