#demo of SPD

source('FMARS.R')

beta_1 <- c(0.1,0.2,0.3,rep(0,17));beta_2 <- c(rep(0.2,8),0.2,0.2,0.2,0.2,rep(0.2,8));beta_3 <- c(rep(0,17),0.1,0.2,0.3)
n <- 500
n_xout <- 100
p <- 20
b_length <-  30
type <- 'matrix'

#generate X
data <- runif(n*p, 0.5, 1) 
X_obs <- matrix(data, n, p) #trainning data
data_out <- runif(n_xout*p, 0.5, 1)
xout <- matrix(data_out , n_xout, p) #testing data

#training sample
m <- 3
Y_obs=array(0,c(n,m,m))
for(j in 1:n){
  pho_1 <- 0.8*cos(as.numeric(t(beta_1)%*%X_obs[j,]*X_obs[j,1])*4*pi)
  pho_2 <- 0.6*cos(as.numeric(t(beta_2)%*%X_obs[j,]*X_obs[j,1])*4*pi)
  pho_3 <- 0.4*cos(as.numeric(t(beta_3)%*%X_obs[j,]*X_obs[j,1])*4*pi)
  D <- matrix(c(1,pho_1,pho_2,pho_1,1,pho_3,pho_2,pho_3,1),3,3)
  Z <- matrix(0,3,3)
  Z[1,1] <- rnorm(1,0,1)
  Z[2,2] <- rnorm(1,0,1)
  Z[3,3] <- rnorm(1,0,1)
  Z[1,2]=Z[2,1] <- rnorm(1,0,1/sqrt(2))
  Z[1,3]=Z[3,1] <- rnorm(1,0,1/sqrt(2))
  Z[2,3]=Z[3,2] <- rnorm(1,0,1/sqrt(2))
  logY <- 0.2*Z+D
  Y_obs[j,,] <-matrix.exp(logY)
}

#test sample
M_true=array(0,c(n_xout,m,m))
for(i in 1:n_xout){
  pho_1 <- 0.8*cos(as.numeric(t(beta_1)%*%X_obs[j,]*X_obs[j,1])*4*pi)
  pho_2 <- 0.6*cos(as.numeric(t(beta_2)%*%X_obs[j,]*X_obs[j,1])*4*pi)
  pho_3 <- 0.4*cos(as.numeric(t(beta_3)%*%X_obs[j,]*X_obs[j,1])*4*pi)
  D <- matrix(c(1,pho_1,pho_2,pho_1,1,pho_3,pho_2,pho_3,1),3,3)
  M_true[i,,] <- matrix.exp(D)
}

cl <- makeCluster(2) # The number of cores used in parallel computing

#global frechet regression
sample_num <- 0
mars.fit <- mars_train_c(Y_obs,X_obs,2,b_length,sample_num,type)
M_predict <- mars_predict_list(xout,mars.fit,Y_obs,X_obs,type)
model_loss <- test_loss(M_true,M_predict,type)
model_loss

#lcoal frechet regression
sample_num <- 5
mars.fit <- mars_train_c(Y_obs,X_obs,2,b_length,sample_num,type)
M_predict <- mars_predict_list(xout,mars.fit,Y_obs,X_obs,type)
model_loss <- test_loss(M_true,M_predict,type)
model_loss

























