#demo of density

source('FMARS.R')

beta <- c(0.1,0.2,0.3,0.4,rep(0,12),0.1,0.2,0.3,0.4)/2
n <- 500
n_xout <- 100
p <- 20
b_length <-  30
type <- 'distribution'

#generate X
data <- runif(n*p, 0.5, 1) 
X_obs <- matrix(data, n, p) #trainning data
data_out <- runif(n_xout*p, 0.5, 1)
xout <- matrix(data_out , n_xout, p) #testing data

#training sample
nqSup <- 21     #分位数个数
qSup <- seq(0,1,length.out = nqSup)
mu <- rep(0,n)
sigma <- rep(0,n)
Y_obs<- matrix(0,n,nqSup) 
for (i in 1:n) {
  a <- as.numeric(t(beta)%*%X_obs[i,])
  mu[i] <- rnorm(1,8*X_obs[i,1]^2*(2*a-1),0.2)
  sigma[i] <- 1
  Y_obs[i,] <- qnorm(qSup,mu[i],sigma[i])
}
Y_obs[,1] <- Y_obs[,2]
Y_obs[,nqSup] <- Y_obs[,nqSup-1] 

#test sample
mu_true <- rep(0,n_xout) 
sigma_true <- rep(0,n_xout)
M_true<- matrix(0,n_xout,nqSup)
for (i in 1:n_xout) {
  a <- as.numeric(t(beta)%*%xout[i,])
  mu_true[i] <- 8*xout[i,1]^2*(2*a-1)
  sigma_true[i] <- 1
  M_true[i,] <- qnorm(qSup, mu_true[i], sigma_true[i])
}
M_true[,1] <- M_true[,2]
M_true[,nqSup] <- M_true[,nqSup-1]

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









