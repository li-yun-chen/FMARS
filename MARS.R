forward <- function(y,x,degree,mx,alpha_list)
{
  n <- nrow(x)  
  p <- ncol(x) 
  B <- as.matrix(rep(1,n),nrow=n)
  splits = list(data.frame(
    m = 0,
    v = 0,
    s = NA,
    t = NA
  ))
  M <- 1
  
  while (M < mx) {
    lof_best = Inf
    for (m in 1:M) {   
      print(m)
      splits_loc <- splits[[m]]
      if (dim(splits_loc)[1]>degree)
      {
        next
      }
      
      for (v in c(1:p)) {      
        tt = split_points(x[, v], B[, m])    
        for (t in tt) {
          Bnew =data.frame(B, Btem1 = B[, m] * h(+1, x[, v], t), Btem2 = B[, m] * h(-1, x[, v], t))
          B_test <- as.matrix(Bnew[-1])
          
          fit <- try(solve(var(B_test)),silent = TRUE)
          if ('try-error' %in% class(fit)){
            next
          }
          
          lof = loss(y,B_test,alpha_list)  
          if (lof < lof_best) {
            lof_best = lof
            split_best = c(m = m,v = v,t = t)
          } 
        } 
      } 
    }
    m = split_best[['m']]
    v = split_best[['v']]
    t = split_best[['t']]
    B = cbind(B, B[, m] * h(+1, x[, v], t))
    B = cbind(B, B[, m] * h(-1, x[, v], t))     
    
    right = rbind(splits[[m]], c(m, v, -1, t))   
    left = rbind(splits[[m]], c(m, v, 1, t))
    splits = c(splits, list(left), list(right))   
    print(M)
    M = M + 2
  } 
  
  splits[[1]] <- NULL   
  B <- B[,2:ncol(B)]

  
  print(111111111)
  
  return(list(y = y, B = B, splits = splits,x = x))
  
}






