# All the functions needs to run the EMS
#  
#Last updated 3/7/2025

#"make an edit here

# https://pmc.ncbi.nlm.nih.gov/articles/PMC3404522/#sec16
# https://bmcresnotes.biomedcentral.com/articles/10.1186/s13104-022-06017-z
# https://www.r-bloggers.com/2011/08/qtl-analysis-in-r/

get.x.index = function(x.index.i){
  pos.x.i = which(x.index.i == 1)
  pos.x.i.num = length(pos.x.i)
  index.u = unique(t(combn(rep(c(0,1),pos.x.i.num),pos.x.i.num)))
  list(pos.x.i,index.u)
}

inheritance = function(chr,theta,marker.chr){
  unlist(lapply(rep(theta,chr), markov.c, marker.chr = marker.chr))
}

markov.c = function(theta,marker.chr){
  M = numeric(marker.chr)
  M[1] = sample(c(0,1), 1, replace = TRUE, prob= c(1/2,1/2))
  for (j in 2:marker.chr)
  {
    a = M[j-1]*theta+(1-M[j-1])*(1-theta)
    b = (1-M[j-1])*theta+M[j-1]*(1-theta)
    
    M[j] = sample(c(0,1), 1, replace = TRUE, prob = c(a,b))
  }
  M
} #recodes the markers as ons or zeros

get.prep.theta.hat = function(index.j, index.j.plus){
  marker.j = M[,index.j]
  marker.j.plus = M[,index.j.plus]
  x.miss = sort(unique(c(marker.j,marker.j.plus)))
  x.comple = subj.ind[-x.miss]
  X.sub = X[x.comple,c(index.j,index.j.plus)]
  c(length(x.comple),sum(abs(X.sub[,1] - X.sub[,2])))
}

get.s.x.i=function(x.design.i){
  count.s.x.i = 0
  for(k in 1:chr){
    for(j in 1:(marker.chr -1)){
      index.j = (k - 1)*marker.chr + j
      index.j.plus = index.j +1
      count.s.x.i = count.s.x.i + (abs(x.design.i[index.j] - x.design.i[index.j.plus]) == 1)
    }
  }
  count.s.x.i
}

get.theta.comp = function(X.index.i,x.design.i,gamma.est){
  pos.x.i = X.index.i[[1]]
  len.index.u = nrow(X.index.i[[2]])
  if(length(pos.x.i) == 0){
    theta.comp = get.s.x.i(x.design.i)
  }else{
    num.i = rep(NA,len.index.u)
    den.i = rep(NA,len.index.u)
    for(o in 1:len.index.u){
      x.design.temp = x.design.i
      x.design.temp[pos.x.i] = X.index.i[[2]][o,]
      s.x.i = get.s.x.i(x.design.temp)
      temp.o = gamma.est^s.x.i
      num.i[o] = s.x.i*temp.o
      den.i[o] = temp.o
    }
    theta.comp= sum(num.i)/sum(den.i)
  }
  theta.comp
}

get.beta.comp = function(index.c, index.m, X.index.i,x.design.i,gamma.est, y.i, var.est, beta.est){
  pos.x.i = X.index.i[[1]]
  len.index.u = nrow(X.index.i[[2]])
  len.index.m = length(index.m)
  len.c = length(index.c)
  if(length(pos.x.i) == 0){
    x.c.i = x.design.i[index.c]
    x.m.i = x.design.i[index.m]
    temp.w = x.m.i%*%t(x.m.i)
    first.term = temp.w
    second.term = temp.w
    third.term = x.m.i*y.i
    fourth.term = x.m.i%*%t(x.c.i)
  }else{
    w.o.i.num = array(NA, c(len.index.m,len.index.m,len.index.u))
    w.o.i.den = rep(NA, len.index.u)
    w.o.i.num.sec = matrix(NA, nrow = len.index.u, ncol = len.index.m)
    w.m.i.num = array(NA, c(len.index.m,len.index.m,len.index.u))
    w.m.i.den = rep(NA, len.index.u)
    w.m.i.num.sec = array(NA,c(len.index.m,len.c,len.index.u))
    for(o in 1:len.index.u){
      x.design.temp = x.design.i
      x.design.temp[pos.x.i] = X.index.i[[2]][o,]
      s.x.i = get.s.x.i(x.design.temp)
      temp.o = gamma.est^s.x.i
      x.m.i = x.design.temp[index.m]
      x.c.i = x.design.temp[index.c]
      temp.o.i = exp(-(2*var.est)^(-1)*(y.i - x.c.i%*%beta.est)^2)*temp.o
      w.o.i.num[, ,o] = as.numeric(temp.o.i)*x.m.i%*%t(x.m.i)
      w.o.i.den[o] = as.numeric(temp.o.i)
      w.o.i.num.sec[o,] = as.numeric(temp.o.i)*x.m.i
      w.m.i.num[, ,o] = temp.o*x.m.i%*%t(x.m.i) 
      w.m.i.den[o] = temp.o
      w.m.i.num.sec[, ,o] = temp.o*as.matrix(x.m.i)%*%t(as.matrix(x.c.i)) 
    }
    w.o.num.sum = apply(w.o.i.num,1:2,sum)
    w.o.den.sum = sum(w.o.i.den)
    w.m.num.sum = apply(w.m.i.num,1:2,sum)	
    w.m.den.sum = sum(w.m.i.den)
    w.m.num.sec.sum = apply(w.m.i.num.sec,1:2,sum)
    w.o.xy = apply(w.o.i.num.sec,2,sum)
    first.term = w.o.num.sum/w.o.den.sum
    second.term = w.m.num.sum/w.m.den.sum
    third.term = (w.o.xy/w.o.den.sum)*y.i 
    fourth.term = w.m.num.sec.sum/w.m.den.sum 
  }
  list(first.term,second.term,third.term,fourth.term)
}

get.beta.hat = function(index.c,index.m,X.index,X,gamma.est,y,var.est,beta.est){
  
  temp.prep = vector("list", n)
  for(i in 1:n){
    temp.prep[[i]] = get.beta.comp(index.c, index.m, X.index[[i]],X[i,],gamma.est, y[i], var.est, beta.est)
  }
  temp.1 = lapply(temp.prep,extract,1)
  temp.2 = lapply(temp.prep,extract,2)
  temp.3 = sapply(temp.prep,extract,3)
  temp.4 = lapply(temp.prep,extract,4)
  
  temp.1.s.t = lapply(1:len.o,extract.1,A=temp.1, B=I.o)
  temp.2.s.t = lapply(1:len.m,extract.1,A=temp.2, B=I.m)
  if(length(index.m) == 1){
    temp.3.s.t = temp.3[I.o]
  }else{
    temp.3.s.t = temp.3[,I.o]
  }
  # print(temp.3)
  temp.4.s.t = lapply(1:len.m,extract.1,A=temp.4, B=I.m)
  
  arrayA = abind(temp.1.s.t, along=3)
  temp.a = apply(arrayA,1:2,sum)
  #temp.a.1 = Reduce('+', temp.1.s.t) 
  arrayB = abind(temp.2.s.t, along=3)
  temp.b = apply(arrayB,1:2,sum)
  if(length(index.m) == 1){
    temp.c = sum(temp.3.s.t)
  }else{
    temp.c = apply(temp.3.s.t,1,sum)
  }
  arrayD = abind(temp.4.s.t, along=3)
  temp.d = apply(arrayD,1:2,sum)
  solve(temp.a + temp.b)%*%(temp.c +(temp.d%*%beta.c))
}

get.RSS = function(index.c, index.m, X.index,x.design,gamma.est, y, var.est, beta.est, beta.hat){
  first.term.sum = rep(NA, n)
  second.term.sum = rep(NA,n)
  for(i in 1:n){
    X.index.i = X.index[[i]]
    pos.x.i = X.index.i[[1]]
    len.index.u = nrow(X.index.i[[2]])
    len.index.m = length(index.m)
    x.design.i = x.design[i,]
    y.i = y[i]
    if(length(pos.x.i) == 0){
      x.c.i = x.design.i[index.c]
      x.m.i = x.design.i[index.m]
      temp.y.i = as.numeric(t(as.matrix(x.m.i))%*%beta.hat)
      first.term.sum[i] = (y.i - temp.y.i)^2
      second.term.sum[i] = (as.numeric(t(as.matrix(x.c.i))%*%beta.est) - temp.y.i)^2 
    }else{
      first.term.num = rep(NA, len.index.u)
      first.term.den = rep(NA, len.index.u)
      second.term.num = rep(NA, len.index.u)
      second.term.den = rep(NA, len.index.u)
      for(o in 1:len.index.u){
        x.design.temp = x.design.i
        x.design.temp[pos.x.i] = X.index.i[[2]][o,]
        s.x.i = get.s.x.i(x.design.temp)
        temp.o = gamma.est^s.x.i
        x.m.i = x.design.temp[index.m]
        x.c.i = x.design.temp[index.c]
        temp.o.i = as.numeric(exp(-(2*var.est)^(-1)*(y.i - x.c.i%*%beta.est)^2)*temp.o)
        y.i.hat.temp = t(as.numeric(x.m.i))%*%beta.hat
        first.term.num[o] = temp.o.i*sum((y.i - y.i.hat.temp)^2)
        first.term.den[o] = temp.o.i
        y.i.c.hat.temp = t(as.numeric(x.c.i))%*%beta.est 
        second.term.num[o] = temp.o*sum((y.i.c.hat.temp - y.i.hat.temp)^2)
        second.term.den[o] = temp.o
      }
      first.term.sum[i] = sum(first.term.num)/sum(first.term.den)
      second.term.sum[i] = sum(second.term.num)/sum(second.term.den)
    }
  }
  sum(first.term.sum[I.o]) + sum(second.term.sum[I.m])
}

get.model.dic = function(index.c, gamma.est,var.c,beta.c){
  iter = 1:p
  marker.fb = vector("list", 29)
  temp = NULL
  for(k in 1:15){
    len.iter = length(iter)
    RSS = numeric(len.iter)
    for(i in 1:len.iter){
      index.m = sort(c(temp,iter[i]))
      beta.hat = get.beta.hat(index.c,index.m,X.index,X,gamma.est,y,var.c,beta.c)
      RSS[i] = get.RSS(index.c,index.m,X.index,X,gamma.est,y,var.c,beta.c,beta.hat)
    }
    RSS.fb = min(RSS)
    m = iter[which(RSS == min(RSS))]
    iter = iter[-which(iter == m)]
    temp = c(temp, m)
    marker.fb[[k]] =temp
  }
  iter = sort(temp)
  
  for(k in 16:29){
    len.iter = length(iter)
    RSS = numeric(len.iter)
    for(i in 1:len.iter){
      index.m = iter[-which(iter ==iter[i])]
      beta.hat = get.beta.hat(index.c,index.m,X.index,X,gamma.est,y,var.c,beta.c)
      RSS[i] = get.RSS(index.c,index.m,X.index,X,gamma.est,y,var.c,beta.c,beta.hat)
    }
    m = iter[which(RSS == min(RSS))]
    iter = iter[-which(iter==m)]
    marker.fb[[k]]= iter
  }
  marker.fb
}

get.BIC = function(index.c, index.m, X.index,x.design,gamma.est, y, var.est, beta.est){
  beta.m = get.beta.hat(index.c,index.m,X.index,x.design,gamma.est,y,var.est,beta.est)
  first.term.sum = rep(NA, n)
  second.term.sum = rep(NA,n)
  for(i in 1:n){
    X.index.i = X.index[[i]]
    pos.x.i = X.index.i[[1]]
    len.index.u = nrow(X.index.i[[2]])
    len.index.m = length(index.m)
    x.design.i = x.design[i,]
    y.i = y[i]
    if(length(pos.x.i) == 0){
      x.c.i = x.design.i[index.c]
      x.m.i = x.design.i[index.m]
      temp.y.i = as.numeric(t(as.matrix(x.m.i))%*%beta.m)
      first.term.sum[i] = (y.i - temp.y.i)^2
      second.term.sum[i] = (as.numeric(t(as.matrix(x.c.i))%*%beta.est) - temp.y.i)^2 
    }else{
      first.term.num = rep(NA, len.index.u)
      first.term.den = rep(NA, len.index.u)
      second.term.num = rep(NA, len.index.u)
      second.term.den = rep(NA, len.index.u)
      for(o in 1:len.index.u){
        x.design.temp = x.design.i
        x.design.temp[pos.x.i] = X.index.i[[2]][o,]
        s.x.i = get.s.x.i(x.design.temp)
        temp.o = gamma.est^s.x.i
        x.m.i = x.design.temp[index.m]
        x.c.i = x.design.temp[index.c]
        temp.o.i = as.numeric(exp(-(2*var.est)^(-1)*(y.i - x.c.i%*%beta.est)^2)*temp.o)
        y.i.hat.temp = t(as.numeric(x.m.i))%*%beta.m
        first.term.num[o] = temp.o.i*sum((y.i - y.i.hat.temp)^2)
        first.term.den[o] = temp.o.i
        y.i.c.hat.temp = t(as.numeric(x.c.i))%*%beta.est
        second.term.num[o] = temp.o*sum((y.i.c.hat.temp - y.i.hat.temp)^2)
        second.term.den[o] = temp.o
      }
      first.term.sum[i] = sum(first.term.num)/sum(first.term.den)
      second.term.sum[i] = sum(second.term.num)/sum(second.term.den)
    }
  }
  n.m = length(I.m)
  n*log(n.m*var.est + sum(first.term.sum[I.o]) + sum(second.term.sum[I.m])) + length(index.m)*log(n)
}

get.var = function(index.c, index.m, X.index,x.design,gamma.est, y, var.est, beta.est){
  beta.m = get.beta.hat(index.c,index.m,X.index,x.design,gamma.est,y,var.est,beta.est)
  first.term.sum = rep(NA, n)
  second.term.sum = rep(NA,n)
  
  for(i in 1:n){
    X.index.i = X.index[[i]]
    pos.x.i = X.index.i[[1]]
    len.index.u = nrow(X.index.i[[2]])
    len.index.m = length(index.m)
    x.design.i = x.design[i,]
    y.i = y[i]
    if(length(pos.x.i) == 0){
      x.c.i = x.design.i[index.c]
      x.m.i = x.design.i[index.m]
      temp.y.i = as.numeric(t(as.matrix(x.m.i))%*%beta.m)
      first.term.sum[i] = (y.i - temp.y.i)^2
      second.term.sum[i] = (as.numeric(t(as.matrix(x.c.i))%*%beta.est) - temp.y.i)^2 
    }else{
      first.term.num = rep(NA, len.index.u)
      first.term.den = rep(NA, len.index.u)
      second.term.num = rep(NA, len.index.u)
      second.term.den = rep(NA, len.index.u)
      for(o in 1:len.index.u){
        x.design.temp = x.design.i
        x.design.temp[pos.x.i] = X.index.i[[2]][o,]
        s.x.i = get.s.x.i(x.design.temp)
        temp.o = gamma.est^s.x.i
        x.m.i = x.design.temp[index.m]
        x.c.i = x.design.temp[index.c]
        temp.o.i = as.numeric(exp(-(2*var.est)^(-1)*(y.i - x.c.i%*%beta.est)^2)*temp.o)
        y.i.hat.temp = t(as.numeric(x.m.i))%*%beta.m
        first.term.num[o] = temp.o.i*sum((y.i - y.i.hat.temp)^2)
        first.term.den[o] = temp.o.i
        y.i.c.hat.temp = t(as.numeric(x.c.i))%*%beta.est
        second.term.num[o] = temp.o*sum((y.i.c.hat.temp - y.i.hat.temp)^2)
        second.term.den[o] = temp.o
      }
      first.term.sum[i] = sum(first.term.num)/sum(first.term.den)
      second.term.sum[i] = sum(second.term.num)/sum(second.term.den)
    }
  }
  n.m = length(I.m)
  n^(-1)*(n.m*var.est + sum(first.term.sum[I.o]) + sum(second.term.sum[I.m]))
}

FB.selection = function(p, marker,y, twentynine){
  BIC = numeric(29)
  iter = 1:p
  marker.fb = vector("list", twentynine)
  temp = NULL
  for(k in 1:15){
    len.iter = length(iter)
    RSS = numeric(len.iter)
    for(i in 1:len.iter){
      a = sort(c(temp,iter[i])) # start of something
      X.sub = marker[,a]
      beta.hat = solve(t(X.sub)%*%X.sub)%*%t(X.sub)%*%y
      RSS[i] = t(y - X.sub%*%beta.hat)%*%(y - X.sub%*%beta.hat)
    }
    RSS.fb = min(RSS)
    m = iter[which(RSS == RSS.fb)]
    BIC[k] = log(RSS.fb) + length(beta.hat)*(log(n)/n)
    iter = iter[-which(iter == m)]
    temp = c(temp, m)
    marker.fb[[k]] =temp
  }
  iter = sort(temp)
  
  for(k in 16:29){
    len.iter = length(iter)
    RSS = numeric(len.iter)
    for(i in 1:len.iter){
      a = iter[-which(iter ==iter[i])] # this changes between loops
      X.sub = marker[,a]
      beta.hat = solve(t(X.sub)%*%X.sub)%*%t(X.sub)%*%y
      RSS[i] = t(y - X.sub%*%beta.hat)%*%(y - X.sub%*%beta.hat)
    }
    RSS.fb = min(RSS)
    m = iter[which(RSS == RSS.fb)]
    BIC[k] = log(RSS.fb) + length(beta.hat)*(log(n)/n)
    iter = iter[-which(iter==m)]
    marker.fb[[k]]= iter
  }
  
  BIC.fb = min(BIC)
  index = which(BIC == BIC.fb)
  if (length(index) >1){
    index = sample(index, 1, replace = FALSE, prob= c(1/length(index),length(index))) 
  }
  sort(marker.fb[[index]])
}


extract = function(X,a){
  X[[a]]
}

extract.1 = function(A,B,index){
  A[[B[index]]]
}


Generate_X_Index <- function(length,seed) {
  set.seed(seed)
  # List to store the results
  result <- list()
  
  # Loop to generate 100 elements
  for (i in 1:length) {
    # Generate random numbers for the vector (first element of each list element)
    num_elements <- sample(0:10, 1)  # Randomly choose a length of the vector 
    num_cols <- sample(2:20, 1)  # Randomly choose the number of columns for the matrix 
    vector_data <- sample(1:10, num_elements)  # Randomly choose numbers between 1 and 30 for the vector
    if (num_elements== 0 || num_cols == 0 || num_elements == 0) {
      num_cols = 0
      num_elements = 0
      vector_data = 0
    }
    vector_data <- sample(1:10, num_elements)  # Randomly choose numbers between 1 and 30 for the vector
    
    # Generate random matrix (second element of each list element)
    # num_rows <- num_elements  # Randomly choose the number of rows for the matrix (between 4 and 8)
    
    matrix_data <- matrix(sample(c(0,1), num_cols * num_elements, replace = TRUE), nrow = num_cols, ncol = num_elements)
    
    # Add the list element to the result
    result[[i]] <- list(vector_data, matrix_data)
  }
  
  return(result)
}

Generate_X_Index_Missing <- function(x, pctNA = 0.1, P = 0, seed) {
  set.seed(seed)
  n <- nrow(x)
  p <- ncol(x)
  
  if (P > n) stop("P cannot be greater than the number of rows.")
  
  # Select P rows to be fully complete (no missing values)
  complete_rows <- sample(1:n, P)
  
  # All indices in the matrix
  all_indices <- expand.grid(row = 1:n, col = 1:p)
  
  # Remove indices corresponding to rows that must be complete
  eligible_indices <- all_indices[!all_indices$row %in% complete_rows, ]
  
  # Total number of values that can be set to NA
  total_values <- n * p
  allowed_indices_count <- nrow(eligible_indices)
  num_NA <- min(floor(total_values * pctNA), allowed_indices_count)
  
  # Sample indices to insert NA
  NA_indices <- eligible_indices[sample(1:allowed_indices_count, num_NA), ]
  
  # Introduce NAs
  for (i in seq_len(nrow(NA_indices))) {
    x[NA_indices$row[i], NA_indices$col[i]] <- NA
  }
  
  return(x)
}

Generate_X_Row_Missing <- function(x, pctRows = 0.1, P = 0, seed) {
  if (!is.null(seed)) set.seed(seed)
  
  n <- nrow(x)
  
  # Number of rows to be made missing
  num_missing_rows <- floor(n * pctRows)
  
  # Choose P rows that must remain complete
  if (P > n) stop("P cannot be greater than the number of rows.")
  complete_rows <- sample(1:n, P)
  
  # Eligible rows (not in P)
  eligible_rows <- setdiff(1:n, complete_rows)
  
  # Sample rows to become missing
  missing_rows <- sample(eligible_rows, min(num_missing_rows, length(eligible_rows)))
  
  # Set entire rows to NA
  x[missing_rows, ] <- NA
  
  return(x)
}

generate_matrix <- function(from, to, rows, columns, seed) {
  # Generate values from 1 to 100
  set.seed(seed)
  values <- from:to
  
  # Repeat the values as needed to fill the 10x30 matrix
  values_repeated <- rep(values, length.out = rows * columns)
  random_values <- sample(from:to, rows * columns, replace = TRUE)
  
  # Create the matrix with 10 rows and 30 columns
  matrix_data <- matrix(random_values, nrow = rows, ncol = columns, byrow = TRUE)
  
  return(matrix_data)
}

generate_missing <- function(from, to, size, seed) {
  # Generate values from 1 to 100
  set.seed(seed)
  values <- from:to
  if (size == 0) {
    data_10<- values[1]
    data_90<- values
  } else {
  
    # # Repeat the values as needed to fill the 10x30 matrix
    # values_repeated <- rep(values, length.out = rows * columns)
    # random_values <- sample(from:to, rows * columns, replace = TRUE)
    
    # Create the matrix with 10 rows and 30 columns
    # matrix_data <- matrix(random_values, nrow = rows, ncol = columns, byrow = TRUE)
    
    indices <- sample(1:length(values), size = size)
    
    # Split data into 10% and 90% portions
    data_10 <- values[indices]  # 10% portion
    data_90 <- values[-indices]  # 90% portion
  }
  
  return(list(data_10, data_90))
}

generate_EMS_Simulation <- function(n, chr, marker.chr, missing, coef, seed, save) {
  
  source("EMS_Functions.R")
  library(MASS)
  library('plyr')
  library(abind)
  set.seed(seed)
  # original data for EMS example Sunil Rao
  n <<- n
  chr <<- chr
  marker.chr <<- marker.chr
  p <<- chr*marker.chr # total dim of regression model
  m.ind <<- 1:p  # 
  subj.ind <<- 1:n  # Define subject indicator (1 to n)
  # coef = 3#.1 to 3 works for the seed 3
  
  # Call the function to generate and print the matrix
  M <<- generate_matrix(1, n, 10, p, seed =seed)
  
  X.index <<- Generate_X_Index(n,seed = seed)
  save(X.index, file="EMS/X.index_notwork.rda")
  
  
  c.m <<- (marker.chr - 1)*chr
  
  
  
  theta <<- .2
  
  
  beta <- numeric(p)
  beta[1] <- coef[1]
  beta[2] <- coef[2]
  beta[3] <- coef[3]
  beta[7] <- coef[4]
  beta[8] <- coef[5]
  beta[15] <- coef[6]
  
  index <<- 1:p
  
  start.time <- Sys.time()
  result = vector("list", 1)
  r=1
  set.seed(seed)
  X <<- t(sapply(rep(chr,n),inheritance,theta = theta, marker.chr = marker.chr))
  
  if (missing == 0) {
    # tenpercent <- round(0.1 * n)
    # X.miss <- generate_missing(1, n, tenpercent, seed = seed)
    # I.m <<- X.miss[[1]]
    # I.o <<- X.miss[[2]]
    # print(I.o)
    # print(I.m)
    X.miss <- Generate_X_Row_Missing(X, P = p, seed = seed)
    vals <- 1:n
    I.o <<- vals[which(complete.cases(X.miss))]
    I.m <<- vals[-which(complete.cases(X.miss))]
    
  } else if (missing < 0) {
    X.miss <- generate_missing(1, n, 0, seed = seed)
    I.m <<- X.miss[[1]]
    I.o <<- X.miss[[2]]
  } else {

    X.miss <- Generate_X_Index_Missing(X, missing, P = p, seed = seed)
    vals <- 1:n
    I.o <<- vals[which(complete.cases(X.miss))]
    I.m <<- vals[-which(complete.cases(X.miss))]
    
  }
  
  len.m <<- length(I.m) # number of missing
  len.o <<- length(I.o) # number of observed
  I.index <<- numeric(n)
  I.index[I.m] = 1
  
  epsilon = 0.1*rnorm(n, mean = 0, sd =1)
  y <<- X%*%beta + epsilon
  index.c = index
  count.prep = matrix(NA, nrow = c.m, ncol =2)
  iter = 0
  for(k in 1:chr){
    for(j in 1:(marker.chr -1)){
      iter = iter +1
      index.j = (k - 1)*marker.chr + j
      index.j.plus = index.j +1
      count.prep[iter,] = get.prep.theta.hat(index.j, index.j.plus)
    }
  }
  N.count = apply(count.prep,2,sum)
  theta.pre = N.count[2]/N.count[1]
  gamma.est = theta.pre/(1-theta.pre)	
  
  eps = 1
  iteration = 1	
  while(eps > 10^(-4)& iteration <= 100){
    theta.i = rep(NA,n)
    for(i in 1:n){
      theta.i[i] = get.theta.comp(X.index[[i]],X[i,],gamma.est)
    }
    theta.hat = sum(theta.i)/(n*chr*(marker.chr-1))
    eps = abs(theta.hat - theta.pre)
    theta.pre = theta.hat
    gamma.est = theta.pre/(1-theta.pre)	
    iteration = iteration + 1
  }
  
  # print(X)
  y.o = y[I.o] 
  X.sub = X[I.o,]
  beta.c <<- solve(t(X.sub)%*%X.sub)%*%t(X.sub)%*%y.o
  y.o.hat <<- X.sub%*%beta.c
  var.c <<- as.numeric(t(y.o - y.o.hat)%*%(y.o - y.o.hat))/(len.o - p)
  check = 1
  
  
  
  while(check != 0){
    model.dic = get.model.dic(index.c, gamma.est,var.c,beta.c)
    model.dic = unique(lapply(model.dic,sort))
    len.model = length(model.dic)
    
    BIC.v = sapply(model.dic,get.BIC,index.c=index.c,X.index=X.index,x.design=X,gamma.est=gamma.est,y=y,var.est=var.c,beta.est=beta.c)
    #This if statement was added to stop specfic errors for N=50 with missing = .05
    if (length(BIC.v) == 0 || all(is.na(BIC.v))) {
      stop("BIC.v is empty or contains only NA values.")
      model.sel = model.sel
    } else {
      model.sel = model.dic[[which(BIC.v == min(BIC.v))]]
    }
    # model.sel = model.dic[[which(BIC.v == min(BIC.v))]]
    var.c <<- get.var(index.c,model.sel,X.index,X,gamma.est,y,var.c,beta.c)
    beta.c <<- get.beta.hat(index.c,model.sel,X.index,X,gamma.est,y,var.c,beta.c)
    
    if(length(model.sel) == length(index.c)){
      if(sum(model.sel %in% index.c) == length(index.c)){
        check = 0
      }else{
        index.c = model.sel
        check = check + 1
      }
    }else{
      index.c = model.sel
      check = check + 1
    }
  }
  
  result.1 = model.sel
  result.2 = FB.selection(p,X,y,29) #is the EMS method
  result[[r]] = list(result.1,result.2)
  
  result.1
  result.2
  saveseed = as.character(seed)
  # save(result, file="EMS/Sim_EMS_Seed_5_CoefChange")
  # save(result, file = paste0("EMS/Sims/", save, ".rda"))
  
  
  end.time <- Sys.time()
  time.taken <- end.time - start.time
  time.taken
  
  # fulldata <- list(X, y, I.m, result)
  fulldata <- list(X, y, I.m, result.1, result.2)
  
  # return(list(xlist = X, Ylist = y, I.m = I.m, results =result))
  print("### Run Complete ###")
  # cat("### Run Complete ###")
  AllResults <- list(xlist = X, Ylist = y, I.m = I.m, model.sel =result.1, EMS =result.2, FullMiss =X.miss)
  save(AllResults, file = paste0("Sims/", save, ".rda"))
  return(list(xlist = X, Ylist = y, I.m = I.m, model.sel =result.1, EMS =result.2, FullMiss =X.miss))
}


evaluate_methods <- function(top_selection, bottom_selection, correct_values, top_name = "Top", bottom_name = "Bottom") {
  calculate_metrics <- function(selections, correct_values) {
    selection_counts <- table(unlist(selections))
    correct_counts <- table(unlist(correct_values))
    
    all_items <- union(names(selection_counts), names(correct_counts))
    
    tp <- 0
    fp <- 0
    fn <- 0
    
    for (item in all_items) {
      sel_count <- ifelse(item %in% names(selection_counts), selection_counts[[item]], 0)
      cor_count <- ifelse(item %in% names(correct_counts), correct_counts[[item]], 0)
      
      tp <- tp + min(sel_count, cor_count)
      fp <- fp + max(sel_count - cor_count, 0)
      fn <- fn + max(cor_count - sel_count, 0)
    }
    
    precision <- if ((tp + fp) > 0) tp / (tp + fp) else 0
    recall <- if ((tp + fn) > 0) tp / (tp + fn) else 0
    f1 <- if ((precision + recall) > 0) 2 * (precision * recall) / (precision + recall) else 0
    
    return(list(
      true_positives = tp,
      false_positives = fp,
      false_negatives = fn,
      precision = precision,
      recall = recall,
      f1_score = f1
    ))
  }
  
  # Calculate metrics for each method
  top_metrics <- calculate_metrics(top_selection, correct_values)
  # print(paste("Top: ", top_selection, correct_values))
  bottom_metrics <- calculate_metrics(bottom_selection, correct_values)
  # print(paste("Bottom: ", bottom_selection, correct_values))

  # Build results data frame
  results <- data.frame(
    Method = c(top_name, bottom_name),
    True_Positives = c(top_metrics$true_positives, bottom_metrics$true_positives),
    False_Positives = c(top_metrics$false_positives, bottom_metrics$false_positives),
    False_Negatives = c(top_metrics$false_negatives, bottom_metrics$false_negatives),
    Precision = c(top_metrics$precision, bottom_metrics$precision),
    Recall = c(top_metrics$recall, bottom_metrics$recall),
    F1_Score = c(top_metrics$f1_score, bottom_metrics$f1_score)
  )

  # print(results)
  return(results)
}
