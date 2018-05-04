# setup work directiry
#setwd("~/GitHub/MultiVarFinalProject")

# load training data
load("STAT715_ink_data_training.rData")


head(ink.training.dat)

# The data is a 4 dim table which is not same as what we had
# we need 'down degree' it into a 2 dim table, first dim is observations(22*13), second is label and data(1+200*31)
my.table <- NULL
for(i in 1:13){
  for(j in 1:22){
    rawMatrix <- ink.training.dat[,,j,i]
    rawVector <- as.vector(rawMatrix)
    vectorWithLabel <- c(i,rawVector)
    my.table <- cbind(my.table, vectorWithLabel)
  }
}
colnames(my.table) <- NULL
dim(my.table) 
my.table<-t(my.table)
my.table[,1]

# Now we have table as what we using during class

# 1 mcLDA

# 1.1 split data and label
label <- my.table[,1]
data <- as.matrix(my.table[,-1])
dim(label)
dim(data)

# 1.2 compute mean and cov for each class
mcLDA.mean <- by(data, label, colMeans)
mcLDA.class.cov <- by(data,label, cov)
mcLDA.allClassMean <- colMeans(data)

# 1.3 compute between and within variance for each class
tmp <- lapply(1:13,
              function(x){13 * (mcLDA.mean[[x]]-mcLDA.allClassMean) %*% t(mcLDA.mean[[x]]-mcLDA.allClassMean)})
mcLDA.between.cov <- tmp[[1]]
dim(tmp)
for( i in 2:13 ){
  mcLDA.between.cov <- mcLDA.between.cov + tmp[[i]]
}
rm(tmp)

tmp <- lapply(1:13, 
       function(x){12 * mcLDA.class.cov[[x]]})
dim(tmp)
mcLDA.within.cov <- tmp[[1]]
for( i in 2:13 ){
  mcLDA.within.cov <- mcLDA.within.cov + tmp[[i]]
}
rm(tmp)

# 1.4 compute e-value and e-vector
mcLDA.eigen <- eigen(solve(mcLDA.within.cov) %*% mcLDA.between.cov)
#Error in solve.default(mcLDA.within.cov) : 
#system is computationally singular: reciprocal condition number = 1.59585e-24
# we cannot move on in this direction. 

# try apply kernel function
# 1.5 project data



# 1.x clean variables
rm(mcLDA.mean,mcLDA.class.cov,mcLDA.allClassMean,mcLDA.between.cov,mcLDA.within.cov,label,data)



# # 2 mcLDA with kernel
# 
# # 2.0 functions
# # 2.0.1 kernel: gaussian with sd = 2
# kernel_1 <- function(v1, v2){
#   exp(-1 * sum( (v1-v2)^2) ) / 2
# }
# 
# # 2.0.2 Mj function
# calcMj <- function(data, f){
#   return( Reduce('+',
#                  lapply(1:13,
#                         function(x){})))
# }
# 
# # 2.1 same as 1.1
# label <- my.table[,1]
# data <- as.matrix(my.table[,-1])
# dim(label)
# dim(data)
# 
# # 2.2 calc overall mean and
# mcLDA.mean <- by(data, label, colMeans)
# mcLDA.overallmean <- colMeans(data)
# 
# # 2.3 calc M matrix and N matrix
# M_j <- Reduce('+',
#               lapply(1:13, 
#                      function(x){ 1/13 }))
# 



# Based on paper


# 1.Support function
# libraries:
library(mvtnorm)


# generate matrix given dim n
get.P.Matrix <- function(n){
  result <- NULL
  for(i in 1:n){
    subDim = n-i
    subMatrix_0 = matrix(0,ncol = n-1-subDim, nrow = subDim)
    subMatrix_1 = matrix(1,ncol = 1, nrow=subDim)
    subMatrix_Diag = diag(subDim)
    subMatrix <- cbind(subMatrix_0,subMatrix_1,subMatrix_Diag)
    result <- rbind(result, subMatrix)
  }
  return(result)
}

# get score Matrix based on a kernel function
# obsevation index is first index
get.Score <- function(data, kernelFun){
  n = dim(data)[1]
  v <- sapply(1:n, function(x){
    sapply(1:n, function(y){
      kernelFun(data[x,], data[y,])
    })
  })
  return(v)
}

# calculate SSa
get.SSa <- function(Vec_Sn, Matrix_Vk, n){
    return(
      Reduce('+',
           sapply(2:n, function(x){
             ( t(Matrix_Vk[,x]) %*% Vec_Sn )^2
           })
           )
  )
}

# calculate SSe
get.SSe <- function(Vec_Sn, Matrix_Vk, n){
  return(
    Reduce('+',
           sapply((n+1):dim(Matrix_Vk)[1], function(x){
             ( t(Matrix_Vk[,x]) %*% Vec_Sn )^2
           })
    )
  )
}

# MSa
get.MSa <- function(SSa, n){
  return(SSa/(n-1))
}

# MSe
get.MSe <- function(SSe, n){
  N <- n*(n-1)/2
  return(SSe/(N-n))
}

# MSt 
get.MSt <- function(SSa, SSe, n){
  N <- n*(n-1)/2
  return((SSa+SSe)/(N-1))
}

# Calculate Capital Sigma
get.Cov <- function(var.a, var.e, n){
  P <- get.P.Matrix(n)
  return(P%*%t(P)*var.a + var.e * diag(n*(n-1)/2))
}

# Calc parms by given data and kernel function
get.Parms <- function(dat, kernel_fun){
  # observation count
  n <- dim(dat)[1]
  
  # score matrix
  score <- get.Score(dat, kernel_fun)
  
  # Sn
  Sn <- NULL
  for(j in 1:(n-1)){
    Sn <- c(Sn, score[j,][c((j+1):n)])
  }
  
  # P matrix
  P_matrix <- get.P.Matrix(n)
  tmp <- P_matrix %*% t(P_matrix)
  
  # eigen of P*Pt
  eigen <- eigen(tmp)
  
  # SSa 
  SSa <- get.SSa(Sn, eigen$vectors, n)
  
  # SSe
  SSe <- get.SSe(Sn, eigen$vectors, n)
  
  # MSa
  MSa <- get.MSa(SSa, n)
  
  # MSe
  Mse <- get.MSe(SSe, n)
  
  # mean hat
  mean <- mean(Sn)
  
  # variance
  var.a <- ( MSa - Mse ) / ( n - 2 )
  if(var.a < 0)
  {
    var.a = 0
    var.e = get.MSt(SSa, SSe, n)
  }
  else
  {
    var.e = Mse
  }
  
  result <- NULL
  result$a <- var.a
  result$e <- var.e
  result$mu <- mean
  result$sn <-Sn
  
  return(result)
}

# 2. Kernel function
ker_1 <- function(x,y){
  return(cor(x,y))
}


# 3. Main function
# 3.1 Training
#   training will give us a model with:
#   a.  13 means, 
#   b.  13 a_var, 
#   c.  13 e_var

trainModel <- function(traning.data){
  # clean data
  my.table <- NULL
  for(i in 1:13){
    for(j in 1:22){
      rawMatrix <- traning.data[,,j,i]
      rawVector <- as.vector(rawMatrix)
      vectorWithLabel <- c(i,rawVector)
      my.table <- cbind(my.table, vectorWithLabel)
    }
  }
  colnames(my.table) <- NULL
  #dim(my.table) 
  my.table<-t(my.table)
  my.table[,1]
  
  mu.array <- rep(1, 13)
  var.a.array <- rep(1, 13)
  var.e.array <- rep(1, 13)
  P.array <- rep(1, 13)
  
  # calc mean, var a, var e for each class
  # use this calc P(Sn)
  for(i in 1:13){
    # data in class i
    class.dat <- my.table[which(my.table[,1] == i),-1]
    
    class.result <- get.Parms(class.dat, ker_1)
    
    # assign result to matrix
    mu.array[i] = class.result$mu
    var.a.array[i] = class.result$a
    var.e.array[i] = class.result$e
    
    Sigma_Sn <- get.Cov(class.result$a, class.result$e, dim(class.dat)[1])
    
    P.array[i] <- dmvnorm(class.result$sn, mean=rep(class.result$mu, length(class.result$sn)), sigma=Sigma_Sn, log=TRUE)
    
  }  
  
  cat(var.e.array)
  
  result <- NULL
  result$mu <- mu.array
  result$a <- var.a.array
  result$e <- var.e.array
  result$p <- P.array
  result$table <- my.table
  
  return(result)
}

# 3.2 Classification
#   classification function will get 200*31 matrix
#   return which class it should belong to
classificationResult <- function(testData, model){
  # check testData
  if(
    dim(testData)[1] != 200 ||
    dim(testData)[2] != 31 ||
    dim(testData)[4] != 13
  ){
    cat('Test Data does not fit model. We need 200*31*p*13 sized test data.')
    cat('\r\n')
    return(NULL)
  }
  
  # clean test data, change into 6300 col ,  (13*p) row matrix
  testingSampleinEachGroup <- dim(testData)[3]
  testingTable <- NULL
  for(i in 1:13){
    for(j in 1:testingSampleinEachGroup){
      rawMatrix <- testData[,,j,i]
      rawVector <- as.vector(rawMatrix)
      vectorWithLabel <- c(i,rawVector)
      # Here we add 'label' to indicate these sample are belong to same group. But, we don't know which pen this group is belong to.
      # Be careful when using this label
      testingTable <- cbind(testingTable, vectorWithLabel)
    }
  }
  colnames(testingTable) <- NULL
  testingTable <- t(testingTable)
  
  # read model data
  mu.array <- model$mu
  var.a.array <- model$a
  var.e.array <- model$e
  p.array <- model$p
  trainingTable <- model$table
  
  
  # for each testing group i, 
  #   1.  Combine its data with training data with label j
  #   2.  Calculate Cap_Sigma_s using var_a, var_e from class j
  #   3.  Calculate the joint prob Pij( c( testingGroup_i, trainingGroup_j ) ) 
  #       And pull the prob of Pj( trainingGroup_j ) from model
  #   4.  Then conditional prob Pi|j = Pij / Pj
  #   5.  Find the j, such that Pi|j has the biggest value
  #   6.  j is the pen group i come from
  resultLabel <- rep(0, 13)
  highestProb <- rep(0, 13)
  for(i in 1:13){
    # testing data on group i, without group label
    groupData <- testingTable[which(testingTable[,1] == i), -1]
    for( j in 1:13){
      
      var.a <- var.a.array[j]
      var.e <- var.e.array[j]
      mu <- mu.array[j]
      Pj <- p.array[j]
      
      # training data for pen j, without pen label
      trainingData <- trainingTable[which(trainingTable[,1] == j), -1]
      
      # combine data
      dat <- rbind(groupData, trainingData)
      dat.n <- dim(dat)[1]
      
      # Calculate Cap_Sigma_s
      Sigma_s <- get.Cov(var.a, var.e, dim(dat)[1])
      
      # Calculate s
      # score matrix
      score <- get.Score(dat, kernel_fun)
      
      # S
      S <- NULL
      for(k in 1:(dat.n-1)){
        S <- c(S, score[k,][c((k+1):dat.n)])
      }
      
      # joint prob:
      Pij <- dmvnorm(S, mean=rep(mu, length(S)), sigma=Sigma_s, log = TRUE)
      
      # Pi|j
      P_cond <- Pij - Pj
      
      if( ( P_cond > highestProb[i]) ){
        #we find a higher Pi|j, replace existing
        resultLabel[i] = j
        highestProb[i] = P_cond
      }
    }
    
    
  }
  
  
  return(resultLabel)
}



















