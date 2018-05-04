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

# 2. Kernel function
ker_1 <- function(x,y){
  return(abs(sum(x,y)))
}


# 3. Main function
# 3.1 Training
#   training will give us a model with:
#   13 means and 13 variance 
#   original data
#   kernel function (maybe)
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
  dim(my.table) 
  my.table<-t(my.table)
  my.table[,1]
  
  mu.array <- seq(from = 1, to = 1, length.out = 13)
  var.a.array <- seq(from = 1, to = 1, length.out = 13)
  var.e.array <- seq(from = 1, to = 1, length.out = 13)
  
  # calc mean, var a, var e for each class
  for(i in 1:13){
    # data in class i
    class.dat <- my.table[which(my.table[,1] == i),-1]
    
    # observation count
    class.n <- dim(class.dat)[1]
    
    # score matrix
    class.score <- get.Score(class.dat, ker_1)
    
    # Sn
    class.Sn <- NULL
    for(j in 1:(class.n-1)){
      class.Sn <- c(class.Sn, class.score[j,][c((j+1):class.n)])
    }
    
    # P matrix
    P_matrix <- get.P.Matrix(class.n)
    tmp <- P_matrix %*% t(P_matrix)
    
    # eigen of P*Pt
    class.eigen <- eigen(tmp)
    
    # SSa 
    class.SSa <- get.SSa(class.Sn, class.eigen$vectors, class.n)
    
    # SSe
    class.SSe <- get.SSe(class.Sn, class.eigen$vectors, class.n)
    
    # MSa
    class.MSa <- get.MSa(class.SSa, class.n)
    
    # MSe
    class.Mse <- get.MSe(class.SSe, class.n)
    
    # mean hat
    class.mean <- mean(Sn)
    
    # variance
    class.var.a <- ( class.MSa - class.Mse ) / ( class.n - 2 )
    if(class.var.a < 0)
    {
      class.var.a = 0
      class.var.e = get.MSt(class.SSa, class.SSe, class.n)
    }
    else
    {
      class.var.e = class.Mse
    }
    
    # assign result to matrix
    mu.array[j] <- mean(class.Sn)
    var.a.array[j] <- class.var.a
    var.e.array[j] <- class.var.e
  }  
  
  
}

# 3.2 Classification
#   classification function will get 200*31 matrix
#   return which class it should belong to
classificationResult <- function(data){
  
}



















