
remove(list=ls())

data(iris)
dat <- iris
head(dat)

class.labs <- dat[,5]
dat.no.labs <- as.matrix(dat[,-5])
n.classes <- length(unique(class.labs))
class.labs <- as.numeric(class.labs)

### get mean vectors
mean.vects <- by(dat.no.labs, class.labs, colMeans)
grand.mean <- colMeans(dat.no.labs)

### get class covariance matrices
cov.mats <- by (dat.no.labs, class.labs, cov)

### get within covariance
within.cov <- (1/length(class.labs))*
  Reduce("+", 
         lapply(1:n.classes, 
                function(x, class.labs, cov.mats)
                      {(sum(class.labs==x)-1)*cov.mats[[x]]}, 
                class.labs, 
                cov.mats
                )
         )

### get between covariance
between.cov <- Reduce("+",
                      lapply(1:n.classes, function(x, class.labs, mean.vects, grand.mean){
                                                          sum(class.labs==x)/length(class.labs)*(mean.vects[[x]]-grand.mean)%*%t(mean.vects[[x]]-grand.mean)
                                                  },
                                class.labs, mean.vects, grand.mean
                             )
                      )

### do eigenvalue decomposition
eigen.vects <- eigen(solve(within.cov)%*%between.cov)$vectors
eigen.vals <- eigen(solve(within.cov)%*%between.cov)$values

### transform data
proj.dat <- as.matrix(dat.no.labs)%*%matrix(as.numeric(eigen.vects[,1:(n.classes-1)]), dim(eigen.vects)[1], n.classes-1)
plot(proj.dat[,1], proj.dat[,2], col=class.labs, pch=16)

### predict classes
proj.means <- by(proj.dat, class.labs, colMeans)
proj.means <- rbind(proj.means[[1]], proj.means[[2]], proj.means[[3]])

pred.labs <- apply(proj.dat, 1, function(x, proj.means){which.min(diag((x-proj.means)%*%t(x-proj.means)))}, proj.means)

table(pred.labs, class.labs)
sum(pred.labs==class.labs)/length(class.labs)
