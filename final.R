# load training data
load("~/GitHub/MultiVarFinalProject/STAT715_ink_data_training.rData")
# setup work directiry
setwd("~/GitHub/MultiVarFinalProject")


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
# 1.0 pre functions

# 1.1 split data and label
label <- my.table[,1]
data <- my.table[,-1]

mcLDA.mean <- by(data, label, mean)
mcLDA.cov <- by(data,label, cov)
mcLDA.grandMean <- colMeans(data)
mcLDA.withinCov <- 



# 1.x clean variables
rm(label,data)