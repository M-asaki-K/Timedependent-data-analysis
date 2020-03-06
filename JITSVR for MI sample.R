library(readr) # データ読み込み
library(dplyr) # データ操作一般
library(assertr) # データのチェック
library(rsample)

#Packages installation
library(genalg)
library(pls)
library(e1071)
library(kernlab)
library(iml)
library(devtools)
library(parallelDist)
library(bigmemory)
library(rBayesianOptimization)

pkgs <- c('foreach', 'doParallel')
lapply(pkgs, require, character.only = T)
#if you want to change the number of threads for the calculation, please change the value "detectCores()"
registerDoParallel(makeCluster(detectCores()))


#-----------pick up the file path--------------
path <- file.choose()
path

#-----------read csv file as compounds--------------
compounds <- read.csv(path)
#View(compounds)

#-----------remove some columns if needed--------------
trimed.compounds <- compounds[,]

#-----------select rows without empty cells---------
is.completes <- complete.cases(trimed.compounds)
is.completes

complete.compounds <- trimed.compounds[is.completes,]
#View(complete.compounds)

#row 201-203 should be neglected if you do not care about Mahalanobis distance
Gauss_holdout <- function(dd){
  
#-----------select y from the dataset------------------
y <- complete.compounds[,c(2)]
y

#-----------select x from the dataset-----------------
x <- x.s

#-----------calculate standard distribution of x------
x.sds <- apply(x, 2, sd)

#-----------remove columns of 0 distribution from x----
sd.is.not.0 <- x.sds != 0
x <- x[, sd.is.not.0]

#-----------remove columns with too high cov-----------
library(caret)

df2 <- as.matrix(cor(x))
hc = findCorrelation(df2, cutoff = .8, exact = FALSE) # putt any value as a "cutoff" 
hc = sort(hc)
reduced_Data = x[,-c(hc)]

x <- reduced_Data
rm(reduced_Data)

#-----------pick up columns if needed---------------------------
multi.regression.x <- x[ , ]

#-----------definition of multi.regression.compounds (used for MLR)--------
multi.regression.compounds <- as.data.frame(cbind(y, x))
multi.regression.compounds.d <- distinct(multi.regression.compounds,y, .keep_all = TRUE)

#--------------------divide into test and training data----------------------
set.seed(2016)
train_size = 0.9

n = nrow(multi.regression.compounds.d)
perm = sample(n, size = round(n * train_size))
#-------------------training data----------------------------------------
multi.regression.compounds.train <- multi.regression.compounds.d[perm, ]
y.train <- multi.regression.compounds.train[,c(1)]
multi.regression.x.train <- multi.regression.compounds.train[,-c(1)]
#-----------------------test data----------------------------------------
multi.regression.compounds.test <-multi.regression.compounds.d[-perm, ]
y.test <- multi.regression.compounds.test[,c(1)]
multi.regression.x.test <- multi.regression.compounds.test[,-c(1)]

#-----------transform into data frame--------------------------
multi.regression.compounds.train <- as.data.frame(multi.regression.compounds.train)
#--------------------------variables elimination by importance threshold in SVM-------------------------------
multi.regression.compounds.train.s.t <- cbind(y.train, multi.regression.x.train[,])
multi.regression.x.train.s.t <- multi.regression.x.train[,]
multi.regression.compounds.test.s.t <- cbind(y.test,multi.regression.x.test)
multi.regression.x.test.s.t <- multi.regression.x.test[,]

multi.regression.x.train.s.t.sc <- scale(multi.regression.x.train.s.t[,])
y.train.sc <- scale(y.train)
col_means_train <- attr(multi.regression.x.train.s.t.sc, "scaled:center") 
col_stddevs_train <- attr(multi.regression.x.train.s.t.sc, "scaled:scale")
multi.regression.x.test.s.t.sc <- scale(multi.regression.x.test.s.t, center = col_means_train, scale = col_stddevs_train)
y.test.sc <- (y.test - mean(y.train)) / sd(y.train)

multi.regression.x.sc <- rbind(multi.regression.x.train.s.t.sc, multi.regression.x.test.s.t.sc)

# Euclidean distances
my_distances1 <- parDist(multi.regression.x.sc, method = "euclidean")
m <- as.matrix(my_distances1)
rm(my_distances1)
gc()

nsamp <- round(dd*nrow(multi.regression.x.train.s.t))
query <- matrix(0, nrow = nsamp, ncol = nrow(multi.regression.x.test.s.t))
#View(query)
low = nrow(multi.regression.x.train.s.t)
upper = nrow(multi.regression.compounds.d)

foreach(k = 1:nrow(multi.regression.x.test.s.t), .combine = cbind, .packages = c("foreach", "doParallel", "readr", "dplyr", "assertr", "rsample", "pls", "distances", "stats", "bigmemory")) %do% {
ss <- low+1
kkk <- low+k
query[,c(k)] <- as.matrix(order(m[-c(ss:upper),c(kkk)]), nrow=1)[1:nsamp]
}
#View(query)

ggit <- foreach(j = 1:nrow(multi.regression.compounds.test.s.t), .combine = rbind, .packages = c("foreach", "doParallel", "readr", "dplyr", "assertr", "rsample", "pls", "distances", "stats", "e1071", "kernlab")) %dopar% {
dat <- as.matrix(multi.regression.x.test.s.t[c(j),], ncol = ncol(multi.regression.x.test.s.t)) #蝗ｺ螳壹☆繧区擅莉ｶ縺後≠繧句?ｴ蜷医?ｯ縺薙％縺ｧ髯､蜴ｻ

multi.regression.x.train.s.t.git <- multi.regression.x.train.s.t[query[, c(j)],]
x.sds <- apply(multi.regression.x.train.s.t.git, 2, sd)
sd.is.not.0 <- x.sds != 0
multi.regression.x.train.s.t.git <- multi.regression.x.train.s.t.git[, sd.is.not.0]
dat <- matrix(dat[sd.is.not.0], ncol = ncol(multi.regression.x.train.s.t.git))
multi.regression.x.train.s.t.git <- scale(multi.regression.x.train.s.t.git) 

# Use means and standard deviations from training set to normalize test set
col_means_train <- attr(multi.regression.x.train.s.t.git, "scaled:center") 
col_stddevs_train <- attr(multi.regression.x.train.s.t.git, "scaled:scale")
dat <- scale(dat, center = col_means_train, scale = col_stddevs_train)

preprocessed.y.train.git <- y.train[query[, c(j)]]
centy <- mean(preprocessed.y.train.git)
sdy <- sd(preprocessed.y.train.git)
preprocessed.y.train.git <-(preprocessed.y.train.git - centy) / sdy
multi.regression.compounds.train.s.t.git <- as.data.frame(cbind(preprocessed.y.train.git, multi.regression.x.train.s.t.git))

#generating SVM model with selected variables
gam <- matrix(data = 0, nrow = 31, ncol = 1)
for(k in -20:10){
   rbf <- rbfdot(sigma = 2^k)
   rbf
   
   asmat <- as.matrix(multi.regression.x.train.s.t.git)
   asmat
   
   kern <- kernelMatrix(rbf, asmat)
   sd(kern)
   gam[c(k + 21),] <- sd(kern)
}

hakata <- which.max(gam)

obj.se.t.t <- tune.svm(preprocessed.y.train.git~., data = multi.regression.compounds.train.s.t.git, gamma = 2^(hakata - 21), cost = 3, epsilon = 2^(-10:0), scale = FALSE)
obj.sc.t.t <- tune.svm(preprocessed.y.train.git~., data = multi.regression.compounds.train.s.t.git, gamma = 2^(hakata - 21), cost = 2^(-5:10), epsilon = obj.se.t.t$best.parameters[,c(3)], scale = FALSE)
obj.s.t.t <- tune.svm(preprocessed.y.train.git~., data = multi.regression.compounds.train.s.t.git, gamma = 2^(-20:10), cost = obj.sc.t.t$best.parameters[,c(2)], epsilon = obj.se.t.t$best.parameters[,c(3)], scale = FALSE)
model <- svm(preprocessed.y.train.git~., data = multi.regression.compounds.train.s.t.git, gammma = obj.s.t.t$best.parameters[,c(1)], cost = obj.s.t.t$best.parameters[,c(2)], epsilon = obj.s.t.t$best.parameters[,c(3)])

#model <- plsr(preprocessed.y.train.git~., data = multi.regression.compounds.train.s.t.git, validation = "CV")
#ncompopt <- (order(as.matrix(model$validation$PRESS)))[1]
#predic <- sdy*predict(model, newdata = dat)[,,ncompopt] + centy
predic <- sdy*predict(model, newdata = dat) + centy
data.frame(y.test[j], predic)
}
plot(ggit[,c(1)], ggit[,c(2)])

r2jit <- cor(ggit)[2]^2

#r2jit <- 1 - sum((ggit[,c(1)] - ggit[,c(2)])^2) / sum((mean(ggit[,c(1)]) - ggit[,c(2)])^2)
#r2jit
Pred <- r2jit#最小化問題なら予測値に-1をかける

list(Score=Pred, Pred=Pred)}

opt_svm <- BayesianOptimization(Gauss_holdout,bounds=list(dd=c(0.05,1)), init_points=5, n_iter=10, acq='ei', kappa=2.576,eps=0.0, verbose=TRUE)
write.csv(opt_svm$History, "C:/Users/uni21/OneDrive/デスクトップ/Data science/JITPLSMIBayeshistory.csv")
