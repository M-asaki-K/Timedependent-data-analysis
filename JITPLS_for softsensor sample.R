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
Gauss_holdout <- function(bb, cc, dd){
bb <- 0.0038
cc <- 0.2316
dd <- 0.3387
   
#-----------select y from the dataset------------------
y <- complete.compounds[,c(1)]
y

#-----------select x from the dataset-----------------
x <- complete.compounds[,c(2:16)]

#-----------calculate standard distribution of x------
x.sds <- apply(x, 2, sd)

#-----------remove columns of 0 distribution from x----
sd.is.not.0 <- x.sds != 0
x <- x[, sd.is.not.0]

#-----------remove columns with too high cov-----------
#library(caret)

#df2 <- as.matrix(cor(x))
#hc = findCorrelation(df2, cutoff = .8, exact = FALSE) # putt any value as a "cutoff" 
#hc = sort(hc)
#reduced_Data = x[,-c(hc)]

#x <- reduced_Data
#rm(reduced_Data)

#---------時系列データ-----------
delay <- round(5*bb)
range <- round(5*cc)

if (delay*range > 0) {

stv <- ncol(x)
cd <- matrix(0, ncol = ncol(x)*range, nrow = nrow(x))

cd <- foreach(i = 1:stv, .combine = cbind, .packages = c("foreach", "doParallel", "readr", "dplyr", "assertr", "rsample", "pls", "distances", "stats")) %dopar% {
 foreach(j = 1:range, .combine = cbind, .packages = c("foreach", "doParallel", "readr", "dplyr", "assertr", "rsample", "pls", "distances", "stats")) %dopar% { 
   cd[, c(i + ncol(x)*(j - 1))] <- mutate(x, s = lag(x[,c(i)], delay*(j-1)))[,c(ncol(x) + 1)]
 }
}

colnames(cd) <- paste0("col_", 1:ncol(cd))
   } else cd <- as.matrix(x)

#------------recollect data with delay---------------------
compounds.rev <- cbind(y, cd)
is.completes.rev <- complete.cases(compounds.rev)
complete.compounds <- compounds.rev[is.completes.rev, ]
rm(compounds.rev)
#View(complete.compounds)

#-----------select y from the dataset------------------
y <- complete.compounds[,c(1)]

#-----------select x from the dataset-----------------
x <- complete.compounds[,c(2:ncol(complete.compounds))]
#View(x)

#-----------calculate standard distribution of x------
x.sds <- apply(x, 2, sd)

#-----------remove columns of 0 distribution from x----
sd.is.not.0 <- x.sds != 0
x <- x[, sd.is.not.0]
#View(x)

  
#-----------pick up columns if needed---------------------------
multi.regression.x <- x[ , ]

#-----------definition of multi.regression.compounds (used for MLR)--------
multi.regression.compounds <- cbind(y, x)

#--------------------divide into test and training data----------------------
train_size = 0.7

n = nrow(multi.regression.compounds)
#------------collect the data with n*train_size from the dataset------------
perm = sample(n, size = round(n * train_size))

#-------------------training data----------------------------------------
multi.regression.compounds.train <- multi.regression.compounds[c(1:round(n*train_size)), ]
preprocessed.y.train <- multi.regression.compounds.train[,c(1)]
multi.regression.x.train <- multi.regression.compounds.train[,-c(1)]
#-----------------------test data----------------------------------------
multi.regression.compounds.test <-multi.regression.compounds[-c(1:round(n*train_size)), ]
preprocessed.y.test <- multi.regression.compounds.test[,c(1)]
multi.regression.x.test <- multi.regression.compounds.test[,-c(1)]

#-----------transform into data frame--------------------------
multi.regression.compounds.train <- as.data.frame(multi.regression.compounds.train)
#multi.regression.compounds.train <- as.data.frame(as.data.frame(multi.regression.compounds.train) %>% distinct(y,.keep_all=TRUE))

#--------------------------variables elimination by importance threshold in SVM-------------------------------
multi.regression.compounds.train.s.t <- cbind(preprocessed.y.train, multi.regression.x.train[,])
#multi.regression.compounds.train.s.t <- as.data.frame(as.data.frame(multi.regression.compounds.train.s.t) %>% distinct(preprocessed.y.train,.keep_all=TRUE))
multi.regression.x.train.s.t <- multi.regression.compounds.train.s.t[,-c(1)]
multi.regression.compounds.test.s.t <- cbind(preprocessed.y.test,multi.regression.x.test)
multi.regression.x.test.s.t <- multi.regression.x.test[,]

multi.regression.x.r <- rbind(multi.regression.x.train.s.t, multi.regression.x.test.s.t)

# Euclidean distances
my_distances1 <- parDist(as.matrix(multi.regression.x.r), method = "euclidean")
m <- as.big.matrix(as.matrix(my_distances1))
rm(my_distances1)
gc()

nsamp <- round(nrow(multi.regression.compounds.train.s.t)*dd)
query <- matrix(0, nrow = nsamp, ncol = nrow(multi.regression.x.test.s.t))
#View(query)
low = nrow(multi.regression.x.train.s.t)
upper = nrow(multi.regression.x.r)

foreach(k = 1:nrow(multi.regression.x.test.s.t), .combine = cbind, .packages = c("foreach", "doParallel", "readr", "dplyr", "assertr", "rsample", "pls", "distances", "stats", "bigmemory")) %do% {
ss <- low+1
kkk <- low+k
query[,c(k)] <- as.matrix(order(m[-c(ss:upper),c(kkk)]), nrow=1)[1:nsamp]
}

git <- foreach(j = 1:nrow(multi.regression.compounds.test.s.t), .combine = rbind, .packages = c("foreach", "doParallel", "readr", "dplyr", "assertr", "rsample", "pls", "distances", "stats", "randomForest")) %dopar% {
   dat <- matrix(multi.regression.x.test.s.t[c(j),], ncol = ncol(multi.regression.x.test.s.t)) #蝗ｺ螳壹☆繧区擅莉ｶ縺後≠繧句?ｴ蜷医?ｯ縺薙％縺ｧ髯､蜴ｻ
#View(dat)
   
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

preprocessed.y.train.git <- preprocessed.y.train[query[, c(j)]]
centy <- mean(preprocessed.y.train.git)
sdy <- sd(preprocessed.y.train.git)
preprocessed.y.train.git <-(preprocessed.y.train.git - centy) / sdy
multi.regression.compounds.train.s.t.git <- as.data.frame(cbind(preprocessed.y.train.git, multi.regression.x.train.s.t.git))

model <- plsr(preprocessed.y.train.git~., data = multi.regression.compounds.train.s.t.git, validation = "CV")
ncompopt <- (order(as.matrix(model$validation$PRESS)))[1]
predic <- sdy*predict(model, newdata = dat)[,,ncompopt] + centy
#predic <- sdy*predict(model, newdata = dat) + centy
data.frame(preprocessed.y.test[j], predic)
}
plot(git[,c(1)], git[,c(2)])
#plot(git %>% distinct(preprocessed.y.test.j.,.keep_all=TRUE))
plot(git[, c(1)], xlim = c(0,300), ylim = c(1,5), ylab = "")
par(new = T)
plot(git[, c(2)], xlim = c(0,300), ylim = c(1,5), col = "blue")

#r2jit <- cor(git)[2]^2 * sign(cor(git)[2])
r2jit <- 1 - sum((git[,c(1)] - git[,c(2)])^2) / sum((mean(git[,c(1)]) - git[,c(2)])^2)
r2jit

Pred <- r2jit#最小化問題なら予測値に-1をかける

list(Score=Pred, Pred=Pred)}

opt_svm <- BayesianOptimization(Gauss_holdout,bounds=list(bb=c(0,0.5), cc=c(0,0.5), dd=c(0.1,0.5)), init_points=20, n_iter=50, acq='ei', kappa=2.576,eps=0.0, verbose=TRUE)

write.csv(opt_svm$History, "C:/Users/uni21/OneDrive/デスクトップ/Data science/JITPLSMIBayeshistoryss.csv")
