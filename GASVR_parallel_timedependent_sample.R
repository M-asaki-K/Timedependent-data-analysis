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

#-----------select x from the dataset-----------------
x <- complete.compounds[,c(2:16)]
#View(x)

#-----------calculate standard distribution of x------
x.sds <- apply(x, 2, sd)
x.sds[1]

x.mean <- apply(x, 2, mean)
x.mean[2]
x.mean

#-----------remove columns of 0 distribution from x----
sd.is.not.0 <- x.sds != 0
x <- x[, sd.is.not.0]
#View(x)

#-----------select y from the dataset------------------
y <- complete.compounds[,c(1)]
y

#-----------standarization of y------------------------
preprocessed.y <- (y - mean(y)) / sd(y)
mean(preprocessed.y)
sd(preprocessed.y)

#-----------standarization of x------------------------
apply(x, 2, mean)
apply(x, 2, sd)
preprocessed.x <- apply(x, 2, function(x) {(x - mean(x)) / sd(x)})
#View(preprocessed.x)

#-----------x converted into data frame type for machine learning-----------
#class(preprocessed.x)
preprocessed.x <- data.frame(preprocessed.x)

#-----------compare the number of columns and rows--------------
ncol(preprocessed.x)
nrow(preprocessed.x)

library(caret)

df2 <- cor(preprocessed.x)
hc = findCorrelation(df2, cutoff=0.8) # putt any value as a "cutoff" 
hc = sort(hc)
reduced_Data = preprocessed.x[,-c(hc)]
preprocessed.x <- reduced_Data

#-----------pick up columns if needed---------------------------
multi.regression.x <- preprocessed.x[ , ]

#-----------definition of multi.regression.compounds (used for MLR)--------
multi.regression.compounds <- cbind(preprocessed.y, multi.regression.x)

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
#--------------------------variables elimination by importance threshold in SVM-------------------------------
multi.regression.compounds.train.s.t <- cbind(preprocessed.y.train, multi.regression.x.train[,])
multi.regression.x.train.s.t <- multi.regression.x.train[,]
multi.regression.compounds.test.s.t <- cbind(preprocessed.y.test,multi.regression.x.test)
multi.regression.x.test.s.t <- multi.regression.x.test[,]

Xtrain <- multi.regression.x.train.s.t
ytrain <- preprocessed.y.train
datasum <- as.data.frame(cbind(ytrain,Xtrain))

numberOfWavelenghts = ncol(Xtrain) - 1


evaluateNIR <- function(chromosome=c()) {
  n.validation.int <- 10
  n.validation <- n.validation.int + 1
  ratio.train <- 3
  
  split.num <- n.validation + ratio.train
  
  returnVal = 100
  minLV = 2
  if (sum(chromosome) < minLV) {
    returnVal
  } else {
    xtrain = Xtrain[,chromosome == 1];
    datasum = as.data.frame(cbind(ytrain, xtrain))
    
    df2 <- datasum
    ### SPLIT DATA INTO K FOLDS ###
    set.seed(2016)
    
    df2fold <- rolling_origin(multi.regression.compounds, initial = round(nrow(multi.regression.compounds) / split.num * ratio.train), 
                              assess = round(nrow(multi.regression.compounds) / split.num), 
                              skip = round(nrow(multi.regression.compounds) / split.num), 
                              cumulative = TRUE)
    ### PARAMETER LIST ###
    cost <- 3
    epsilon <- c(-10,-9,-8,-7,-6,-5,-4,-3,-2,-1, 0)

      gam <- foreach(k = -20:10, .combine = rbind,.packages = c("kernlab"))%dopar%{
      rbf <- rbfdot(sigma = 2^k)
      rbf
      
      asmat <- as.matrix(xtrain)
      asmat
      
      kern <- kernelMatrix(rbf, asmat)
      sd(kern)
      data.frame((k +21), sd(kern))
    }

    hakata <- which.max(gam$sd.kern.)

    gamma <- hakata - 21
    parms <- expand.grid(epsilon = epsilon, cost = cost, gamma = gamma)
    ### LOOP THROUGH PARAMETER VALUES ###
    result <- foreach(i = 1:nrow(parms), .combine = rbind, .packages = c("foreach", "doParallel", "readr", "dplyr", "assertr", "rsample")) %dopar% {
      c <- parms[i, ]$cost
      g <- parms[i, ]$gamma
      e <- parms[i, ]$epsilon
      ### K-FOLD VALIDATION ###
      out <- foreach(j = 1:(n.validation -1), .combine = rbind,.packages = c("foreach", "doParallel", "readr", "dplyr", "assertr", "rsample"), .inorder = FALSE) %dopar% {
        deve <- df2[df2fold$splits[[j]]$in_id, ]
        test <- na.omit(df2[df2fold$splits[[j]]$out_id, ])
        mdl <- e1071::svm(ytrain~., data = deve, cost = c, gamma = 2^g, epsilon = 2^e)
        pred <- predict(mdl, test)
        data.frame(test[, c(1)], pred)
        
      }
      ### CALCULATE SVM PERFORMANCE ###
      roc <- sum((out[, c(1)] - out[, c(2)])^2) / nrow(out)
      data.frame(parms[i, ], roc)
    }

    epsilon <- min(result[result[, c(4)] <= (min(result[,c(4)])), c(1)])
    cost <- c(-5,-4,-3,-2,-1,0,1,2,3,4,5,6,7,8,9,10)
    parms <- expand.grid(epsilon = epsilon, cost = cost, gamma = gamma)
    parms
    ### LOOP THROUGH PARAMETER VALUES ###
    result <- foreach(i = 1:nrow(parms), .combine = rbind, .packages = c("foreach", "doParallel", "readr", "dplyr", "assertr", "rsample")) %dopar% {
      c <- parms[i, ]$cost
      g <- parms[i, ]$gamma
      e <- parms[i, ]$epsilon
      ### K-FOLD VALIDATION ###
      out <- foreach(j = 1:(n.validation -1), .combine = rbind,.packages = c("foreach", "doParallel", "readr", "dplyr", "assertr", "rsample"), .inorder = FALSE) %dopar% {
        deve <- df2[df2fold$splits[[j]]$in_id, ]
        test <- na.omit(df2[df2fold$splits[[j]]$out_id, ])
        mdl <- e1071::svm(ytrain~., data = deve, cost = 2^c, gamma = 2^g, epsilon = 2^e)
        pred <- predict(mdl, test)
        data.frame(test[, c(1)], pred)
        
      }
      ### CALCULATE SVM PERFORMANCE ###
      roc <- sum((out[, c(1)] - out[, c(2)])^2) / nrow(out)
      data.frame(parms[i, ], roc)
    }
    
    cost <- min(result[(result[, c(4)] <= (min(result[,c(4)]))), c(2)])
    gamma <- c(-20,-19,-18,-17,-16,-15,-14,-13,-12,-11,-10,-9,-8,-7,-6,-5,-4,-3,-2,-1,0,1,2,3,4,5,6,7,8,9,10)
    parms <- expand.grid(epsilon = epsilon, cost = cost, gamma = gamma)
    ### LOOP THROUGH PARAMETER VALUES ###
    result <- foreach(i = 1:nrow(parms), .combine = rbind, .packages = c("foreach", "doParallel", "readr", "dplyr", "assertr", "rsample")) %dopar% {
      c <- parms[i, ]$cost
      g <- parms[i, ]$gamma
      e <- parms[i, ]$epsilon
      ### K-FOLD VALIDATION ###
      out <- foreach(j = 1:(n.validation -1), .combine = rbind,.packages = c("foreach", "doParallel", "readr", "dplyr", "assertr", "rsample"), .inorder = FALSE) %dopar% {
        deve <- df2[df2fold$splits[[j]]$in_id, ]
        test <- na.omit(df2[df2fold$splits[[j]]$out_id, ])
        mdl <- e1071::svm(ytrain~., data = deve, cost = 2^c, gamma = 2^g, epsilon = 2^e)
        pred <- predict(mdl, test)
        data.frame(test[, c(1)], pred)
        
      }
      ### CALCULATE SVM PERFORMANCE ###
      roc <- sum((out[, c(1)] - out[, c(2)])^2) / nrow(out)
      data.frame(parms[i, ], roc)
    }

    gamma <- min(result[(result[, c(4)] <= (min(result[,c(4)]))), c(3)])
    bestperformance <- min(result[, c(4)])
    returnVal = bestperformance
      }
}

monitor <- function(obj) {
  minEval = min(obj$evaluations);
  filter = obj$evaluations == minEval;
  bestObjectCount = sum(rep(1, obj$popSize)[filter]);
  # ok, deal with the situation that more than one object is best
  if (bestObjectCount > 1) {
    bestSolution = obj$population[filter,][1,];
  } else {
    bestSolution = obj$population[filter,];
  }
  outputBest = paste(obj$iter, " #selected=", sum(bestSolution),
                     " Best (Error=", minEval, "): ", sep="");
  for (var in 1:length(bestSolution)) {
    outputBest = paste(outputBest,
                       bestSolution[var], " ",
                       sep="");
  }
  outputBest = paste(outputBest, "\n", sep="");
  cat(outputBest);
}

#GA
itersinput = 20

nir.results = rbga.bin(size=numberOfWavelenghts, zeroToOneRatio=10,
                       evalFunc=evaluateNIR, monitorFunc=monitor,
                       popSize=100, iters= itersinput, mutationChance = 0.2, verbose=TRUE)

plot(nir.results$best, ylim = c(0,1), ylab = "", lty = 2, type = "b")
par(new=T)
plot(nir.results$mean, axes = FALSE, ylab = "MSECV", col = "blue", type = "b", pch = 2)
axis(4)

#最も成績の良かった変数群を抽出し設定
res <- as.matrix(nir.results$population)
res <- nir.results$population[nir.results$evaluations == min(nir.results$evaluations),][c(1),]
#res <- nir.results$population[nrow(min(nir.results$evaluations)),]
View(res)
res.is.not.0 <- res != 0
res.is.not.0

#GA selected features
multi.regression.compounds.train.s.t <- cbind(preprocessed.y.train, multi.regression.x.train[, res.is.not.0])
multi.regression.x.train.s.t <- multi.regression.x.train[, res.is.not.0]
multi.regression.compounds.test.s.t <- cbind(preprocessed.y.test,multi.regression.x.test[, res.is.not.0])
multi.regression.x.test.s.t <- multi.regression.x.test[, res.is.not.0]

View(multi.regression.compounds.train.s.t)

write.csv(multi.regression.x.train.s.t, file = "C:/Users/u00140035/Desktop/timedependentga0129.csv")

Xtrain <- multi.regression.x.train.s.t
ytrain <- preprocessed.y.train
datasum <- as.data.frame(cbind(ytrain,Xtrain))

    xtrain = Xtrain
    datasum = as.data.frame(cbind(ytrain, xtrain))
    
    df2 <- datasum
    n.validation.int <- 5
    n.validation <- n.validation.int + 1
    ratio.train <- 2
    
    split.num <- n.validation + ratio.train
    
    ### SPLIT DATA INTO K FOLDS ###
    set.seed(2016)
    
    df2fold <- rolling_origin(multi.regression.compounds, initial = round(nrow(multi.regression.compounds) / split.num * ratio.train), 
                              assess = round(nrow(multi.regression.compounds) / split.num), 
                              skip = round(nrow(multi.regression.compounds) / split.num), 
                              cumulative = FALSE)
    ### PARAMETER LIST ###
    cost <- 3
    epsilon <- c(-10,-9,-8,-7,-6,-5,-4,-3,-2,-1, 0)
    
    gam <- foreach(k = -20:10, .combine = rbind,.packages = c("kernlab"))%dopar%{
      rbf <- rbfdot(sigma = 2^k)
      rbf
      
      asmat <- as.matrix(xtrain)
      asmat
      
      kern <- kernelMatrix(rbf, asmat)
      sd(kern)
      data.frame((k +21), sd(kern))
    }
    
    hakata <- which.max(gam$sd.kern.)
    
    gamma <- hakata - 21
    parms <- expand.grid(epsilon = epsilon, cost = cost, gamma = gamma)
    ### LOOP THROUGH PARAMETER VALUES ###
    result <- foreach(i = 1:nrow(parms), .combine = rbind, .packages = c("foreach", "doParallel", "readr", "dplyr", "assertr", "rsample")) %dopar% {
      c <- parms[i, ]$cost
      g <- parms[i, ]$gamma
      e <- parms[i, ]$epsilon
      ### K-FOLD VALIDATION ###
      out <- foreach(j = 1:(n.validation -1), .combine = rbind,.packages = c("foreach", "doParallel", "readr", "dplyr", "assertr", "rsample"), .inorder = FALSE) %dopar% {
        deve <- df2[df2fold$splits[[j]]$in_id, ]
        test <- na.omit(df2[df2fold$splits[[j]]$out_id, ])
        mdl <- e1071::svm(ytrain~., data = deve, cost = c, gamma = 2^g, epsilon = 2^e)
        pred <- predict(mdl, test)
        data.frame(test[, c(1)], pred)
        
      }
      ### CALCULATE SVM PERFORMANCE ###
      roc <- sum((out[, c(1)] - out[, c(2)])^2) / nrow(out)
      data.frame(parms[i, ], roc)
    }
    
    epsilon <- min(result[result[, c(4)] <= (min(result[,c(4)])), c(1)])
    cost <- c(-5,-4,-3,-2,-1,0,1,2,3,4,5,6,7,8,9,10)
    parms <- expand.grid(epsilon = epsilon, cost = cost, gamma = gamma)
    ### LOOP THROUGH PARAMETER VALUES ###
    result <- foreach(i = 1:nrow(parms), .combine = rbind, .packages = c("foreach", "doParallel", "readr", "dplyr", "assertr", "rsample")) %dopar% {
      c <- parms[i, ]$cost
      g <- parms[i, ]$gamma
      e <- parms[i, ]$epsilon
      ### K-FOLD VALIDATION ###
      out <- foreach(j = 1:(n.validation -1), .combine = rbind,.packages = c("foreach", "doParallel", "readr", "dplyr", "assertr", "rsample"), .inorder = FALSE) %dopar% {
        deve <- df2[df2fold$splits[[j]]$in_id, ]
        test <- na.omit(df2[df2fold$splits[[j]]$out_id, ])
        mdl <- e1071::svm(ytrain~., data = deve, cost = 2^c, gamma = 2^g, epsilon = 2^e)
        pred <- predict(mdl, test)
        data.frame(test[, c(1)], pred)
        
      }
      ### CALCULATE SVM PERFORMANCE ###
      roc <- sum((out[, c(1)] - out[, c(2)])^2) / nrow(out)
      data.frame(parms[i, ], roc)
    }
    
    cost <- min(result[(result[, c(4)] <= (min(result[,c(4)]))), c(2)])
    gamma <- c(-20,-19,-18,-17,-16,-15,-14,-13,-12,-11,-10,-9,-8,-7,-6,-5,-4,-3,-2,-1,0,1,2,3,4,5,6,7,8,9,10)
    parms <- expand.grid(epsilon = epsilon, cost = cost, gamma = gamma)
    ### LOOP THROUGH PARAMETER VALUES ###
    result <- foreach(i = 1:nrow(parms), .combine = rbind, .packages = c("foreach", "doParallel", "readr", "dplyr", "assertr", "rsample")) %dopar% {
      c <- parms[i, ]$cost
      g <- parms[i, ]$gamma
      e <- parms[i, ]$epsilon
      ### K-FOLD VALIDATION ###
      out <- foreach(j = 1:(n.validation -1), .combine = rbind,.packages = c("foreach", "doParallel", "readr", "dplyr", "assertr", "rsample"), .inorder = FALSE) %dopar% {
        deve <- df2[df2fold$splits[[j]]$in_id, ]
        test <- na.omit(df2[df2fold$splits[[j]]$out_id, ])
        mdl <- e1071::svm(ytrain~., data = deve, cost = 2^c, gamma = 2^g, epsilon = 2^e)
        pred <- predict(mdl, test)
        data.frame(test[, c(1)], pred)
        
      }
      ### CALCULATE SVM PERFORMANCE ###
      roc <- sum((out[, c(1)] - out[, c(2)])^2) / nrow(out)
      data.frame(parms[i, ], roc)
    }
    
    gamma <- min(result[(result[, c(4)] <= (min(result[,c(4)]))), c(3)])
    bestperformance <- min(result[, c(4)])
    bestperformance
    
    gamma
    cost
    epsilon

    compounds.svr.s.t <- svm(multi.regression.x.train.s.t, preprocessed.y.train, gamma = 2^gamma, cost = 2^cost, epsilon = 2^epsilon)
    
    ttest <- predict(compounds.svr.s.t, newdata = multi.regression.x.test.s.t)
    ttrain <- predict(compounds.svr.s.t)
    
    plot(0, 0, type = "n", xlim = c(-3, 3), ylim = c(-3,3),xlab = "Observed Value", ylab = "Predicted Value")
    
    points(preprocessed.y.test, ttest, col = "orange", pch = 2)
    points(preprocessed.y.train, ttrain, col = "darkgray", pch = 3)
    abline(a=0, b=1)
    
    r2ttest <- cor(preprocessed.y.test, ttest)^2
    r2ttest
    
    r2ttrain <- cor(preprocessed.y.train, ttrain)^2
    r2ttrain
    