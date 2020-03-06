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
#----------------------------

x <- as.matrix(gasoline)
x <- as.data.frame(x[, c(1:5)]) #列数多すぎたので減らします

#-------------どれだけ遅らせるか、何パターン作るか---------------------
bb <- 1
cc <- 1

#---------時系列データ-----------
delay <- round(5*bb) #今回は最大5＋5行分だけ遅らせます（ここは任意）
range <- round(5*cc)

if (delay*range > 0) {
  
  stv <- ncol(x)
  cd <- matrix(0, ncol = ncol(x)*range, nrow = nrow(x))
  cd
  cd <- foreach(i = 1:stv, .combine = cbind, .packages = c("foreach", "doParallel", "readr", "dplyr", "assertr", "rsample", "pls", "distances", "stats")) %dopar% {
    foreach(j = 1:range, .combine = cbind, .packages = c("foreach", "doParallel", "readr", "dplyr", "assertr", "rsample", "pls", "distances", "stats")) %dopar% { 
      cd[, c(i + ncol(x)*(j - 1))] <- mutate(x, s = lag(x[,c(i)], delay*(j-1)))[,c(ncol(x) + 1)]
    }
  }
  
  colnames(cd) <- paste0("col_", 1:ncol(cd))
} else cd <- as.matrix(x)

#------------recollect data with delay---------------------
compounds.rev <- cd
is.completes.rev <- complete.cases(compounds.rev)
is.completes.rev
complete.compounds <- compounds.rev[is.completes.rev, ]

View(complete.compounds)
