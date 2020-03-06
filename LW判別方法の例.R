#-----------------------必要なパッケージのインストール-----------------------
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

#-------------------サンプルデータの設定--------------------------------------------
multi.regression.x <- as.matrix(gasoline)[, ]
multi.regression.x.train <- as.matrix(gasoline)[c(1:50), ]
multi.regression.x.test <- as.matrix(gasoline)[-c(1:50), ]

#-------------------------サンプル間ユークリッド距離の計算-------------------------------------------
my_distances1 <- parDist(as.matrix(multi.regression.x), method = "euclidean")
m <- as.big.matrix(as.matrix(my_distances1))
rm(my_distances1)
gc()

#-------------------------テストデータに対し、相対的に「近い」データの抽出-----------------------------
dd <- 0.3 #上位何割を抽出するか

nsamp <- round(nrow(multi.regression.x)*dd)
query <- matrix(0, nrow = nsamp, ncol = nrow(multi.regression.x.test))
#View(query)
low = nrow(multi.regression.x.train)
upper = nrow(multi.regression.x)

foreach(k = 1:nrow(multi.regression.x.test), .combine = cbind, .packages = c("foreach", "doParallel", "readr", "dplyr", "assertr", "rsample", "pls", "distances", "stats", "bigmemory")) %do% {
  ss <- low+1
  kkk <- low+k
  query[,c(k)] <- as.matrix(order(m[-c(ss:upper),c(kkk)]), nrow=1)[1:nsamp]
}

#------------------------各テストデータに対し「近い」行の表示-------------------------------------
View(query)