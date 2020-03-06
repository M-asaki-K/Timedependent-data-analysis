#----------------nothing to do with EDA（ガソリンのデータを引っ張ってくるだけです。時系列じゃねーじゃんとかそういうのやめて）--------------------
library(pls)

preprocessed.x <- as.matrix(gasoline)

#---------------時系列データ的CV-----------------------------------------
library(rsample)

multi.regression.compounds <- as.matrix(gasoline)

n.validation.int <- 5
n.validation <- n.validation.int + 1
ratio.train <- 3

split.num <- n.validation + ratio.train

### SPLIT DATA INTO K FOLDS ###
set.seed(2016)

df2fold <- rolling_origin(multi.regression.compounds, initial = round(nrow(multi.regression.compounds) / split.num * ratio.train), 
                          assess = round(nrow(multi.regression.compounds) / split.num), 
                          skip = round(nrow(multi.regression.compounds) / split.num), 
                          cumulative = TRUE)

#----------------------結果を見てみる----------------------------
multi.regression.compounds[df2fold$splits[[3]]$in_id, ]
multi.regression.compounds[df2fold$splits[[3]]$out_id, ]
