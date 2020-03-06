#----------------nothing to do with EDA（ガソリンのデータを引っ張ってくるだけです）--------------------
library(pls)

preprocessed.x <- as.matrix(gasoline)

#------------------removing the columns with too high correlation-----------------
library(caret)

df2 <- cor(preprocessed.x)
hc = findCorrelation(df2, cutoff=0.8) # putt any value as a "cutoff" 
hc = sort(hc)
reduced_Data = preprocessed.x[,-c(hc)]
preprocessed.x <- reduced_Data

#-----------pick up columns if needed---------------------------
multi.regression.x <- preprocessed.x[ , ]

View(multi.regression.x)
