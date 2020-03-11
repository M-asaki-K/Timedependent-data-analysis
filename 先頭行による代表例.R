library(dplyr)

multi.regression.compounds.train <- matrix(c(0,0,0,1,1,1,2,2,2,2,3,4,5,6,7,8,9,10), ncol = 2)
View(multi.regression.compounds.train)

#--------------------------1列目に関し、値が変化した先頭行のみで代表------------------------------------
multi.regression.compounds.train <- as.data.frame(as.data.frame(multi.regression.compounds.train) %>% distinct(V1,.keep_all=TRUE))
View(multi.regression.compounds.train)
