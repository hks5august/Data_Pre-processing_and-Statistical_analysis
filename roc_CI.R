
library(pROC)
#####load file containing actual and predicted vlaues with header as actual and pred

ab <-read.table("ab",sep="\t",header=T)
#### calculate confidence intervals  or CI only ####
ci(ab$actual, ab$predicted)
##### calculate AUC and CI both #### 
rocobj <- roc(ab$actual, ab$predicted, ci=TRUE, of="auc")
rocobj$ci
