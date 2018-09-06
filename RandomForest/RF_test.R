#install.packages("pROC")
library(pROC)
library(randomForest)
library(ROCR)
require(verification)
args <- commandArgs(TRUE)
traindata_file=args[1]
load(traindata_file)
test_data_file=args[2]
testd=read.table(test_data_file,header=TRUE)
Res=testd[,2]
AA=testd[,1]
testd=testd[,-1]
testd=testd[,-1]
Score=predict(AllFit,newdata=testd,type="prob")[,2]
output=data.frame(Res,AA,Score)
outputfile=args[3]
write.csv(output,outputfile)
