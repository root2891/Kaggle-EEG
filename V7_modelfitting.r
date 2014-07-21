library("fNonlinear")
library("R.matlab")
library("tseriesChaos")
library("fractal")
library("fractaldim")
library("tseriesChaos")
library("fields")
library("FNN")
library("synchrony")
library("caret")
library("RSNNS")
library("e1071")

noOfNeighbors<-5
noOfSeconds <- 5
seriesNo<- 2
noOftrainingValues <-250
noOfIctalFiles <- 170
noOfInterictalFiles <-1148
frequency <- 400
totalNoofInterIctalFiles <- 1148
nooffiles <- noOfIctalFiles + noOfInterictalFiles

fileName<- paste("D:/Work Related/Kaggle EEG/Volumes/Seagate/seizure_detection/competition_data/clips/Dog_2/Dog_2_ictal_segment_",1:172, sep = "")
fileName <- paste(fileName, ".mat", sep = "")
patient7_ictal <-lapply(fileName,readMat)
patient7_ictal.data <- patient7_ictal[[1]]$latency

for(i in 2:172) 
{
  patient7_ictal.data <- cbind(patient7_ictal.data,patient7_ictal[[i]]$latency)
}

temp2 = data.frame(NULL)

inputData <- read.table("result.txt", header = FALSE,  sep = " ")
finalData <- data.frame(NULL)
for(k in 1:1318)
{
  temp = matrix(inputData[k,], ncol = 15, byrow = TRUE)
  temp3 = data.frame(NULL)
  for(i in 1:15)
  {
    for(j in 1:15)
    {
      if(j >=i)
      {
        temp3 <- rbind(temp3, temp[i,j])
      }
    }
  }
  temp3 <- t(temp3)
  finalData <- rbind(finalData, temp3)
}
write.table(finalData, file = "finalNonLinearINterdependence", row.names = FALSE, col.names=FALSE)
inputData <- finalData
#Classification between ictal and interictal
targetDataSet <- data.frame(NULL)
for(i in 1:noOfInterictalFiles)
{
  targetDataSet <- rbind(targetDataSet, 0)
}

for(i in 1:noOfIctalFiles)
{
  targetDataSet <- rbind(targetDataSet, 2)
}

#classification of early and later ictal
targetDataSet2 <- data.frame(NULL)
k<-1
while(k<=noOfIctalFiles)
{
  prevLatency <- -1
  i = 1
  currentLatency <- patient7_ictal.data[k]
  while((currentLatency > prevLatency)&(k<=noOfIctalFiles))
  {
    if(i<=15)
    {
      targetDataSet2 <- rbind(targetDataSet2, 0)
    }else
    {
      targetDataSet2 <- rbind(targetDataSet2, 1)
    }
    i = (i+1)
    k=(k+1)
    prevLatency = currentLatency
    currentLatency = patient7_ictal.data[k]
  }
}

plot(targetDataSet2[,1], type = "l")
lines(patient7_ictal.data[1,], type ="l", col = "red")

prinComponents <- princomp(result, cor= TRUE)
temp1 <- predict(prinComponents, result)

#classify preictal and interictal
ictalDataSet<- inputData[1149:1318,]
trainTest<- createDataPartition(targetDataSet2[,1], p = .75, list = FALSE)
trainInput <- ictalDataSet[trainTest,]
testInput <- ictalDataSet[-trainTest,]
trainTarget <- targetDataSet2[trainTest]
testTarget <- targetDataSet2[-trainTest]


svmmodel <- svm(trainInput, trainTarget, kernel = "radial",type = "nu-classification",nu = .2,probability = TRUE)
plot(trainTarget, type = "l")
lines(data.matrix(svmmodel$fitted), type = "l", col = "red")
prediction <- predict(svmmodel, testInput)
plot(testTarget)
lines(data.matrix(prediction), col = "red")

jmodel <- jordan(trainInput, trainTarget, size = c(5), inputsTest = testInput, targetsTest = testTarget)
plotIterativeError(jmodel)
plot(trainTarget, type="l")
lines(jmodel$fitted.values , col = "red")
plot(testTarget, type="l")
lines(jmodel$fittedTestValues , col = "red")

emodel <- elman(trainInput, trainTarget, size = c(50), inputsTest = testInput, targetsTest = testTarget)
plotIterativeError(emodel)
plot(trainTarget, type="l")
lines(emodel$fitted.values , col = "red")
plot(testTarget, type="l")
lines(emodel$fittedTestValues , col = "red")

rbfmodel <- rbf(trainInput, trainTarget, size = c(20), inputsTest = testInput, targetsTest = testTarget)
plotIterativeError(rbfmodel)
plot(trainTarget, type="l")
lines(rbfmodel$fitted.values , col = "red")
plot(testTarget, type="l")
lines(rbfmodel$fittedTestValues , col = "red")

mlpmodel<- mlp(trainInput, trainTarget, size = c(10), inputsTest = testInput, targetsTest = testTarget)
plotIterativeError(mlpmodel)
plot(trainTarget, type="l")
lines(mlpmodel$fitted.values , col = "red")
plot(testTarget, type="l")
lines(mlpmodel$fittedTestValues , col = "red")

art2model <- art2(trainInput, f2Units = 3)
plotIterativeError(art2model)
plot(trainTarget, type="l")
lines(art2model$fitted.values , col = "red")
plot(testTarget, type="l")
lines(mlpmodel$fittedTestValues , col = "red")


#classify ictal and interictal
interictalDataSet <- result[1:1148,]
targetDataSet <- data.frame(NULL)
for(i in 1:(noOfInterictalFiles))
{
  targetDataSet <- rbind(targetDataSet, 0)
}

for(i in 1:noOfIctalFiles)
{
  targetDataSet <- rbind(targetDataSet, 1)
}
trainTest<- createDataPartition(targetDataSet$X0, p = .75, list = FALSE)
trainTarget <- targetDataSet[trainTest]
testTarget <- targetDataSet[-trainTest]
trainInput <- inputData[trainTest,]
testInput <- inputData[-trainTest,]

svmmodel <- svm(trainInput, trainTarget, kernel = "radial", type = "nu-classification", nu = 0.1)
plot(trainTarget, type = "l")
lines(data.matrix(svmmodel$fitted), type = "l", col = "red")
prediction <- predict(svmmodel, testInput)
plot(testTarget, type = "l")
lines(data.matrix(prediction), type = "l", col = "red")

jmodel <- jordan(trainInput, trainTarget, size = c(5), inputsTest = testInput, targetsTest = testTarget)
plotIterativeError(jmodel)
plot(trainTarget, type="l")
lines(jmodel$fitted.values , col = "red")
plot(testTarget, type="l")
lines(jmodel$fittedTestValues , col = "red")

emodel <- elman(trainInput, trainTarget, size = c(50), inputsTest = testInput, targetsTest = testTarget)
plotIterativeError(emodel)
plot(trainTarget, type="l")
lines(emodel$fitted.values , col = "red")
plot(testTarget, type="l")
lines(emodel$fittedTestValues , col = "red")

rbfmodel <- rbf(trainInput, trainTarget, size = c(20), inputsTest = testInput, targetsTest = testTarget)
plotIterativeError(rbfmodel)
plot(trainTarget, type="l")
lines(rbfmodel$fitted.values , col = "red")
plot(testTarget, type="l")
lines(rbfmodel$fittedTestValues , col = "red")

mlpmodel<- mlp(trainInput, trainTarget, size = c(10), inputsTest = testInput, targetsTest = testTarget)
plotIterativeError(mlpmodel)
plot(trainTarget, type="l")
lines(mlpmodel$fitted.values , col = "red")
plot(testTarget, type="l")
lines(mlpmodel$fittedTestValues , col = "red")

art2model <- art2(trainInput, f2Units = 3)
plotIterativeError(art2model)
plot(trainTarget, type="l")
lines(art2model$fitted.values , col = "red")
plot(testTarget, type="l")
lines(mlpmodel$fittedTestValues , col = "red")

