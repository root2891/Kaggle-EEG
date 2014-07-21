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


calclda <- function(variables,loadings)
{
  # find the number of samples in the data set
  as.data.frame(variables)
  numsamples <- nrow(variables)
  # make a vector to store the discriminant function
  ld <- numeric(numsamples)
  # find the number of variables
  numvariables <- length(variables)
  # calculate the value of the discriminant function for each sample
  for (i in 1:numsamples)
  {
    valuei <- 0
    for (j in 1:numvariables)
    {
      valueij <- variables[i,j]
      loadingj <- loadings[j]
      valuei <- valuei + (valueij * loadingj)
    }
    ld[i] <- valuei
  }
  # standardise the discriminant function so that its mean value is 0:
  ld <- as.data.frame(scale(ld, center=TRUE, scale=FALSE))
  ld <- ld[[1]]
  return(ld)
}

calculateR1<- function(neighbour1, index)
{
  #distMatrix <- neighbour1[neighbour1$original == index,]["distance"]
  distMatrix <- neighbour1$nn.dist[index,]
  return (mean(distMatrix))
}

calculateR2 <- function(embedded1, embedded2, neighbour2, index)
{
  #neighborIndex <- neighbour2[neighbour2$original == index,]["neighbor"]
  #distmatrix <- rdist(embedded1[index], embedded2[neighborIndex$neighbor])
  neighborIndex <- neighbour2$nn.index[index,]
  distmatrix <- rdist(embedded1[index], embedded2[neighborIndex])
  return (mean(distmatrix[1,]))
}

calculateS<- function(embedded1, embedded2, neighbour1, neighbour2)
{
  sValue <- 0
  for(index in 1:length(unique(neighbour1)))
  {
    
    #distMatrix <- neighbour1[neighbour1$original == index,]["distance"]
    #R1<- mean(distMatrix$distance)
    #neighborIndex <- neighbour2[neighbour2$original == index,]["neighbor"]
    #distmatrix <- rdist(embedded1[index], embedded2[neighborIndex$neighbor])
    #R2<- (mean(distmatrix[1,]))
    R1<- calculateR1(neighbour1, index)
    R2<- calculateR2(embedded1, embedded2,neighbour2, index)
    sValue<- sValue+R1/R2
  }
  return (sValue)
}

result <- data.frame(NULL)
corrValue <- data.frame(NULL)
#synchrony <- data.frame(NULL)
startValue <- totalNoofInterIctalFiles - noOfInterictalFiles+1
k<- 455
while(k <= totalNoofInterIctalFiles)
{
  print(k)
  if(noOfSeconds >1)
  {
    fileName<- paste("D:/Work Related/Kaggle EEG/Volumes/Seagate/seizure_detection/competition_data/clips/Dog_2/Dog_2_interictal_segment_",k:(k+noOfSeconds), sep = "")
  }
  else
  {
    fileName<- paste("D:/Work Related/Kaggle EEG/Volumes/Seagate/seizure_detection/competition_data/clips/Dog_2/Dog_2_interictal_segment_",k, sep = "")
  }
  fileName <- paste(fileName, ".mat", sep = "")
  patient7_ictal <-lapply(fileName,readMat)
  patient7_ictal.data <- patient7_ictal[[1]]$data
  if(noOfSeconds >1)
  {
    for(i in 2:noOfSeconds)
    {
      patient7_ictal.data <- cbind(patient7_ictal.data,patient7_ictal[[i]]$data)
    }
  }
  patient7_ictal.data <- as.ts(patient7_ictal.data)
  plot(patient7_ictal.data[seriesNo,], type = 'l')
  
  noOfDimension <-10
  output2.lag <- timeLag(patient7_ictal.data[seriesNo,],method = "acfdecor",plot.data = TRUE)
  output.lag <- output2.lag[1]
  output2.dimension <- FNN(patient7_ictal.data[3, ], dimension = noOfDimension)
  minfraction <- 100
  for(j in 1:noOfDimension)
  {
    if(output2.dimension[3,j] <= minfraction)
    {
      output.dimension <-j
      minfraction = output2.dimension[3,j]
    }
  }
  tempResult <- data.frame(NULL)
  tempCorr <- data.frame(NULL)
  #tempsynchrony<- data.frame(NULL)
  for(series1 in 1:16)
  {
    for(series2 in 1:16)
    {
      if(series1 != series2)
      {
        embeddedSeries1 <- embedd(patient7_ictal.data[series1,], output.dimension, output.lag)
        embeddedSeries2 <- embedd(patient7_ictal.data[series2,], output.dimension, output.lag)
        
        neighbour1<- get.knn(embeddedSeries1, k=5)
        #neighbour1 <- data.frame(neighbour = neighbour1$nn.index, distance= neighbour1$nn.dist)
        neighbour2<- get.knn(embeddedSeries2, k= 5)
        #neighbour2 <- data.frame(neighbour = neighbour2$nn.index, distance= neighbour2$nn.dist)
        
        sValue1 <- calculateS(embeddedSeries1,embeddedSeries2,neighbour1,neighbour2)
        sValue2<- calculateS(embeddedSeries2, embeddedSeries1,neighbour2,neighbour1)
        sValue <- (sValue1+sValue2)/2
        tempResult <- rbind(tempResult, sValue)
        corrVal  <- ccf(patient7_ictal.data[series1,],patient7_ictal.data[series2,], plot = FALSE)
        corrVal<- max(corrVal$acf)
        tempCorr <- rbind(tempCorr, corrVal)
        remove("embeddedSeries1","embeddedSeries2","neighbour1","neighbour2")
        #sync<- phase.sync(embeddedSeries1, embeddedSeries2)
        #tempsynchrony <- rbind(tempsynchrony, sync$pval)
      }
    }
  }
  tempResult <- t(tempResult)
  tempCorr <- t(tempCorr)
  # tempsynchrony<- t(tempsynchrony)
  names(tempResult) <- paste("V", 1:length(tempResult), sep = "")
  #names(tempsynchrony) <- <- paste("sync", 1:length(tempResult), sep = "")
  result <- rbind(result, tempResult)
  corrValue <- rbind(corrValue, tempCorr)
  print(paste("resultLength:",nrow(result)))
  #synchrony <- rbind(synchrony, tempsynchrony)
  k<- k+noOfSeconds
}

write.table(result, file = "resultData")
length(result)

k<-117
while(k <= noOfIctalFiles)
{
  if(noOfSeconds >1)
  {
    fileName<- paste("D:/Work Related/Kaggle EEG/Volumes/Seagate/seizure_detection/competition_data/clips/Dog_2/Dog_2_ictal_segment_",k:(k+noOfSeconds), sep = "")
  }
  else
  {
    fileName<- paste("D:/Work Related/Kaggle EEG/Volumes/Seagate/seizure_detection/competition_data/clips/Dog_2/Dog_2_ictal_segment_",k, sep = "")
    
  }
  fileName <- paste(fileName, ".mat", sep = "")
  patient7_ictal <-lapply(fileName,readMat)
  patient7_ictal.data <- patient7_ictal[[1]]$data
  if(noOfSeconds>1)
  {
    for(i in 2:noOfSeconds) 
    {
      patient7_ictal.data <- cbind(patient7_ictal.data,patient7_ictal[[i]]$data)
    }
  }
  
  patient7_ictal.data <- as.ts(patient7_ictal.data)
  plot(patient7_ictal.data[seriesNo,], type = 'l')
  
  seriesNo <- 1
  noOfDimension <-10
  output2.lag <- timeLag(patient7_ictal.data[seriesNo,],method = "acfdecor",plot.data = TRUE)
  output.lag <- output2.lag[1]
  output2.dimension <- FNN(patient7_ictal.data[3, ], dimension = noOfDimension)
  minfraction <- 100
  for(j in 1:noOfDimension)
  {
    if(output2.dimension[3,j] <= minfraction)
    {
      output.dimension <-j
      minfraction = output2.dimension[3,j]
    }
  }
  tempResult <- data.frame(NULL)
  tempCorr <- data.frame(NULL)
  #tempsynchrony<- data.frame(NULL)
  for(series1 in 1:16)
  {
    for(series2 in 1:16)
    {
      if(series1 != series2)
      {
        embeddedSeries1 <- embedd(patient7_ictal.data[series1,], output.dimension, output.lag)
        embeddedSeries2 <- embedd(patient7_ictal.data[series2,], output.dimension, output.lag)
        
        neighbour1<- get.knn(embeddedSeries1, k=5)
        #neighbour1 <- data.frame(neighbour = neighbour1$nn.index, distance= neighbour1$nn.dist)
        neighbour2<- get.knn(embeddedSeries2, k= 5)
        #neighbour2 <- data.frame(neighbour = neighbour2$nn.index, distance= neighbour2$nn.dist)
        
        sValue1 <- calculateS(embeddedSeries1,embeddedSeries2,neighbour1,neighbour2)
        sValue2<- calculateS(embeddedSeries2, embeddedSeries1,neighbour2,neighbour1)
        sValue <- (sValue1+sValue2)/2
        tempResult <- rbind(tempResult, sValue)
        
        corrVal  <- ccf(patient7_ictal.data[series1,],patient7_ictal.data[series2,], plot = FALSE)
        corrVal<- max(corrVal$acf)
        tempCorr <- rbind(tempCorr, corrVal)
        remove("embeddedSeries1","embeddedSeries2","neighbour1","neighbour2")
        #sync<- phase.sync(embeddedSeries1, embeddedSeries2)
        #tempsynchrony <- rbind(tempsynchrony, sync$pval)
      }
    }
  }
  tempResult <- t(tempResult)
  tempCorr <- t(tempCorr)
  #tempsynchrony<- t(tempsynchrony)
  names(tempResult) <- paste("V", 1:length(tempResult), sep = "")
  names(tempCorr) <- paste("V", 1:length(tempCorr), sep = "")
  #names(tempsynchrony) <- <- paste("sync", 1:length(tempResult), sep = "")
  result <- rbind(result, tempResult)
  corrValue <- rbind(corrValue, tempCorr)
  print(paste("result lenght:", nrow(result)))
  #synchrony <- rbind(synchrony, tempsynchrony)
  k<- k+noOfSeconds  
}

row.names(result) <- paste("record",1:nrow(result),sep = "")
write.table(result, file = "resultData")
write.table(result, file = "result.txt", row.names = FALSE, col.names=FALSE)
write.table(corrValue, file = "corrValue", row.names = FALSE, col.names=FALSE)

prinComponents <- princomp(result, cor= TRUE)
temp1 <- predict(prinComponents, result)
#create target data set
targetDataSet <- data.frame(NULL)
for(i in 1:(noOfInterictalFiles-15))
{
  targetDataSet <- rbind(targetDataSet, 0)
}

for(i in 1:15)
{
  targetDataSet <- rbind(targetDataSet, 1)
}

for(i in 1:noOfIctalFiles)
{
  targetDataSet <- rbind(targetDataSet, 2)
}

#classify preictal and interictal
interictalDataSet<- temp1[1:1148,]
interictalOutputDataSet <- targetDataSet[1:1148,]
trainTest<- createDataPartition(interictalOutputDataSet, p = .6, list = FALSE)
trainTarget <- interictalOutputDataSet[trainTest]
testTarget <- interictalOutputDataSet[-trainTest]
trainInput <- interictalDataSet[trainTest,]
testInput <- interictalDataSet[-trainTest,]

for(i in 1:240)
{
  plot(result[1000:1300,i], xlab = paste("number", i), type = 'l')
}
svmmodel <- svm(trainInput[,1], trainTarget, kernel = "radial")
plot(trainTarget, type = "l")
lines(svmmodel$fitted, type = "l", col = "red")
prediction <- predict(svmmodel, testInput)
plot(testTarget, type = "l")
lines(prediction, type = "l", col = "red")

#classify ictal and interictal
ictalDataSet <- result[1149:1318,]
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
trainInput <- result[trainTest,]
testInput <- result[-trainTest,]

svmmodel <- svm(trainInput, trainTarget, kernel = "radial")
plot(trainTarget, type = "l")
lines(svmmodel$fitted, type = "l", col = "red")
prediction <- predict(svmmodel, testInput)
plot(testTarget, type = "l")
lines(prediction, type = "l", col = "red")

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


