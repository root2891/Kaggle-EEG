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
library("wavelets")
library("scatterplot3d")
library("plot3D")
library("igraph")
library("rgl")

noOfSeconds <- 1
seriesNo<- 2
noOfIctalFiles <- 170
noOfInterictalFiles <-1148
frequency <- 400
totalNoofInterIctalFiles <- 1148
nooffiles <- noOfIctalFiles + noOfInterictalFiles
noOfTestFiles <- 100

calc_splv_spec <- function(wav_transform1, wav_transform2)
{
  sum <- complex(imaginary  = 0, real = 0)
  for(i in 1:nrow(wav_transform1))
  {
    tempVal <- exp(complex(imaginary=(wav_transform1[i,] - wav_transform2[i,])))
    realPart <- Re(sum)+Re(tempVal)
    imaginaryPart <- Im(sum)+Im(tempVal)
    sum <- complex(real = realPart, imaginary = imaginaryPart)
  }
  val <- abs(sum)
  return (val/nrow(wav_transform1))
}

calc_splv <- function(wav_transform1, wav_transform2)
{
  splv <- data.frame(NULL)
  for(i in 1:length(wav_transform1@W))
  {
    tempVal <- calc_splv_spec(data.frame(wav_transform1@W[i]), data.frame(wav_transform2@W[i]))    
    splv <- rbind(splv, tempVal)
  }
  return (splv)
}

result <- list(NULL)
startValue <- totalNoofInterIctalFiles - noOfInterictalFiles+1
k<- startValue
while(k <= totalNoofInterIctalFiles)
{
  print(k)
  if(noOfSeconds >1)
  {
    fileName<- paste("D:/Work Related/Kaggle EEG/Volumes/Seagate/seizure_detection/competition_data/clips/Dog_2/Dog_2_interictal_segment_",k:(k+noOfSeconds), sep = "")
  }else
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
  #plot(patient7_ictal.data[seriesNo,], type = 'l')
  
  tempResult <- data.frame(NULL)
  for(series1 in 1:16)
  {
    for(series2 in 1:16)
    {
      if(series1 > series2)
      {
        wavelet_transform1 <- dwt(patient7_ictal.data[series1,], n.levels = 8)
        wavelet_transform2 <- dwt(patient7_ictal.data[series2,], n.levels = 8)
        tempResult <-rbind(tempResult, t(calc_splv(wavelet_transform1, wavelet_transform2)))
      }
    }
  }
  names(tempResult) <- paste("V", 1:length(tempResult), sep = "")
  result[[k]]<- tempResult
  k<- k+noOfSeconds
}

k<- 1
while(k <= noOfIctalFiles)
{
  print(k+totalNoofInterIctalFiles)
  if(noOfSeconds >1)
  {
    fileName<- paste("D:/Work Related/Kaggle EEG/Volumes/Seagate/seizure_detection/competition_data/clips/Dog_2/Dog_2_ictal_segment_",k:(k+noOfSeconds), sep = "")
  }else
  {
    fileName<- paste("D:/Work Related/Kaggle EEG/Volumes/Seagate/seizure_detection/competition_data/clips/Dog_2/Dog_2_ictal_segment_",k, sep = "")
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
  
  tempResult <- data.frame(NULL)
  for(series1 in 1:16)
  {
    for(series2 in 1:16)
    {
      if(series1 > series2)
      {
        wavelet_transform1 <- dwt(patient7_ictal.data[series1,], n.levels = 8)
        wavelet_transform2 <- dwt(patient7_ictal.data[series2,], n.levels = 8)
        tempResult <-rbind(tempResult, t(calc_splv(wavelet_transform1, wavelet_transform2)))
      }
    }
  }
  names(tempResult) <- paste("V", 1:length(tempResult), sep = "")
  result[[k+totalNoofInterIctalFiles]]<- tempResult
  k<- k+noOfSeconds
}

finalList <- list(NULL)
for(i in 1:200)
{
  finalList[i] <- result[998+i]
}
length(result)
plot(result[[1]]$V1)
temp <- array(unlist(result), c(120,8,1318))
tempgraph<-scatterplot3d(rep(1:120,200),rep(1:200, 120),as.vector(temp[,2,]))
persp3d(1:120,1:200,temp[,2,])
image(1:120,1:200,temp[,2,])
rgl.surface(1:120, 1:200, as.vector(temp[,2,]))

#classify preictal and interictal
trainTest<- createDataPartition(targetDataSet[,1], p = .75, list = FALSE)
trainInput <- temp[1,,trainTest]
testInput <- temp[1,,-trainTest]
trainInput <- t(trainInput)
testInput <- t(testInput)
trainTarget <- targetDataSet2[trainTest]
testTarget <- targetDataSet2[-trainTest]


svmmodel_ictal <- svm(trainInput, trainTarget, kernel = "radial",type = "nu-classification",nu = .2,probability = TRUE)
plot(trainTarget, type = "l")
lines(data.matrix(svmmodel$fitted), type = "l", col = "red")
prediction <- predict(svmmodel_ictal, testInput)
plot(testTarget)
lines(data.matrix(prediction), col = "red")

testResult <- data.frame(NULL)
k<- 1
while(k <= noOfTestFiles)
{
    fileName<- paste("D:/Work Related/Kaggle EEG/Volumes/Seagate/seizure_detection/competition_data/clips/Dog_2/Dog_2_test_segment_",k, sep = "")
  fileName <- paste(fileName, ".mat", sep = "")
  patient7_ictal <- readMat(fileName)
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
  
  tempResult <- data.frame(NULL)
  for(series1 in 1:16)
  {
    for(series2 in 1:16)
    {
      if(series1 > series2)
      {
        wavelet_transform1 <- dwt(patient7_ictal.data[series1,], n.levels = 8)
        wavelet_transform2 <- dwt(patient7_ictal.data[series2,], n.levels = 8)
        tempResult <-rbind(tempResult, t(calc_splv(wavelet_transform1, wavelet_transform2)))
      }
    }
  }
  names(tempResult) <- paste("V", 1:length(tempResult), sep = "")
  isIctal <- predict(svmmodel_ictal, tempResult)
  if(isIctal)
  {
    isEarly <- predict(svmmodel_early, tempResult)
  }
  else
  {
    isEarly <- 0
  }
  testResult <- rbind(testResult, data.frame(clip = paste("Dog_1_test_segment_", k, sep = ""), seizure = isIctal, early = isEarly))
}
