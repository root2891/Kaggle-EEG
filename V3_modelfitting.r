library("R.matlab")

inputData <- readMat("D:/Work Related/Kaggle EEG/Volumes/Seagate/seizure_detection/competition_data/clips/Dog_1/Dog_1_ictal_segment_1.mat")
inputData.data <- inputData$data
for(i in 1:400)
{
  write(inputData.data[1,i], "output.dat")
}
trainData <- data.frame(NULL)
outputValue <- data.frame(NULL)
for(i in 1:298)
{
  for(j in 0:99)
  {
    trainData[i,j] <- inputData.data[1,i+j]
  }
  index <- j+100
  outputValue[j,1] <- inputData[1,index]
}

