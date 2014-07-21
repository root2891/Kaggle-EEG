library("fNonlinear")
library("R.matlab")
library("tseriesChaos")
library("fractal")
library("fractaldim")
library("tseriesChaos")
library("fields")
library("FNN")
library("synchrony")

noOfNeighbors<-5
noOfSeconds <- 5
seriesNo<- 2
noOftrainingValues <-250
noOfIctalFiles <- 160
noOfInterictalFiles <-410
frequency <- 400
totalNoofInterIctalFiles <- 418
nooffiles <- noOfIctalFiles + noOfInterictalFiles

synchrony <- data.frame(NULL)
startValue <- totalNoofInterIctalFiles - noOfInterictalFiles
k<- startValue
while(k < totalNoofInterIctalFiles)
{
  print(k)
  fileName<- paste("D:/Work Related/Kaggle EEG/Volumes/Seagate/seizure_detection/competition_data/clips/Dog_1/Dog_1_interictal_segment_",k:(k+noOfSeconds), sep = "")
  fileName <- paste(fileName, ".mat", sep = "")
  patient7_ictal <-lapply(fileName,readMat)
  patient7_ictal.data <- patient7_ictal[[1]]$data
  for(i in 2:noOfSeconds)
  {
    patient7_ictal.data <- cbind(patient7_ictal.data,patient7_ictal[[i]]$data)
  }
  patient7_ictal.data <- as.ts(patient7_ictal.data)
  plot(patient7_ictal.data[seriesNo,], type = 'l')
  
  wavelet<- dwt(patient7_ictal.data[seriseNo,])
  
  tempsynchrony<- data.frame(NULL)
  for(series1 in 1:16)
  {
    for(series2 in 1:16)
    {
      if(series1 != series2)
      {
        #embeddedSeries1 <- embedd(patient7_ictal.data[series1,], output.dimension, output.lag)
        #embeddedSeries2 <- embedd(patient7_ictal.data[series2,], output.dimension, output.lag)
        
        sync<- phase.sync(patient7_ictal.data[series1,] ,patient7_ictal.data[series2,])
        tempsynchrony <- rbind(tempsynchrony, sync$pval)
      }
    }
  }
  tempsynchrony<- t(tempsynchrony)
  names(tempsynchrony) <-  paste("sync", 1:length(tempsynchrony), sep = "")
  synchrony <- rbind(synchrony, tempsynchrony)
  k<- k+noOfSeconds
}


k<-1
while(k < noOfIctalFiles)
{
  print(k)
  fileName<- paste("D:/Work Related/Kaggle EEG/Volumes/Seagate/seizure_detection/competition_data/clips/Dog_1/Dog_1_ictal_segment_",k:(k+noOfSeconds), sep = "")
  fileName <- paste(fileName, ".mat", sep = "")
  patient7_ictal <-lapply(fileName,readMat)
  patient7_ictal.data <- patient7_ictal[[1]]$data
  for(i in 2:noOfSeconds)
    
  {
    patient7_ictal.data <- cbind(patient7_ictal.data,patient7_ictal[[i]]$data)
  }
  patient7_ictal.data <- as.ts(patient7_ictal.data)
  plot(patient7_ictal.data[seriesNo,], type = 'l')

  tempsynchrony<- data.frame(NULL)
  for(series1 in 1:16)
  {
    for(series2 in 1:16)
    {
      if(series1 != series2)
      {
        
        sync<- phase.sync(patient7_ictal.data[series1,], patient7_ictal.data[series2,])
        tempsynchrony <- rbind(tempsynchrony, sync$pval)
      }
    }
  }
  tempsynchrony<- t(tempsynchrony)
  names(tempsynchrony) <- paste("sync", 1:length(tempsynchrony), sep = "")
  synchrony <- rbind(synchrony, tempsynchrony)
  k<- k+noOfSeconds  
}

row.names(result) <- paste("record",1:nrow(result),sep = "")
for(i in 1:114)
{
  plot(result[,i], type = "l")
}
