library("R.matlab")
library("RSNNS")
library("robfilter")
library("fractal")

seriesNo<- 1
noOftrainingValues <-100
noOfIctalFiles <- 50
noOfInterictalFiles <-418
frequency <- 400
totalNoofInterIctalFiles <- 418

nooffiles <- noOfIctalFiles + noOfInterictalFiles

if(noOfIctalFiles != 0 )
{
  ictalfileList<- paste("D:/Work Related/Kaggle EEG/Volumes/Seagate/seizure_detection/competition_data/clips/Dog_1/Dog_1_ictal_segment_",1:noOfIctalFiles, sep = "")
  ictalfileList <- paste(ictalfileList, ".mat", sep = "")
  ictalinputData <-lapply(ictalfileList,readMat)
  ictalinputData.data <- ictalinputData[[1]]$data
  for(i in 2:noOfIctalFiles)
  {
    ictalinputData.data <- cbind(ictalinputData.data,ictalinputData[[i]]$data)
  }
}

if(noOfInterictalFiles != 0 )
{
  startValue <- totalNoofInterIctalFiles - noOfInterictalFiles +1
  interictalfileList<- paste("D:/Work Related/Kaggle EEG/Volumes/Seagate/seizure_detection/competition_data/clips/Dog_1/Dog_1_interictal_segment_",startValue:totalNoofInterIctalFiles, sep = "")
  interictalfileList <- paste(interictalfileList, ".mat", sep = "")
  interictalinputData <-lapply(interictalfileList,readMat)
  interictalinputData.data <- interictalinputData[[1]]$data
  for(i in 2:noOfInterictalFiles)
  {
    interictalinputData.data <- cbind(interictalinputData.data,interictalinputData[[i]]$data)
  }
}