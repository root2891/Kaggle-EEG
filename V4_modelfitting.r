library("R.matlab")
library("RSNNS")
library("robfilter")
library("fractal")
library("caret")

seriesNo<- 1
frequency <- 400
totalNoofInterIctalFiles <- 418

trainData <- data.frame(NULL)
timeToSeizure <- data.frame(NULL)
for(i in 1:totalNoofInterIctalFiles)
{
  interictalfileList<- paste("D:/Work Related/Kaggle EEG/Volumes/Seagate/seizure_detection/competition_data/clips/Dog_1/Dog_1_interictal_segment_",i, sep = "")
  interictalfileList <- paste(interictalfileList, ".mat", sep = "")
  interictalData <- readMat(interictalfileList)
  trainData <- rbind(trainData, interictalData$data[seriesNo,])
  timeToSeizure <- rbind(timeToSeizure,(totalNoofInterIctalFiles-i))
}


names(trainData) = c(paste("V",1:frequency,sep = ""))
names(timeToSeizure) = c("timeToseizure")

jmodel <- jordan(trainData, timeToSeizure, size = c(10))
plotIterativeError(jmodel)
plotRegressionError(timeToSeizure$timeToseizure, jmodel$fitted.values)
plot(timeToSeizure$timeToseizure, type="l")
lines(jmodel$fitted.values , col = "red")

emodel <- elman(trainData, timeToSeizure, size = c(10))
plotIterativeError(emodel)
plot(timeToSeizure$timeToseizure, type="l")
lines(emodel$fitted.values , col = "red")
