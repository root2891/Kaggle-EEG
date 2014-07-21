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
noOfSeconds <- 1
seriesNo<- 2
noOftrainingValues <-250
noOfIctalFiles <- 171
noOfInterictalFiles <-1147
frequency <- 400
totalNoofInterIctalFiles <- 1148
nooffiles <- noOfIctalFiles + noOfInterictalFiles

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
#synchrony <- data.frame(NULL)
startValue <- totalNoofInterIctalFiles - noOfInterictalFiles
k<- startValue
while(k < totalNoofInterIctalFiles)
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
  output2.dimension <- FNN(patient7_ictal.data[seriesNo, ], dimension = noOfDimension)
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
        remove("embeddedSeries1","embeddedSeries2","neighbour1","neighbour2")
        #sync<- phase.sync(embeddedSeries1, embeddedSeries2)
        #tempsynchrony <- rbind(tempsynchrony, sync$pval)
      }
    }
  }
  tempResult <- t(tempResult)
 # tempsynchrony<- t(tempsynchrony)
  names(tempResult) <- paste("V", 1:length(tempResult), sep = "")
  #names(tempsynchrony) <- <- paste("sync", 1:length(tempResult), sep = "")
  result <- rbind(result, tempResult)
  #synchrony <- rbind(synchrony, tempsynchrony)
  k<- k+noOfSeconds
}

length(result)

k<-1
while(k < noOfIctalFiles)
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
  output2.dimension <- FNN(patient7_ictal.data[seriesNo, ], dimension = noOfDimension)
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
        remove("embeddedSeries1","embeddedSeries2","neighbour1","neighbour2")
        #sync<- phase.sync(embeddedSeries1, embeddedSeries2)
        #tempsynchrony <- rbind(tempsynchrony, sync$pval)
      }
    }
  }
  tempResult <- t(tempResult)
  #tempsynchrony<- t(tempsynchrony)
  names(tempResult) <- paste("V", 1:length(tempResult), sep = "")
  #names(tempsynchrony) <- <- paste("sync", 1:length(tempResult), sep = "")
  result <- rbind(result, tempResult)
  #synchrony <- rbind(synchrony, tempsynchrony)
  k<- k+noOfSeconds  
}

row.names(result) <- paste("record",1:nrow(result),sep = "")
for(i in 1:99)
{
  name <- paste("number",i)
  plot(result[,i], type= "b", ylab = name)
}
