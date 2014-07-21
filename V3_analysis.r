
library("fNonlinear")
library("R.matlab")
library("tseriesChaos")
library("scatterplot3d")
library("crqa")
library("fractal")
library("fractaldim")
library("tseriesChaos")


seriesNo<- 2
noOfIctalFiles <- 170
noOfInterictalFiles <-1145
frequency <- 400
totalNoofInterIctalFiles <- 1148
noOfSeconds <- 1
nooffiles <- noOfIctalFiles + noOfInterictalFiles

result <- data.frame(NULL)

startValue <- totalNoofInterIctalFiles - noOfInterictalFiles
k<- startValue
while(k < totalNoofInterIctalFiles)
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
  #package fractal
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
  
  for(i in 1:(noOfSeconds*frequency))
  {
    write(patient7_ictal.data[seriesNo,i], "output.dat", append =TRUE , sep = "\n")
  }
  output.theiler <-50
  lyapcommand <- paste("\"D:/Work\ Related/Kaggle\ EEG/Tisean/Tisean_3.0.0/bin/lyap_k\" output.dat -M",output.dimension, " -d", output.lag, " -t", output.theiler, sep = "")
  lyapData<- shell (lyapcommand,intern = T )
  stringIndex <- 14
  output.lyap <- as.numeric(strsplit(lyapData[stringIndex], "= ")[[1]][2])
  output2.fractal <- fd.estim.boxcount(patient7_ictal.data[seriesNo,], plot.loglog = FALSE)
  output.fractal <- output2.fractal$fd
  #separationPlot(patient7_ictal.data[seriesNo,], output.dimension, output.lag, 100)
  #output.dimension <- output.dimension -1
  output.cordim <- C2(patient7_ictal.data[seriesNo,], m = output.dimension, d = output.lag, t = output.theiler, eps = 10)
  result <- rbind(result, data.frame(index = k,fractal = output.fractal, lyap = output.lyap,corrdim =  output.cordim))
  unlink("output.dat")
  k<- k+noOfSeconds
}

k <-1
while(k <= noOfIctalFiles)
{
  print(k)
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
  #package fractal
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
  
  for(i in 1:(noOfSeconds*frequency))
  {
    write(patient7_ictal.data[seriesNo,i], "output.dat", append =TRUE , sep = "\n")
  }
  output.theiler <-50
  lyapcommand <- paste("\"D:/Work\ Related/Kaggle\ EEG/Tisean/Tisean_3.0.0/bin/lyap_k\" output.dat -M",output.dimension, " -d", output.lag, " -t", output.theiler, sep = "")
  lyapData<- shell (lyapcommand,intern = T )
  stringIndex <- 14
  output.lyap <- as.numeric(strsplit(lyapData[stringIndex], "= ")[[1]][2])
  output2.fractal <- fd.estim.boxcount(patient7_ictal.data[seriesNo,], plot.loglog = FALSE)
  output.fractal <- output2.fractal$fd
  #separationPlot(patient7_ictal.data[seriesNo,], output.dimension, output.lag, 100)
  #output.dimension <- output.dimension -1
  output.cordim <- C2(patient7_ictal.data[seriesNo,], m = output.dimension, d = output.lag, t = output.theiler, eps = 10)
  result <- rbind(result, data.frame(index = k,fractal = output.fractal, lyap = output.lyap,corrdim =  output.cordim))
  unlink("output.dat")
  k <- k+noOfSeconds
}


