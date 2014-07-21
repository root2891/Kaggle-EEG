
library("fNonlinear")
library("R.matlab")
library("tseriesChaos")
library("scatterplot3d")
library("crqa")
library("fractal")
library("fractaldim")
library("tseriesChaos")


seriesNo<- 1
nooffiles <- 60
result <- data.frame(NULL)
i<-1
while(k < nooffiles)
{
  fileName<- paste("D:/Work Related/Kaggle EEG/Volumes/Seagate/seizure_detection/competition_data/clips/Dog_1/Dog_1_ictal_segment_",k:k+5, sep = "")
  fileName <- paste(fileName, ".mat", sep = "")
  patient7_ictal <-lapply(fileList,readMat)
  patient7_ictal.data <- patient7_ictal[[1]]$data
  for(i in 2:5)
  {
    patient7_ictal.data <- cbind(patient7_ictal.data,patient7_ictal[[i]]$data)
  }
  patient7_ictal.data <- patient7_ictal$data
  patient7_ictal.data <- as.ts(patient7_ictal.data)
  plot(patient7_ictal.data[1,], type = 'l')
  output2.filterseries<- medianFilter(patient7_ictal.data[seriesNo,])
  plot(output2.filterseries, type = 'l')
  #package fractal
  noOfDimension <-10
  output2.lag <- timeLag(patient7_ictal.data[seriesNo,],method = "acfdecor",plot.data = TRUE)
  output.lag <- output2.lag[1]
  output2.dimension <- FNN(patient7_ictal.data[seriesNo, ], dimension = noOfDimension)
  for(j in 1:noOfDimension)
  {
    if(output2.dimension[3,j] == 0)
    {
      output.dimension <-j
      break
    }
  }
  
  for(i in 1:400)
  {
    write(patient7_ictal.data[seriesNo,i], "output.dat", append =TRUE , sep = "\n")
  }
  lyapcommand <- paste("\"D:/Work\ Related/Kaggle\ EEG/Tisean/Tisean_3.0.0/bin/lyap_k\" output.dat -M",output.dimension, " -d", output.lag, " -t", output.theiler, sep = "")
  lyapData<- shell (lyapcommand,intern = T )
  stringIndex <- 9+ output.dimension
  output.lyap <- as.numeric(strsplit(lyapData[stringIndex], "= ")[[1]][2])
  output2.fractal <- fd.estim.boxcount(patient7_ictal.data[seriesNo,], plot.loglog = FALSE)
  output.fractal <- output2.fractal$fd
  #separationPlot(patient7_ictal.data[seriesNo,], output.dimension, output.lag, 100)
  output.theiler <-50
  #output.dimension <- output.dimension -1
  output2.cordim <- C2(patient7_ictal.data[seriesNo,], m = output.dimension, d = output.lag, t = output.theiler, eps = 10)
  #output2.info<- infoDim(patient7_ictal.data[seriesNo,], dimension = output.dimension)
  #output.infodim <- attr(output2.info, "slope")[as.character(output.dimension)]
  result <- rbind(result, data.frame(index = k,fractal = output.fractal, lyap = output.lyap,corrdim =  output.cordim))
  unlink("output.dat")
}
result["fractal"]
plot(result[,"index"], result[, "corrdim"], type = "l")

