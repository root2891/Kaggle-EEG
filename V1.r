
library("fNonlinear")
library("R.matlab")
library("timeSeries")
library("scatterplot3d")
library("crqa")
library("fractal")
library("fractaldim")
library("RTisean")

seriesNo<- 1
nooffiles <- 10
fileList<- paste("D:/Work Related/Kaggle EEG/Volumes/Seagate/seizure_detection/competition_data/clips/Patient_7/Patient_7_ictal_segment_",1:nooffiles, sep = "")
fileList <- paste(fileList, ".mat", sep = "")
patient7_ictal <-lapply(fileList,readMat)
patient7_ictal.data <- patient7_ictal[[1]]$data
for(i in 2:nooffiles)
{
  patient7_ictal.data <- cbind(patient7_ictal.data,patient7_ictal[[i]]$data)
}
patient7_ictal.data <- as.ts(patient7_ictal.data)
plot(patient7_ictal.data[1,], type = 'l')

#filtering
output2.filterseries<- medianFilter(patient7_ictal.data[seriesNo,])
plot(output2.filterseries, type = 'l')

#package fractal
output2.mutual <- timeLag(patient7_ictal.data[seriesNo,],method = "mutual",plot.data = TRUE)
output2.dimension <- FNS(patient7_ictal.data[seriesNo, ], dimension = 10)
output2.dimension
output2.cordim <- corrDim(patient7_ictal.data[seriesNo,])
plot(output2.cordim)
eda.plot(output2.cordim)
output2.info<- infoDim(patient7_ictal.data[seriesNo,], dimension = 10)
output2.determinism <- determinism(patient7_ictal.data[seriesNo,], dimension = 5)
output2.series <- embedSeries(patient7_ictal.data[seriesNo,], dimension = 6, tlag = 127)
plot(output2.series)
output2.lyapunov <- lyapunov(patient7_ictal.data[seriesNo,], dimension = 6, scale = c(1,2,4, 8, 16, 32, 64, 128, 256, 512, 1024, 2048 ))
plot(output2.lyapunov)
output2.spacetime <- spaceTime(patient7_ictal.data[seriesNo,], dimension = 6, olag.max = 1000)
plot(output2.spacetime)
#package fractaldim
output2.fractal <- fd.estim.boxcount(patient7_ictal.data[seriesNo,], plot.loglog = TRUE)

#package tserieschaos
output1.mutual <- mutual(patient7_ictal.data[1,], lag.max = 100)
output1.dimention <- false.nearest(patient7_ictal.data[1,], m = 10, 127, 250)
plot.false.nearest(output1.dimention, type = "h")
output1.reccurence <- recurrencePlot(patient7_ictal.data[1,], 6, 45, 5000, eps=sd(patient7_ictal.data[1,])/10)
output1.embedd <- embedd(patient7_ictal.data[1,],2,127)
scatterplot3d(output1.embedd, type= "l")
separationPlot(patient7_ictal.data[seriesNo,],3,45, 1000)

#package fnonlinier
mutualPlot(patient7_ictal.data[seriesNo,])
falsennPlot(patient7_ictal.data[seriesNo,], 10, 127, 180)
lyapunovPlot(patient7_ictal.data[seriesNo,],6, 127, 150, 1000,50, 5)



#modeling
