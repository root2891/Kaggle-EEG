library("fNonlinear")
library("R.matlab")
library("tseriesChaos")
library("fractal")
library("fractaldim")
library("tseriesChaos")
library("fields")
library("FNN")
library("caret")
library("RSNNS")
library("e1071")
library("wavelets")
library("scatterplot3d")
library("combinat")
library("rasclass")

path <- "C:/Shubham/Volumes/Seagate/seizure_detection/competition_data/clips/"
patients <- c("Dog_1","Dog_2","Dog_3","Dog_4", "Patient_1","Patient_2","Patient_3","Patient_4","Patient_5","Patient_6","Patient_7","Patient_8")
fileTypes_list <- c("ictal", "interictal", "test")
noOfSeconds <- 1

get_no_of_files <- function(path, patient, fileTypes_list)
{
  noOfFiles_list <- list(NULL)
  for(filetype in fileTypes_list)
  {
    fileList <- list.files(paste(path, patient,"/", sep = ""), pattern =paste("_",filetype,sep = ""))
    noOfFiles_list[filetype] <-  length(fileList)
  }
  return (noOfFiles_list)
}

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

calculation_per_combination <- function(inputData.data)
{
  tempResult <- data.frame(NULL)
  for(series1 in 1:16)
  {
    for(series2 in 1:16)
    {
      if(series1 > series2)
      {
        wavelet_transform1 <- dwt(inputData.data[series1,], n.levels = 8)
        wavelet_transform2 <- dwt(inputData.data[series2,], n.levels = 8)
        tempResult <-rbind(tempResult, t(calc_splv(wavelet_transform1, wavelet_transform2)))
      }
    }
  }
  names(tempResult) <- paste("V", 1:length(tempResult), sep = "")
  return (tempResult)
}

get_data <- function(patient, segment_type, noOfSeconds, k)
{
  if(noOfSeconds >1)
  {
    fileName<- paste(path,patient,"/",patient,"_",segment_type,"_segment_",k:(k+noOfSeconds), sep = "")
  }else
  {
    fileName<- paste(path, patient,"/",patient,"_",segment_type,"_segment_",k, sep = "")
  }
  fileName <- paste(fileName, ".mat", sep = "")
  inputData <-lapply(fileName,readMat)
}

results <- data.frame(NULL)
for(patient in patients)
{
  totalNoofInterIctalFiles <- 1148
  noOfFiles_list <- get_no_of_files(path, patient, fileTypes_list)
  noOfIctalFiles <- noOfFiles_list$ictal
  noOfInterictalFiles <-noOfFiles_list$interictal
  noOfTestFiles <- noOfFiles_list$test
  noOfFiles <- noOfIctalFiles + noOfInterictalFiles
  
  #output data set of the splv data
  splv_output <- list(NULL)
  #output data set containing if file is interictal, earlyictal or ictal
  targetDataSet <- data.frame(NULL)
  splv_result  <- list()
  for(filetype in fileTypes_list[1:2])
  {
    k<-1
    while(k <= noOfFiles_list[filetype] )
    {
      print(paste("patient:",patient, "filetype:", filetype, "K:", k, spe = ""))
      
      #get the data in inputData.data
      inputData <- get_data(patient, filetype, noOfSeconds,k)
      
      
      #create input data set for processing
      inputData.data <- inputData[[1]]$data
      
      frequency <- ncol(inputData.data)
      noOfRods <- nrow(inputData.data)
      
      if(noOfSeconds >1)
      {
        for(i in 2:noOfSeconds)
        {
          inputData.data <- cbind(inputData.data,inputData[[i]]$data)
        }
      }
      
      #create the output data set only in case of ictal or interictal files
      if(filetype == "ictal")
      {
        inputData.latency <- inputData[[1]]$latency
        if(inputData.latency < 15)
        {
          targetDataSet <- rbind(targetDataSet, 1)
        }
        else
        {
          targetDataSet <- rbind(targetDataSet, 2)
        }
      }else
      {
        targetDataSet <- rbind(targetDataSet, 0)
      }
      
      inputData.data <- as.ts(inputData.data) 
      splv_result[[length(splv_result)+1]]<- calculation_per_combination(inputData.data)
      k<- k+noOfSeconds
    }
  }
  noOfCombinations <- ncol(combn(noOfRods, 2))
  write(unlist(splv_result), file = paste("splv_", patient,".txt", sep = "", ncolumns = 1))
  splv_result_array <- array(unlist(splv_result), c(noOfCombinations, 8, noOfFiles))
  
temp2 <- readRasterFolder("D:/Work Related/Kaggle EEG/files/", samplename = "temp_1")


output <- classifyRasclass(temp2, method ='neuralNetwork' )
  
  k<- 1
  while(k<=noOfTestFiles)
  {
    inputData <- get_data(patient, "test", noOfSeconds,k)
    inputData.data <- inputData[[1]]$data
    inputData.data <- as.ts(inputData.data)    
    splv_result<- calculation_per_combination(inputData.data)
    count.zero <- 0
    count.one <- 0
    count.two <- 0
    for(i in 1:noOfCombinations)
    {
      prediction <- predict(model_list[[i]], splv_result[i,1:7])
      if(prediction == 0)
      {
        count.zero <- count.zero+1
      }else if(prediction == 1)
      {
        count.one <- count.one+1
      }else if(prediction == 2)
      {
        count.two <- count.two+1
      }
    }
    countVector <- c(count.zero, count.one, count.two)
    maxValue <- max(countVector)
    value <- which(countVector == maxValue)
    if(value[1] == 1)
    {
      seizure <- 0
      early <- 0
    }else if(value[1] == 2)
    {
      seizure <- 1
      early <- 1  
    }else if(value[1] == 3)
    {
      seizure <- 1
      early<- 0
    }
    results <- rbind(results, data.frame(clip = paste(patient, "_test_Segment_", i, ".mat", sep = ""), seizure =seizure, early =early ))
    print(paste(patient, "seizure:", seizure,"early:",early))
    k<- k+noOfSeconds
  }
  save.image(paste("C:/Shubham/", patient,".RData", sep = ""))
}

mysample <- c(rep(rep(c(1,2), each = 25), 25), rep(rep(c(3,4), each = 25), 25))
mysample <- mysample + sample(c(0, NA), 2500, replace = TRUE, prob = c(1, 10))
myvar1 <- rep(1:50, each = 50) + rnorm(2500, 0, 5)
myvar2 <- rep(rep(1:50), 50) + rnorm(2500, 0, 5)
newdata <- data.frame(mysample, myvar1, myvar2)
object <- new('rasclass')
object <- setRasclassData(newdata, ncols = 50, nrows = 50,
                          xllcorner = 0, yllcorner = 0, cellsize = 1, NAvalue = -9999,
                          samplename = 'mysample')

