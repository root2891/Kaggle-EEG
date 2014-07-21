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
library("dtw")
library("parallel")

path <- "D:/Work Related/Kaggle EEG/Volumes/Seagate/seizure_detection/competition_data/clips/"
patients <- c("Patient_5")
fileTypes_list <- c("ictal", "interictal", "test")
noOfSeconds <- 1

v <- dyn.load("D:/Work Related/Kaggle EEG/dtw2.so")
dtw2 <- function(noOfRods, frequency, output)
{
  .C("dtw", noOfRods, frequency,output)
}

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

calculation_per_combination <- function(inputData.data, noOfRods, frequency)
{
  result <- data.frame(NULL)
  write.table(inputData.data, "D:/Work Related/Kaggle EEG/table3.txt", row.names = F, col.names = F)
  outputVal <-1:ncol(combn(noOfRods, 2))
  out <- dtw2(noOfRods,frequency, as.double(outputVal) )
  unlink("C:/Shubham/table3.txt")
  
  return (out[[3]])
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

process_Data <- function(fileNo)
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
      dtw_result<-rbind(dtw_result,calculation_per_combination(inputData.data,noOfRods, frequency))
      k<- k+noOfSeconds
}


results <- data.frame(NULL)
for(patient in patients)
{
  noOfFiles_list <- get_no_of_files(path, patient, fileTypes_list)
  noOfIctalFiles <- noOfFiles_list$ictal
  noOfInterictalFiles <-noOfFiles_list$interictal
  noOfTestFiles <- noOfFiles_list$test
  noOfFiles <- noOfIctalFiles + noOfInterictalFiles
  
  #output data set containing if file is interictal, earlyictal or ictal
  targetDataSet <- data.frame(NULL)
  dtw_result  <- data.frame(NULL)
  noOfRods <- 16
  frequency <- 400
  for(filetype in fileTypes_list[1:2])
  {
    c1<- makeCluster(4)
    parLapply(c1, 1:noOfFiles_list[filetype], process_Data )
    save.image(paste("D:/Work Related/Kaggle EEG/", patient,".RData", sep = ""))
  }
  
  noOfCombinations <- ncol(combn(noOfRods, 2))
  
  trainInput <- dtw_result
  trainTarget <- targetDataSet
  model_list <- svm(trainInput, trainTarget, kernel = "radial", type = "nu-classification", nu = 0.01)
  k<- 1
  
  while(k<=noOfTestFiles)
  {
    inputData <- get_data(patient, "test", noOfSeconds,k)
    inputData.data <- inputData[[1]]$data
    inputData.data <- as.ts(inputData.data)  
    frequency <- ncol(inputData.data)
    noOfRods <- nrow(inputData.data)
    result<- calculation_per_combination(inputData.data,noOfRods,frequency)
    prediction <- predict(model_list, t(result))
    value <- prediction
    if(value == 0)
    {
      seizure <- 0
      early <- 0
    }else if(value== 1)
    {
      seizure <- 1
      early <- 1  
    }else if(value == 2)
    {
      seizure <- 1
      early<- 0
    }
    results <- rbind(results, data.frame(clip = paste(patient, "_test_Segment_", k, ".mat", sep = ""), seizure =seizure, early =early ))
    print(paste(patient, "seizure:", seizure,"early:",early))
    k<- k+noOfSeconds
  }
  save.image(paste("D:/Work Related/Kaggle EEG/", patient,".RData", sep = ""))
}

