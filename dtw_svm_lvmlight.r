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
library("dtw")

path <- "D:/Work Related/Kaggle EEG/Volumes/Seagate/seizure_detection/competition_data/clips/"
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

calculation_per_combination <- function(inputData.data, noOfRods)
{
  result <- data.frame(NULL)
  tempResult <- dtwDist(inputData.data)
  for(i in 1:noOfRods)
  {
    for(j in 1:noOfRods)
    {
      if(i > j)
      {
        result <- rbind(result, tempResult[i,j])
      }
    }
  }
  
  return (t(result))
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

write_to_file_train <- function(dtw_result, targetDataSet)
{
  path <- "D:/Work Related/Kaggle EEG/svm/"
  fileLocation <- paste(path, "trainfile.dat", sep = "")
  for(j in 1:nrow(trainInput))
  {
    if(targetDataSet[j,] == 0)
      outputval <- 2
    else
      outputval <- 1
    
    output <- NULL
    for(i in 1:ncol(trainInput))
    {
      output <- paste(output," ",i, ":",  dtw_result[j,i], sep = "")
    }
    val <- paste(outputval, output)
    write(val, fileLocation, append = TRUE)
  }
}

write_to_file_predict <- function(result)
{
  path <- "D:/Work Related/Kaggle EEG/svm/"
  fileLocation <- paste(path, "testfile.dat", spe = "")
  for(j in 1:nrow(trainInput))
  { 
    output <- NULL
    for(i in 1:ncol(trainInput))
    {
      output <- paste(output," ",i, ":",  result[j,i], sep = "")
    }
    write(output, fileLocation, append = TRUE)
  }
}

results <- data.frame(NULL)
for(patient in patients)
{
  noOfFiles_list <- get_no_of_files(path, patient, fileTypes_list)
  noOfIctalFiles <- noOfFiles_list$ictal
  noOfInterictalFiles <-noOfFiles_list$interictal
  noOfTestFiles <- noOfFiles_list$test
  noOfFiles_list$ictal <- 10
  noOfFiles_list$interictal <-10
  noOfFiles_list$test <- 10
  noOfFiles <- noOfIctalFiles + noOfInterictalFiles
  
  #output data set containing if file is interictal, earlyictal or ictal
  targetDataSet <- data.frame(NULL)
  dtw_result  <- data.frame(NULL)
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
      dtw_result<-rbind(dtw_result,calculation_per_combination(inputData.data, noOfRods))
      k<- k+noOfSeconds
    }
  }
  frequency <- ncol(inputData.data)
  noOfRods <- nrow(inputData.data)
  noOfCombinations <- ncol(combn(noOfRods, 2))
  write_to_file_train(dtw_result, targetDataSet)
  trainCommand <- paste("\"D:/Work Related/Kaggle EEG/svm/svm_multiclass_learn.exe\" -c 1.0 trainfile.dat model.dat", sep = "")
  shell(trainCommand, intern = F)
  k<- 1
  result <- data.frame(NULL)
  while(k<=noOfTestFiles)
  {
    inputData <- get_data(patient, "test", noOfSeconds,k)
    inputData.data <- inputData[[1]]$data
    inputData.data <- as.ts(inputData.data)    
    result <- rbind(dtw_result,calculation_per_combination(inputData.data, noOfRods))
  }
  
    write_to_file_predict(result)
    learnCommand <- paste("\"D:/Work\ Related/Kaggle\ EEG/svm/svm_multiclass_classify.exe\" trainfile.dat model.dat, predictfile_",patient,"_", k , sep = "")
    prediction <- predict(model_list, result)
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
    print(paste(patient, "seizure:", seizure,"early:",early))
    k<- k+noOfSeconds
  
  save.image(paste("C:/Shubham/", patient,".RData", sep = ""))
}


