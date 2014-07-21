import numpy as np 
import neurolab as nl 

frequency  = 400
noOfIctalFiles = 170
noOfInterictalFiles = 1148
totalFiles = noOfIctalFiles +noOfInterictalFiles

noOfTrainValues = 100

totalValues = totalFiles*frequency

data = np.loadtxt("inputData.txt", delimiter = " ", skiprows = 1)

trainData = []
outputValue = []

for(i in 0:(totalFiles -noOfTrainValues -2))
{
	print "i = %d" % i
	tempTrainData = []
	for(j in 0:noOfTrainValues)
	{
		print "j = %d" % j
		tempTrainData[j] <- data[i+j]
	}
	index = i+100
	trainData.append(tempTrainData)
	outputValue[i] =  inputData[1,index]
}