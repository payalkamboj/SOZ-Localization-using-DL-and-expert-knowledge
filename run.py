'''Run this file for the DL part of noise and non noise classification'''

import preprocessing
import NoiseClassificationCNN
import readNoiseFiles

if __name__ == '__main__':
  
    path='../../PCHData/fmridatasets/'
    subjectsID = 'ASUAI_0*'
    datasetHome = 'CNNDataTrainTest/'
    NoiseICLabels = 'CNNNoiselabels.txt'
    preprocess = preprocessing.PreprocessData(path,subjectsID)

    #Noise- non noise classification using DL
    NoiseNonNoiseDLClassification = NoiseClassificationCNN.DLNoiseCLassification(datasetHome)

    TrainingTestingDataExits = 0
    if TrainingTestingDataExits == 0:
        NoiseNonNoiseDLClassification.CreateTrainTestDataset()
    NoiseNonNoiseDLClassification.TrainingTestingCNN()

    readNoiseFiles.extractNoiseLabels(NoiseICLabels)
  


