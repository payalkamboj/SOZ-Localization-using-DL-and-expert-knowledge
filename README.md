# SOZ-Localization-using-DL-and-expert-knowledge

The dataset is expected to be in the following structure:
        `/path/to/fmridatasets/        
                  ASUAI_001/
                    MO/
                      report/
                         -f1.png
                         -IC_1_thresh.png
                         -t1.png
                 ASUAI_002/
                     MO/
                      report/
                         -f1.png
                         -IC_1_thresh.png
                         -t1.png
                 ...
                 labels/
                    -ASUAI_001.csv
                    -ASUAI_002.csv`
      
For Noise and Non noise Independent Components (IC) classification, run "run.py". This will preprocess the dataset, create train and test folders in the following structure: 

`/path/to/CNNDataTrainTest/
  train/
    NoiseClass/
      ICImage1.png
    NonNoiseClass/
      ICImage2.png
  test/
    NoiseClass/
      ICImage3.png
    class/2
      ICImage4.png`

The trained model will output the labels for the test dataset as noise and non-noise. 

For expert knowledge integration, run "VoxelInfoNewWhite.m" in MATLAB for image processing and expert features engineering. For final classification of test ICs as Noise, RSN and SOZ,
run "expertKnowledgeIntegrator.m".
