clc
clear all
close all
FinalResults = [];
numberSubjects = 52

for  testing_Subject = 1:numberSubjects
    OverallDataNumClusters = [];
    OverallDataAssymetryScore = [];
    OverallDataGiniI = [];
    OverallDataGiniS = [];
    OverallDataWhiteOverlap = [];
    OverallDataLabels = [];
    OverallUSALabels = [];
    OverallAge = [];
    OverallGender = [];
    OverallBigClustersCount=[];
  
    TrainDataNumClusters = [];
    TrainDataAssymetryScore = [];
    TrainDataGiniI = [];
    TrainDataGiniS = [];
    TrainDataLabels = [];
    TrainUSALabels = [];
    TrainDataWhiteOverlap=[];
    TrainAge = [];
    TrainGender = [];
    TrainDataBigCluster=[];
    perc = 80;
    TestDataNumClusters = [];
    TestDataAssymetryScore = [];
    TestDataGiniI = [];
    TestDataGiniS = [];
    TestDataLabels = [];
    TestUSALabels = [];
    TestAge = [];
    TestGender = [];
    A1Size = [];
    A2Size = [];
    A3Size = [];
    MSize = [];
    FSize = [];
  
    folder = dir('../../PCHData/fmridatasets/ASUAI_0*');
    for foldIMKV = 1:numberSubjects
    PatID = [1:numberSubjects];
    foldA = [1:numberSubjects];
    foldI = PatID(foldIMKV);
    fold = foldA(foldI); 

    %workspaces created from 'VoxelInforNewWhite.m'
    load(strcat(['Workspace-' folder(foldIMKV).name 'V4.mat'])); 
    
    folder = dir('../../PCHData/fmridatasets/ASUAI_0*');
    
    files = dir(strcat(['../../PCHData/fmridatasets/' folder(foldIMKV).name '/MO/report/*_thresh*']));
    n = length(files) ;
    listing = dir(strcat(['../../PCHData/fmridatasets/' folder(foldIMKV).name '/MO/report/t*.txt']));
    tic
    numberOfClusters = [];
    numberOfClustersP = [];
    assymetryScore = [];
    GiniI = [];
    GiniS = [];
    for fileNum = 1:n
        filename = files(fileNum).name;
        icaData{fileNum} = dlmread(strcat(['../../PCHData/fmridatasets/' folder(foldIMKV).name '/MO/report/' listing(fileNum).name]));
    
        %% BOLD SIGNAL SPARSITY IN ACTIVELET DOMAIN
        %extract activity related signal component from time-courses
        lengthSWT = 256;
        GEx{fileNum} = [];
        GiniI(fileNum) = 0;
        for g = 1:lengthSWT:size(icaData{fileNum},1)-lengthSWT
            signalB = icaData{fileNum}(g:g+lengthSWT-1,1);
            SWC = swt(signalB,3,'bior2.2');
            GExLine = lasso(SWC',signalB,'Lambda',0.5);
            BN = GiniIndex(GExLine);
            if(GiniI(fileNum) < BN)
                GiniI(fileNum) = BN;
                GEx{fileNum} = GExLine;
            end
    
        end
        
    
        %% Sparsity in sine domain
    
        dictionary2 = {'sin'};
        [mpdict2,nbvect2] = wmpdictionary(length(icaData{fileNum}),'lstcpt',dictionary2);
        y2 = wmpalg('MPALG',icaData{fileNum},mpdict2,'itermax',35);
        GiniS(fileNum) = GiniIndex(y2);
    
    
        %%
           [G,H] = find(ClusterSizePat{fileNum} > pixelLimit);
           numberOfClusters(fileNum) = max(size(G));
    
           disp(foldIMKV)
           whiteClusterOverlap=OverallWhiteFile{fileNum};
           sizeOfCluster=ClusterSizePat{fileNum};
           lk =1;
           whiteOverlapCounter = 0;
           BigClusterCounter = 0;
           while lk <= size(sizeOfCluster,2)
               %Number of cluster count
               if(sizeOfCluster(lk) >= pixelLimit)
                   BigClusterCounter=BigClusterCounter+1;
                   
                   %White matter overlap from grey matter towards ventricles count
                   if(whiteClusterOverlap(lk) == 1)
                       whiteOverlapCounter=whiteOverlapCounter+1;
                   end
                   disp(filename)
               end
               lk = lk + 1;
           end
           WhiteOverlapCount{foldIMKV,fileNum} = whiteOverlapCounter;
           WhiteOverlapCnt(fileNum)=whiteOverlapCounter;
           BigClusterCount{foldIMKV,fileNum} = BigClusterCounter;
           BigClusterCnt(fileNum)=BigClusterCounter;
    
    
         filenameString(fileNum) = sscanf(filename,'IC_%d_thresh');
    end
    
        %%
        [B,I] = sort(filenameString');
        labelV6 = label5(I);
        GiniISorted = GiniI(I);
        sum(labelV6-labelV5)
        WhiteOverlapSorted = WhiteOverlapCnt(I);
        BigClusterCntSorted = BigClusterCnt(I);
        OverallBigClustersCount = [OverallBigClustersCount BigClusterCntSorted]
        OverallDataWhiteOverlap = [OverallDataWhiteOverlap WhiteOverlapSorted];
        OverallDataGiniI = [OverallDataGiniI GiniISorted];
        
        TestDataSVMGiniI{foldIMKV} = GiniISorted;
        TestDataSVMGiniS{foldIMKV} = GiniISorted;
        TestDataSVMWhiteOverlap{foldIMKV} = WhiteOverlapSorted;
        TestDataSVMNum{foldIMKV} = numberOfClusters;
        
        GiniSSorted = GiniS(I);
        OverallDataGiniS = [OverallDataGiniS GiniSSorted];
        
        numberOfClusters = numberOfClusters(I);
        OverallDataNumClusters = [OverallDataNumClusters numberOfClusters];
        TrainDataNumClusters = [TrainDataNumClusters numberOfClusters(1:floor((perc/100)*max(size(GiniSSorted))))];
        TestDataNumClusters{foldIMKV} = numberOfClusters(floor((perc/100)*max(size(GiniSSorted)))+1:end);
        
        if (testing_Subject ~= foldIMKV)
            TrainDataWhiteOverlap = [TrainDataWhiteOverlap WhiteOverlapSorted(1:end)];
            TrainDataBigCluster = [TrainDataBigCluster BigClusterCntSorted(1:end)];
            TrainDataGiniS = [TrainDataGiniS GiniSSorted(1:end)]; 
            TrainDataGiniI = [TrainDataGiniI GiniISorted(1:end)];
            TrainDataLabels = [TrainDataLabels num(:,1)'];
    
        else
            TestDataWhiteOverlap = WhiteOverlapSorted(1:end);
            TestDataBigCluster = BigClusterCntSorted(1:end);
            TestDataGiniS = GiniSSorted(1:end);
            TestDataGiniI = GiniISorted(1:end);
            TestDataLabels = num(:,1)';
        end
    
        labelsOverall = num(:,1)';
           
    end
    
    TrainingData = [TrainDataWhiteOverlap', TrainDataBigCluster', TrainDataGiniS', TrainDataGiniI'];
    TrainLabel = TrainDataLabels';
    [G,H] = find(TrainLabel == 2);
    TrainingData(G,:) = [];
    TrainLabel(G) = [];
    [G,H] = find(TrainLabel == 0);
    TrainingData(G,:) = [];
    TrainLabel(G) = [];
    TrainingDataSmote = [TrainDataWhiteOverlap', TrainDataBigCluster', TrainDataGiniS', TrainDataGiniI' TrainDataLabels'];
    TrainLabel = TrainDataLabels';
    [G,H] = find(TrainLabel == 2);
    TrainingDataSmote(G,:) = [];
    TrainLabel(G) = [];
    [G,H] = find(TrainLabel == 0);
    TrainingDataSmote(G,:) = [];
    TrainLabel(G) = [];
    classLabel = [1 3];
    classLabel = classLabel';
    T = array2table(TrainingDataSmote);
    
    %Applying SMOTE on SOZ expert features
    [newdata,visdata] = smote(T,'3');
    SyntheticGeneratedLabels = newdata{:,5}
    SyntheticData=newdata{:,[1,2,3,4]}
    SyntheticLabels=str2num(char(SyntheticGeneratedLabels));
    NewGeneratedData = vertcat( TrainingData,SyntheticData); 
    newGeneratedLabels = vertcat(TrainLabel,SyntheticLabels);
    model = fitcsvm(NewGeneratedData,newGeneratedLabels,'KernelFunction','linear');
    
    Test = [TestDataWhiteOverlap', TestDataBigCluster', TestDataGiniS', TestDataGiniI'];
    [newLabels, score] = predict(model,Test);
    
    overallStats=[];
    
    Stats = horzcat(TestDataLabels', newLabels, score,TestDataWhiteOverlap', TestDataBigCluster',TestDataGiniS', TestDataGiniI');
    overallStats = vertcat(overallStats, Stats);
    
    trueLabels = overallStats(:,1);
    PredictedLabels = overallStats(:,2);
    
    %Integrating expert knowledge with DL classified Noise/Non noise ICs
    subjectName = folder(testing_Subject).name;
    subjectNumber = regexp(subjectName,'\d*','Match');
    subjectNumber= str2double(subjectNumber);
    
    SubjectNoiseICs_Step1 = dir(strcat(['../../PCHData/DeepLearning/code/PatientLevelTestResults/HybridFinalResults/Noise' num2str(subjectNumber) '.txt']));
    SubjectNoiseICs_Step1 = SubjectNoiseICs_Step1.name;
    subjectFileName = strcat(['../../PCHData/DeepLearning/code/PatientLevelTestResults/HybridFinalResults/' SubjectNoiseICs_Step1]);
    fileID = fopen(subjectFileName,'r');
    NoiseICs = fscanf(fileID,'%s')
    NoiseICs = str2num(NoiseICs);
    
    PredictedLabels(NoiseICs)=0
    [G,H] = find(overallStats(:,4) > 0.9 & overallStats(:,2)==3)
    PredictedLabels(G)=3
    PredictedLabels(setdiff(1:end,NoiseICs)) = 3;
    [Precision, Recall, F1Score, Accuracy, Specificity, TP, TN, FP, FN]=SOZDLMLResults(trueLabels,PredictedLabels)
  
    FinalResults=vertcat( FinalResults, [TP, TN, FP, FN])
    filename = 'ResultsWithoutNumClusters.xlsx';
    writematrix(FinalResults,filename,'Sheet',1);
end
