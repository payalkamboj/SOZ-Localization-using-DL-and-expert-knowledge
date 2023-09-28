function [Clusters, whiteM, clusOverStat,whiteMatter] = clusterDetection(image,Offset,endOff,imContour,imagePixels,numRow, numCol, sizeIM)
global minPoints epsilon factorM
newI = [];
for i = Offset:size(image,1)-endOff
    for j = Offset:size(image,2)-endOff
        newI(i,j,:) = image(i,j,:);
        
    end
end
newIB = newI;

for i = 1:size(newI,1)
    for j = 1:size(newI,2)
        if(newI(i,j,1) == newI(i,j,2) && newI(i,j,2)== newI(i,j,3)) % ignores white and gray space
           newI(i,j,:) = [0,0,0]; 
           newIB(i,j,:) = [0,0,0];
        else
            if(max(newI(i,j,1),newI(i,j,3)) == newI(i,j,3)) % ignores blue clusters, may be eliminate this

               newI(i,j,:) = [0,0,0]; 
            end
            if(max(newIB(i,j,1),newIB(i,j,3)) == newIB(i,j,1)) % ignores red clusters, may be eliminate this

               newIB(i,j,:) = [0,0,0]; 
            end
            
        end
    end
end


%Cluster detection
Clusters = dbScan(newI);
Clusters2 = dbScan(newIB);
Clusters = [Clusters, Clusters2]; %number of red and blue clusters ?
%[M,c] = contour(rgb2gray(imContour));
factorX = 1;
maxThresh = min((50/factorX),90); %(5/2)*50;
maxThresh2 = 20/factorX;
thresh2 = 99;%maxThresh*(sizeIM(index)/max(sizeIM));
threshNew = 60;%100 - maxThresh2*sizeIM(index)/max(sizeIM);


%NOISE DETECTION
%

%confirm if grayinside is gray matter presence in a slice or just pixels inside outer contour: 


%Whitem =1 means there is white matter noise, 0 means there is no white
%matter noise
[whiteM, clusOverStat, whitePerc,whiteMatter] = whiteMatterNoiseDetectV2(Clusters,image,epsilon, 40/factorX,imContour,imagePixels);
%% Check extra constraints, if symmetric accross the center line then not white matter noise
% if towards the bottom 30% of the image then not white matter
imshow(newI)
hold on;
contour(rgb2gray(image));
holdoff;
%[M,c] = contour(rgb2gray(newI))