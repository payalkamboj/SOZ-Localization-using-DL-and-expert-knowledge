function [Clusters,clusterWidth, peri, periNum, art, whiteM, maxID, clusOverStat, percWhite] = clusterDetect(image,Offset,endOff,imContour,imagePixels,numRow, numCol, index, sizeIM)
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
factorX = 1;
maxThresh = min((50/factorX),90); %(5/2)*50;
maxThresh2 = 20/factorX;
thresh2 = 99;%maxThresh*(sizeIM(index)/max(sizeIM));
threshNew = 60;%100 - maxThresh2*sizeIM(index)/max(sizeIM);

%NOISE DETECTION

%Periphery Noise detection
[peri, periNum, OuterContour] = peripheryNoiseDetect(Clusters,image,2,threshNew,thresh2,imContour,numRow,numCol,i);

%confirm if grayinside is gray matter presence in a slice or just pixels inside outer contour: 
grayInside = 0;
for ko = 1:size(OuterContour,2)
    pInside = inpoly2(imagePixels, OuterContour{ko}');
    grayInside = grayInside + sum(pInside);
end
art = arteryNoiseDetect(Clusters,image,2,40/factorX,imContour);

%percentage here is 40
%Whitem =1 means there is white matter noise, 0 means there is no white
%matter noise
[whiteM, clusOverStat, whitePerc] = whiteMatterNoiseDetectV2(Clusters,image,epsilon, 40/factorX,imContour,imagePixels);
%% Check extra constraints, if symmetric accross the center line then not white matter noise
% if towards the bottom 30% of the image then not white matter

%%I have deleted one if statement because if was not getting executed

%difference betwenn whitePerc and percwhite????
percWhite = 100*whitePerc/grayInside;
clusterWidth = [];
maxID = 0;
maxSize = 0;
for i = 1:size(Clusters,2)
    if(maxSize < size(Clusters{i},1))
        maxSize = size(Clusters{i},1);
        maxID = i;
    end
    
end
clusterWidth = clusterWidthComputation(newI,newIB,Clusters);
%size(Clusters)
%clusterWidth
imshow(newI)
%[M,c] = contour(rgb2gray(newI))