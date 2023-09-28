function [y, yNum, OuterContour] = peripheryNoiseDetect(Clusters,image,epsilon, perc1, perc2, imContour,numRow, numCol, i)
%evalutaes overlap of the cluster to the brain boundary


%Confirm = 
%y=0 if there is no intersection found in periphery and clusters
%y=-1, if intersection goes over some threshold -> Periphery Noise
%y=1, if there is some intersection
%Perinum: Pixel overlapping

%yNum = ?

y = [];
yNum = [];
if(~isempty(Clusters))
    [M,c] = contour(rgb2gray(imContour));
    perLevel = M(1,1);
    contourData = M;
    index = 1;
    k = 1;
    while(index <= size(M,2))
        if(M(1,index) == perLevel)
            outerPot{k} = M(:,index+1:M(2,index)+index);
            k = k + 1;
        end
        index = index + M(2,index)+1;
    end
%     maxSize = 0;
%     for i = 1:k-1
%        if(maxSize < size(outerPot{i},2))
%            maxSize = size(outerPot{i},2);
%            maxID = i;
%        end
%            
%         
%     end

    for i = 1:k-1
        yM(i) = 0;
        for j = 1:k-1
            if(i~=j)
                p = inpoly2(outerPot{i}',outerPot{j}');
                if(sum(p) > (perc1/100)*size(outerPot{i},2))
                   yM(i) = 1; 
                   break; 
                end
            end
        end
        
    end
    
    OuterContour = [];
    
    kO = 1;
    outerContourSize = 0;
     for i = 1:k-1
         if(yM(i) == 0)
            OuterContour{kO} = outerPot{i};
            outerContourSize = outerContourSize + size(outerPot{i},2);
            kO = kO + 1;
         
         end
         
         
     end
    
    
    %OuterContour = outerPot{maxID}; % very sensitive to dataset, hardcoded value please change this
    
    for i = 1:size(Clusters,2)
        y(i) = 0;
        yNum(i) = 0;
        for k = 1:size(Clusters{i},1)
            minDist = 100;
            for g = 1:size(OuterContour,2)
                for m = 1:size(OuterContour{g},2)
                    %calculate euclidean dist from all 
                    %clusters to outer contour
                    euDist = ((OuterContour{g}(1,m) - Clusters{i}(k,1))^2+(OuterContour{g}(2,m) - Clusters{i}(k,2))^2)^0.5;
                    
                    %get the minimum distance
                    if(euDist < minDist)
                        minDist = euDist;
                    end
                end
            end
            
            if(minDist < epsilon) %if min. distance between cluster pixel and outercontour pixel is less than some threshold
                y(i) = y(i)+1; %there is some intersection
            end    
            
        end
        
        if(~isempty(Clusters{i}))
            trueVal = 0;
            for g = 1:size(OuterContour,2)
                inP = inpoly2(Clusters{i},OuterContour{g}');
                if(100*(sum(inP)/size(Clusters{i},1)) <= perc1)%if ... is less than 60%
                   
                   yNum(i) = y(i); %add y(i) which tells pixel overlapping
                else %if ... greater than 60%
                    yNum(i) = 100*y(i)/size(Clusters{i},1);
                    trueVal = 1;
                    break;
                end
            end
            if(trueVal ~= 1)
                y(i) = -1;
            end
        end
        
        
        if(y(i) ~= -1)
            yNum(i) = 100*y(i)/size(Clusters{i},1);
        else
           yNum(i) = y(i); 
        end
        if(y(i) > (perc2/100)*size(Clusters{i},1) || y(i) > (perc2/100)*outerContourSize)
           y(i) = -1; 
            
        end
    end
else
    
    
end
            