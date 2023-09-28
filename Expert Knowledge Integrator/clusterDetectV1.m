function [Clusters, newI, newIB] = clusterDetectV1(image,Offset,endOff)
global minPoints epsilon
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
            if(max(newI(i,j,1),newI(i,j,3)) == newI(i,j,1)) % ignores red clusters, may be eliminate this

               newIB(i,j,:) = [0,0,0]; 
            end
            
        end
    end
end

Clusters = dbScan(newI);
Clusters2 = dbScan(newIB);
Clusters = [Clusters, Clusters2];