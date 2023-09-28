function Cluster = dbScan(image)
global minPoints epsilon 
minpts = minPoints;

Points = [];
for i = 1:size(image,1)
    for j = 1:size(image,2)
        if(image(i,j,1) > 0 || image(i,j,2) > 0 || image(i,j,3) > 0)
           Points = [Points; [j i]]; 
        end
    end
end
corePoints = [];
remPoints = Points;
coreIndices = [];
for k = 1:size(Points,1)
    neighbor = 0;
    for j = 1:size(Points,1)
       if(k ~= j)
          euDist = sqrt(sum((Points(k,:) - Points(j,:)).^2 ));
          if(euDist <= epsilon)
              neighbor = neighbor+1;
              
          end
           
       end
        
    end
    if(neighbor > minpts)
        corePoints = [corePoints; Points(k,:)];
        coreIndices = [coreIndices k];
    end
    
    
end
remPoints(coreIndices,:) = [];
Cluster{1} = [];
clusNum = 1;
remID = 1:size(corePoints,1);
indexR = [];
if(~isempty(remID))
    k = remID(1);
    Cluster{clusNum} = corePoints(k,:);
    indexR = 1;
    indexClus = 1;
    while(~isempty(remID))

        for j = 1:size(corePoints(remID,:),1)

              euDist = sqrt(sum((Cluster{clusNum}(indexClus,:) - corePoints(remID(j),:)).^2 ));
              if(euDist ~= 0)
                  if(euDist <= epsilon)
                      Cluster{clusNum} = [Cluster{clusNum}; corePoints(remID(j),:)];
                      indexR = [indexR j];

                  end
              end


        end

        %Cluster{clusNum} = [Cluster{clusNum}; corePoints(k,:)];

        if(size(indexR,2)<=1 && indexClus == size(Cluster{clusNum},1))
           clusNum = clusNum + 1;

           remID(indexR) = [];
           indexR = [];
           if(~isempty(remID))
                k = remID(1);
                indexClus = 1;
                indexR = 1;
                Cluster{clusNum} = corePoints(k,:);
           end
        else
           remID(indexR) = [];
           indexR = []; 
           indexClus = indexClus + 1;
        end




    end
    %% Border Points
    borderIndex = [];
    for i = 1:size(remPoints,1)
       for k = 1:clusNum-1
           for p = 1:size(Cluster{k},1)
                euDist = sqrt(sum((remPoints(i,:) - Cluster{k}(p,:)).^2 ));
                if(euDist <= epsilon)
                   Cluster{k} = [Cluster{k}; remPoints(i,:)];
                   borderIndex = [borderIndex, i];
                   k = clusNum;
                   break;
                end
           end
       end
    end
    remPoints(borderIndex,:) = [];
end

%% Displaying
 ClusI = image;
 if(~isempty(Cluster))
 for i = 1:size(Cluster,2)
     for k = 1:size(Cluster{i},1)
        ClusI(Cluster{i}(k,1),Cluster{i}(k,2),:) = [255,0,0]; 
     end
 end
 end

imshow(ClusI);

