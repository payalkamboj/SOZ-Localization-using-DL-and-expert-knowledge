function [y, clusOverStat, whitePerc,whiteMatter] = whiteMatterNoiseDetectV2(Clusters,image,epsilon, perc, imContour,imagePixels)
y = [];
clusOverStat = [];
epsilon2 = 0.002;

if(~isempty(Clusters))
    [M,c] = contour(rgb2gray(imContour));
    
    %% whitematter region
    [G,H] = find(M(1,:) == c.LevelList(end-1));
    whiteMatter = [];
    i = 1;
    k = 1;
    while (i <= size(M,2))
        if(M(1,i) == c.LevelList(end-1))
          whiteMatter{k} =   M(:,i+1:i+M(2,i));
          k = k + 1;  
        else
           
            
        end
        i = i + M(2,i)+1; 
    end
    
    %% white matter in out
    sizeMax = 0;
    maxID = [];
    for i = 1:size(whiteMatter,2)
        if(size(whiteMatter{i},2) > sizeMax)
            sizeMax = size(whiteMatter{i},2);
            maxID = i;
        end
        
    end
    whiteMatter
    %% Figure out which contours are inside the white matter
    yM = zeros(1,size(whiteMatter,2));
   
    for i = 1:size(whiteMatter,2)
        for j = 1:size(whiteMatter,2)
           if(~isempty(whiteMatter{i})&& ~isempty(whiteMatter{j})) 
               if(i ~= j)
                    yM(i) = completeInside(whiteMatter{i},whiteMatter{j});
               end
           end
        end
    end
    
    %% Find if a cluster is in white matter
    
    whitePerc = 0;
    for j = 1:size(whiteMatter,2)
       if(yM(j) ~= 1)
           pInside = inpoly2(imagePixels,whiteMatter{j}');
           whitePerc = whitePerc + sum(pInside);
       end
    end
    
    %whitePerc = 100*whitePerc/size(imagePixels,1);
    
    for i = 1:size(Clusters,2)
        y(i) = 0;
        %for k = 1:size(Clusters{i},1)
        if(~isempty(whiteMatter) && ~isempty(maxID))
            for j = 1:size(whiteMatter,2)
                if(yM(j) ~= 1)
                    if(~isempty(Clusters{i}) && ~isempty(whiteMatter{j}))
                        isInside = inpoly2(Clusters{i},whiteMatter{j}');
                        
                        if(100*(sum(isInside)/size(Clusters{i},1)) > perc)
                           for k = 1:size(yM,2)
                              if(yM(k) == 1)
                                 isInsideG = inpoly2(Clusters{i},whiteMatter{k}');
                                 clusOverStat(i) = 100*sum(isInsideG)/size(Clusters{i},1);
                                 if(100*sum(isInsideG)/size(Clusters{i},1) > perc  )
                                    y(i) = 0; 
                                    break;

                                 else
                                    y(i) = 1; 

                                 end

                              end
                           end

                        end
                        
                       
                        
                    end
                end
            end
            
        end
        if(~isempty(Clusters{i}))
            %inPolyI = inpoly2(whiteMatter{maxID}',Clusters{i});
            inPolyI2 = inpoly2(Clusters{i},whiteMatter{maxID}');
            clusOverStat(i) = 100*sum(inPolyI2)/size(Clusters{i},1);
            if( sum(inPolyI2) > (perc/100)*size(Clusters{i},1)) %% sum(inPolyI) > (perc/100)*size(whiteMatter{maxID},2) ||
               y(i) = 1; 
            end
        end
    end
    
else
    
    
end
