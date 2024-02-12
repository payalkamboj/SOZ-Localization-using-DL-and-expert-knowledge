function y = arteryNoiseDetect(Clusters,image,epsilon, perc, imContour)
y = [];
if(~isempty(Clusters))
    [M,c] = contour(rgb2gray(imContour));
%    perLevel = M(1,1);
%    contourData = M;
%     index = 1;
%     k = 1;
%     while(index <= size(M,2))
%         if(M(1,index) == perLevel)
%             outerPot{k} = M(:,index+1:M(2,index)+index);
%             k = k + 1;
%         end
%         index = index + M(2,index)+1;
%     end
%     maxSize = 0;
%     for i = 1:k-1
%        if(maxSize < size(outerPot{i},2))
%            maxSize = size(outerPot{i},2);
%            maxID = i;
%        end
%            
%         
%     end
%     
%     OuterContour = outerPot{maxID}; % very sensitive to dataset, hardcoded value please change this
    %% artery region
    midX = floor(size(image,2)/2);
    
    index = 1;
    k = 1;
    while(index <= size(M,2))
        if(M(1,index) == c.LevelList(1))%c.LevelList(end-2))
            outerPot{k} = M(:,index+1:M(2,index)+index);
            k = k + 1;
        end
        index = index + M(2,index)+1;
    end
    %maxSize = 0;
    
    for i = 1:k-1
        yM(i) = 0;
        for j = 1:k-1
            if(i~=j)
                p = inpoly2(outerPot{i}',outerPot{j}');
                if(sum(p) > (perc/100)*size(outerPot{i},2))
                   yM(i) = 1; 
                   break; 
                end
            end
        end
        
    end
    
    artOuter = [];
    artInner = [];
    kO = 1;
    kI = 1;
     for i = 1:k-1
         if(yM(i) == 0)
            artOuter{kO} = outerPot{i};
            kO = kO + 1;
         else
             artInner{kI} = outerPot{i};
             kI = kI + 1;
         end
         
         
     end
    if(~isempty(artOuter)) 
        maxSize = 0;
        for i = 1:size(artOuter,1)
           if(maxSize < size(artOuter{i},2))
               maxSize = size(artOuter{i},2);
               maxID = i;
           end


        end
    end
    
    artContour = artOuter{maxID};
    [GX,HX] = find(artContour(1,:)<=midX);
    [GXU,HXU] = find(artContour(1,:)>midX);

    [GY,GYM] = min(artContour(2,HX));
    [GYU,GYMU] = min(artContour(2,HXU));
    arterialRegion = [];
%     artStart = artContour(:,HX(GYM));
%     artEnd = artContour(:,HXU(GYMU));
%     if(~isempty(artContour) && ~isempty(artStart) && ~isempty(artEnd))
%         newArtCont = [];
%         startI = 0;
%         for i = 1:size(artContour,2)
%             dist1 = sqrt(sum((artContour(:,i) - artStart).^2));
%             dist2 = sqrt(sum((artContour(:,i) - artEnd).^2));
%             if(startI == 0)
%                if(dist1 == 0 || dist2 == 0)
%                    newArtCont = [newArtCont i];
%                    startI = 1;
%                else
% 
%                end
%             else
%                 if(dist1 == 0 || dist2 == 0)
%                     newArtCont = [newArtCont i];
%                     break;
%                 else
%                     newArtCont = [newArtCont i];
%                 end
%             end
%         end
%         C = setdiff(1:size(artContour,2),newArtCont);
%         arterialRegion = artContour(:,C);
%     end
    
    
    [G,H] = find(artContour(1,:) == midX);
    minY = min(artContour(2,H));
    
    xL = midX - 5;
    xM = midX + 5;
    yL = minY;
    yM = minY + 25;
    if(~isempty(minY))
        arterialRegion = [[xL; yL] [xL; yM] [xM; yM] [xM; yL]];
    else
        
        arterialRegion = [];
    end
    
    %% cluster in the region
    for i = 1:size(Clusters,2)
        y(i) = 0;
        y0(i) = 0;
        yI(i) = 0;
        
        %for k = 1:size(Clusters{i},1)
        if(~isempty(Clusters{i}))
            %[G,H] = find(Clusters{i}(:,1) >= xL & Clusters{i}(:,1) <= xM );
            %[G1,H1] = find(Clusters{i}(G,2) >= yL & Clusters{i}(G,2) <= yM);
            for mp = 1:kI-1
                if(~isempty(artInner{mp}))
                   p = inpoly2(Clusters{i},artInner{mp}'); 
                   pI = inpoly2(artInner{mp}',Clusters{i}); 
                   y(i) = y(i) + sum(p);
                   if(sum(pI) > (perc/100)*size(artInner{mp},2))
                        yI(i) = 1000;
                   end
                end
            end
            
            for mp = 1:kO-1
                if(~isempty(artOuter{mp}))
                   p = inpoly2(Clusters{i},artOuter{mp}'); 

                   y0(i) = y0(i) + sum(p); 
                end
            end
        
        %end
            if(y(i) > (perc/100)*size(Clusters{i},1) || y0(i) < ((1-(perc/100))*size(Clusters{i},1)) || yI(i) == 1000)
               y(i) = -1; 

            end
            
             %region inside cluster
%              if(~isempty(arterialRegion))
%                  p = inpoly2(arterialRegion',Clusters{i});
%                  if(100*(sum(p)/size(p,2)) > perc)
%                     y(i) = -1; 
%                  end
%              end
        end
    end
    
   
    
    
else
    
    
end
          