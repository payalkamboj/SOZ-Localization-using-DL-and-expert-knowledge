function y = clusterWidthComputation(newI, newIB, Clusters)
y = [];
modI = zeros(size(newI,1),size(newI,2));
for i = 1:size(Clusters,2)
    
     if(~isempty(Clusters{i}))       
         modI(Clusters{i}(:,1),Clusters{i}(:,2)) = 1;


        xWidth = 0;
        for j = 1:size(modI,1)
            [G,H] = find(modI(j,:) > 0);
            if(~isempty(H))
                xWidthTemp = max(H)-min(H);
                if(xWidthTemp > xWidth)
                    xWidth = xWidthTemp;
                end
            end
        end
        yWidth = 0;
        for j = 1:size(modI,2)
            [G,H] = find(modI(:,j) > 0);
            if(~isempty(H))
                yWidthTemp = max(G)-min(G);
                if(yWidthTemp > yWidth)
                    yWidth = yWidthTemp;
                end
            end
        end

        y = [y; [xWidth yWidth]];
     else
         
        y = [y; [0 0]]; 
     end
end

