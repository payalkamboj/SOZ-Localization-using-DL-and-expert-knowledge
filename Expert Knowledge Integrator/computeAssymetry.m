function y = computeAssymetry(VR)
if(isempty(VR))
    y = 1000;
    return;
end
midPoint = floor(size(VR,2)/2)+1;
if(sum(sum(VR-VR(1,1))) == 0)
    y = 1000;
else
    y = 0;
    for i = 1:size(VR,1)
       assymDiff = abs(VR(i,1:midPoint-1)-VR(i,midPoint+1:end)); 
       y = y + sum(assymDiff); 
    end
end