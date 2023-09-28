function y = removeR(image,Offset,endOff)
y = image;
for i = 1:Offset
    for j = 1:Offset
        y(i,j,:) = [0 0 0];
    end
end
for i = size(y,1)-endOff:size(y,1)
    for j = 1:Offset
        y(i,j,:) = [0 0 0];
    end
end