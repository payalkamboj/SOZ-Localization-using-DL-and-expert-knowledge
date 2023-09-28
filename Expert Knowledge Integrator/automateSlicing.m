
function [startX, startY, sizeX, sizeY,numRow,numCol]= automateSlicing(mainImage,templateImage)
    
    %mainImage    = imread('IC_2_thresh.png');
    %templateImage  = imread('R.png');
    %xmin ymin width height.
    %figure
    %imshow(imcrop(mainImage,[1 134 84 114])); correct scalings
    
    % Convert images to grayscale if colored
    if size(mainImage, 3) == 3
        mainImage = rgb2gray(mainImage);
    end
    if size(templateImage, 3) == 3
        templateImage = rgb2gray(templateImage);
    end
    
    % Normalized cross-correlation to find template matches of "R"
    correlationResult = normxcorr2(templateImage, mainImage);
    
    % Threshold for correlation values to be considered as matches
    threshold = 0.7; 
    
    % Coordinates of matches above the threshold
    [rows, cols] = find(correlationResult > threshold);
    
    % Number of rows and columns in the IC
    unique_rows = unique(rows);
    unique_cols = unique(cols);
    
    numRow = numel(unique_rows)-1; % -1 because we don't need top most row of IC
    numCol = numel(unique_cols);
    
    % IC image with bounding boxes around all detected "R"s
    %figure;
    %imshow(mainImage);
    %hold on;
    %for i = 1:numel(rows)
    %    yOffset = rows(i) - size(templateImage, 1);
    %    xOffset = cols(i) - size(templateImage, 2);
    %    rectangle('Position', [xOffset, yOffset, size(templateImage, 2), size(templateImage, 1)], 'EdgeColor', 'r', 'LineWidth', 2);
    %end
    %title('Detected "R"s');
    
    % Output the coordinates of all detected "R"s
    %The detected R coordinates (row and column here) represents top left corner of
    %the template image, so the real R coordinates are little shifted
    %fprintf('Detected "R"s are located at the following coordinates:\n');
    %for i = 1:numel(rows)
    %    y = rows(i) - size(templateImage, 1)/3; 
    %    x = cols(i) - size(templateImage, 2)/3;
    %    fprintf('R%d: (x=%d, y=%d)\n', i, x, y);
    %end
    startY = rows(1) - size(templateImage, 1)/3; 
    startX = cols(1) - size(templateImage, 2)/3;

    % Calculate the width and height between Rs
    if numel(unique_cols) >= 2
        % Sort the unique_cols in ascending order
        sorted_cols = sort(unique_cols);
        
        % Calculate the width between the first two R's
        sizeX = sorted_cols(2) - sorted_cols(1);
    end
    
    if numel(unique_rows) >= 2
        % Sort the unique_rows in ascending order
        sorted_rows = sort(unique_rows);
        
        % Calculate the height between two vertically aligned R's
        sizeY = sorted_rows(2) - sorted_rows(1);
        
    end
end