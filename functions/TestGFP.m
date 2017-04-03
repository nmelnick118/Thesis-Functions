function [cellsOfInterest, stalkedIdx, stalked] = TestGFP(phase)

% I = imread('BB130_LPho_001.nd2 - s=1 - c=2 - z=0 - t=0.tif');
cd ~/Desktop/'Final Data Sets'/TimeCourse/801/
% read in files
contourFiles = dir('*pill_MESH.mat');



for i = 1:length(contourFiles)
    contourImage = contourFiles(i).name;
    C = strsplit(contourImage, '-');
    basename = char(C(1));
    fluorImage = sprintf('%s-gfp.tif', basename);
    image = imread(fluorImage); %gfp image
    image = uint8(image);
    M = image;
numcols = size(M, 2);

% pick an intensity threshold for your GFP signal
threshold = 125;

pixels = find(M > threshold);

% initial elimination of nearby fluorescence to reduce data
newPixels = pixels(1);  % assign the first one automatically
for i = 2:length(pixels)
    if abs(pixels(i) - pixels(i-1)) > 1
        newPixels = [newPixels pixels(i)];
    end
end


%-------------------------------------------------------------------------%
%{
% Generate Test Image
rgbImage = cat(3, M, M, M);

red = rgbImage(:,:,1);
green = rgbImage(:,:,2);
blue = rgbImage(:,:,3);

% [green_image, map] = imread('green.png') for color
red(newPixels) = 148;
green(newPixels) = 253;
blue(newPixels) = 48;

newRGB = cat(3, red, green, blue);
imshow(newRGB);

%}

% matrix1 = [1:4; 5:8; 9:12];
numrows = size(M, 1);
numcols = size(M, 2);
% IND = [newPixels];
s = [numrows,numcols];      % dimensions of the image
[I,J] = ind2sub(s,newPixels); % convert pixel indices to coordinates

%-------------------------------------------------------------------------%
% Pick cells that have a signal
cells = load(contourImage);
numCells = length(cells.frame.object);
fluorcells = zeros(1, numCells);
fluorI = zeros(1, numCells);
fluorJ = zeros(1, numCells);
for iCell = 1:numCells
    Xperim = cells.frame.object(iCell).Xperim;
    Yperim = cells.frame.object(iCell).Yperim;
    
    in = inpolygon(J,I,Xperim,Yperim);
    all_I = I(find(in == 1));
    all_J = J(find(in == 1));
    if (max(all_I) - min(all_I)) > 5 | (max(all_J) - min(all_J) > 5)
        count = count + 1;
        continue
    end
    
    fluorcells(iCell) = max(in);
    if (max(in) == 1)
        fluorI(iCell) = I(find(in, 1, 'first'));
        fluorJ(iCell) = J(find(in, 1, 'first'));
    else
        fluorI(iCell) = 0;
        fluorJ(iCell) = 0;
    end
end
end

%display(numCells);
%display(length(fluorcells));
%display(sum(fluorcells));

% cells with fluorescence to be examined
cellsOfInterest = cells.frame.object(find(fluorcells));
fluorI = fluorI(fluorI > 0);
fluorJ = fluorJ(fluorJ > 0);

%-------------------------------------------------------------------------%
numCellsOfInterest = length(cellsOfInterest);
pole1distToFluor = zeros(1, numCellsOfInterest);
pole2distToFluor = zeros(1, numCellsOfInterest);
stalked = zeros(1, numCellsOfInterest);
stalkedIdx = zeros(1, numCellsOfInterest);

for iCell = 1:numCellsOfInterest

% identify the pole of the cells (2 points on contour closest to centerline)
Xcont = cellsOfInterest(iCell).Xcont;
Ycont = cellsOfInterest(iCell).Ycont;
centerline = cellsOfInterest(iCell).centerline;
if (isEmpty(centerline))
    continue;
end
contourDist1 = sqrt((Xcont-centerline(1,1)).^2+(Ycont-centerline(1,2)).^2);
[dist1, pole1Idx] = min(contourDist1);

contourDist2 = sqrt((Xcont-centerline(end,1)).^2+(Ycont-centerline(end,2)).^2);
[dist2, pole2Idx] = min(contourDist2);

% which pole is closer to fluorescence?
pole1distToFluor(iCell) = sqrt((Xcont(pole1Idx) - fluorJ(iCell))^2 + (Ycont(pole1Idx) - fluorI(iCell))^2);
pole2distToFluor(iCell) = sqrt((Xcont(pole2Idx) - fluorJ(iCell))^2 + (Ycont(pole2Idx) - fluorI(iCell))^2);

% assign stalked pole based on distance to fluorescence
if pole1distToFluor(iCell) <= pole2distToFluor(iCell)
    stalkedIdx(iCell) = pole1Idx;
    stalked(iCell) = 1;
else
    stalkedIdx(iCell) = pole2Idx;
    stalked(iCell) = 2;
end

end


% return cellsOfInterest and (stalkedIdx and/or stalked)

