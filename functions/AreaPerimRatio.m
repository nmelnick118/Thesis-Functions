function max_min = AreaPerimRatio(cells)
numCells = length(cells.frame.object);
pole1value = zeros(1, numCells);
pole2value = zeros(1, numCells);
for iCell = 1:numCells
% identify the pole of the cells
Xcont = cells.frame.object(iCell).Xcont;
Ycont = cells.frame.object(iCell).Ycont;
centerline = cells.frame.object(iCell).centerline;
if (isempty(centerline))
    continue;
end
contourDist1 = sqrt((Xcont-centerline(1,1)).^2+(Ycont-centerline(1,2)).^2);
[dist1, pole1Idx] = min(contourDist1);

contourDist2 = sqrt((Xcont-centerline(end,1)).^2+(Ycont-centerline(end,2)).^2);
[dist2, pole2Idx] = min(contourDist2);

% count 3 spots away to each side of the poles
L = length(Xcont);

n = 2;

pole1Idl = pole1Idx - n;
if (pole1Idl < 1) pole1Idl = pole1Idl + L; end
pole2Idl = pole2Idx - n;
if (pole2Idl < 1) pole2Idl = pole2Idl + L; end

% pole1
loopFix = 6 - (L - pole1Idl);
if loopFix > 0
    X = Xcont([1:loopFix, pole1Idl:end]);
    Y = Ycont([1:loopFix, pole1Idl:end]);
else
    X = Xcont(pole1Idl:pole1Idl + 6);
    Y = Ycont(pole1Idl:pole1Idl + 6);
end

perim = polyperim(X, Y);
area = polyarea(X, Y);
% final output pole1
pole1value(iCell) = area/perim;

% pole2
loopFix = 6 - (L - pole2Idl);
if loopFix > 0
    X = Xcont([1:loopFix, pole2Idl:end]);
    Y = Ycont([1:loopFix, pole2Idl:end]);
else
    X = Xcont(pole2Idl:pole2Idl + 6);
    Y = Ycont(pole2Idl:pole2Idl + 6);
end

perim = polyperim(X, Y);
area = polyarea(X, Y);
% final output pole2
pole2value(iCell) = area/perim;
end

maximum = max(pole1value, pole2value);
minimum = min(pole1value, pole2value);

max_min = [minimum; maximum];
end