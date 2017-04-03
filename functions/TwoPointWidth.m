function max_min = TwoPointWidth(cells)
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
% pole1Idl = mod(pole1Idx + L - 3, L) + 1;
% pole1Idr = mod(pole1Idx + 3, L);
% pole2Idl = mod(pole2Idx + L - 3, L) + 1;
% pole2Idr = mod(pole2Idx + 3, L);

n = 2;

pole1Idl = pole1Idx - n;
if (pole1Idl < 1) pole1Idl = pole1Idl + L; end
pole1Idr = pole1Idx + n;
if (pole1Idr > L) pole1Idr = pole1Idr - L; end
pole2Idl = pole2Idx - n;
if (pole2Idl < 1) pole2Idl = pole2Idl + L; end
pole2Idr = pole2Idx + n;
if (pole2Idr > L) pole2Idr = pole2Idr - L; end

% final outputs
pole1value(iCell) = sqrt((Xcont(pole1Idl) - Xcont(pole1Idr)).^2 + (Ycont(pole1Idl)-Ycont(pole1Idr)).^2);
pole2value(iCell) = sqrt((Xcont(pole2Idl) - Xcont(pole2Idr)).^2 + (Ycont(pole2Idl)-Ycont(pole2Idr)).^2);

end

maximum = max(pole1value, pole2value);
minimum = min(pole1value, pole2value);

max_min = [minimum; maximum];
end