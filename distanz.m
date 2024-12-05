function d = distanz(XCoords,YCoords)
% Euclidean distance between each pair of points
% E.g: If 3 points => X (x1,x2,x3) and Y (y1,y2,y3) => d = 3x3 matrix
    N = numel(XCoords);
    matrix_X = repmat(XCoords,1,N);
    matrix_Y = repmat(YCoords,1,N);
    d = sqrt((matrix_X-matrix_X').^2 + (matrix_Y-matrix_Y').^2);
end
