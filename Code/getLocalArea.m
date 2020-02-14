function areaInfo = getLocalArea(cellIndex, localInfo)
% function areaInfo = getLocalArea(cellIndex, localInfo)
%
% Compute the area around the given active cell index.
%
% Inputs
% cellIndex:        The active index of the cell to compute area around.
% localInfo:        Struct containing structural information bout the domain.
%           .dim:   Dimension of the domain.
%           .activeToCoord: Active cell to coordinate map.
%           .radius: The radius of the local domain.
%
% Outputs
% areaInfo
%         .localInfo: Local info used to constrct this area.
%         .minArea: Minimum coordinate of the local domain for each axis.
%         .minArea: Maximum coordinate of the local domain for each axis.
%         .distances: The distance to the current cell for each cell in the local
%                     domain. Has the shape of the local domain.
%
% Copyright (c) 2016-2017 IRIS, All Rights Reserved.
% $Id: //depot/rfmatlab/anma/mlabScripts/getLocalArea.m#1 $
% $DateTime: 2017/12/19 16:10:07 $

ndim = size(localInfo.activeToCoord, 2);
coord = localInfo.activeToCoord(cellIndex, :);

radius = localInfo.radius;
s = max(1, coord - radius);
e = min(localInfo.dim, coord + radius);

coord_slices = cell(ndim, 1);
for I=1:ndim
    range = s(I):e(I);
    coord_slices{I} = range(:);
end

grid = cell(ndim, 1);
[ grid{:} ] = ndgrid(coord_slices{:});

distances = zeros(e - s + ones(size(s)));
for I=1:ndim
    distances = distances + (grid{I} - coord(I)).^2;
end

areaInfo.localInfo = localInfo;
areaInfo.minArea = s;
areaInfo.maxArea = e;
areaInfo.distances = sqrt(distances);

end
