function localInfo = getLocalInfo(options, kalmanOptions)
% function localInfo = getLocalInfo(options, kalmanOptions)
%
% Calculate local info.
%
% Helper functions for local analysis need some common data structures. This
% prevents those structures from being calculated several times by instead
% storing them in a struct.
%
% Outputs
% localInfo
%           .dim: Dimension of the grid.
%           .activeToCoord: Map for converting active index to coordinate.
%           .coordToActive: Map for converting coordinate to active index.
%           .radius: Requested radius of local area.
%           .LMap: Map of measurement indices in each cell.
%           .taperFunc: Function to convert distance to weight.
%
% Copyright (c) 2016-2017 IRIS, All Rights Reserved.
% $Id: //depot/rfmatlab/anma/mlabScripts/getLocalInfo.m#1 $
% $DateTime: 2017/12/19 16:10:07 $

actnum = logical(options.actnum);
actnum = reshape(actnum, options.dim);

localInfo.dim = options.dim;
localInfo.activeToCoord = getActiveToCoord(actnum);
localInfo.coordToActive = getCoordToActive(actnum);
localInfo.radius = kalmanOptions.localRadius;

if isfield(kalmanOptions, 'LMap')
    localInfo.LMap = kalmanOptions.LMap;
elseif isfield(kalmanOptions, 'LMapFile')
    lmf = load(kalmanOptions.LMapFile);
    localInfo.LMap = lmf.LMap;
end

if ischar(kalmanOptions.localTaperFunc)
    if strcmp(kalmanOptions.localTaperFunc, 'linear')
        localInfo.taperFunc = @(r) (1 - r) .* (r < 1);
    elseif strcmp(kalmanOptions.localTaperFunc, 'sphere')
        localInfo.taperFunc = @(r) sqrt((1 - r.^2) .* (r < 1));
    elseif strcmp(kalmanOptions.localTaperFunc, 'parzen')
        localInfo.taperFunc = @(r) ...
            (1 - 6*r.^2.*(1-r)) .* (r < 0.5) +...
            2*(1 - r).^3 .* (r >= 0.5) .* (r < 1);
    elseif strcmp(kalmanOptions.localTaperFunc, 'sine')
        localInfo.taperFunc = @(r) cos(pi/2 * r) .* (r < 1);
    else
        error('Unknown taper function ''%s''', kalmanOptions.localTaperFunc);
    end
else
    localInfo.taperFunc = kalmanOptions.localTaperFunc;
end
end

function activeToCoord = getActiveToCoord(actnum)
% function activeToCoord = getActiveToCoord(actnum)
%
% Inputs
% actnum:           A 3-dimensional array of active cells.
%
% Outputs
% activeToCoord:    A `N x 3` dimensional array where `N` is the number of
%                   active cells. `activeToCoord(i, :)` is the 3-dimensional
%                   coordinate of active cell index `i`.

indexMap = cell(1, 3);
linIndex = find(actnum);
[ indexMap{:} ] = ind2sub(size(actnum), linIndex);
activeToCoord = cell2mat(indexMap);
end

function coordToActive = getCoordToActive(actnum)
% function coordToActive = getCoordToActive(actnum)
%
% Inputs
% actnum:           A `d`-dimensional array of active cells.
%
% Outputs
% coordToActive:    A `d`-dimensional array containing the active index of
%                   each cell, or 0 if the cell is inactive.

coordToActive = reshape(cumsum(actnum(:)) .* actnum(:), size(actnum));
end
