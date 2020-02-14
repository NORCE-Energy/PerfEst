function localMeas = getLocalMeasurements(meas, areaInfo)
% function localMeas = getLocalMeasurements(meas, areaInfo)
%
% Given a global measurement array and a local area, returns local measurements
% scaled according to tapering function.
%
% Copyright (c) 2016-2017 IRIS, All Rights Reserved.
% $Id: //depot/rfmatlab/anma/mlabScripts/getLocalMeasurements.m#1 $
% $DateTime: 2017/12/19 16:10:07 $

    ndim = length(areaInfo.localInfo.dim);
    weights = areaInfo.localInfo.taperFunc(areaInfo.distances / areaInfo.localInfo.radius);

    coord_slices = cell(ndim, 1);
    for I=1:ndim
        range = areaInfo.minArea(I):areaInfo.maxArea(I);
        coord_slices{I} = range(:);
    end

    localMeasI = areaInfo.localInfo.LMap(coord_slices{:});

    % Look up measurement values
    localMeas = cellfun(@(indices) meas(indices, :), localMeasI, 'UniformOutput', false);

    % Apply distance weights
    localMeas = cellfun(@times, localMeas, num2cell(weights), 'UniformOutput', false);
    localMeas = vertcat(localMeas{:});
end
