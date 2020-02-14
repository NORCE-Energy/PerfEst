function W = adjustMeasurementNoise

% function W = adjustMeasurementNoise
%
% Adjust measurement noise based on distance to simulated ensemble.
%
% Copyright (c) 2010-2016 IRIS, All Rights Reserved.
% $Id: //depot/rfmatlab/main/preProcess/adjustMeasurementNoise.m#1 $
% $DateTime: 2018/12/05 14:14:11 $

load trueSolutionSmoother.mat;
load inputData.mat;
load simulatedDataIter0.mat;

N = getOption(kalmanOptions,'ensembleSize',size(simData,2)-1);
for I = 1:length(measurement)
   
    if measurement(I) < min(simData(I,1:N))
        D = abs(measurement(I) - min(simData(I,1:N)));
        W(I) = max(W(I),D^2); %#ok<*AGROW>
    elseif measurement(I) > max(simData(I,1:N))           
        D = abs(measurement(I) - max(simData(I,1:N)));
        W(I) = max(W(I),D^2);
    end
    
end