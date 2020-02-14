function generateLocationMap(options,kalmanOptions)
% function generateLocationMap(options,kalmanOptions)
%
% Generate measurement location map from well and seismic location
% registries. Saves the location map to `LMap.mat`.
%
% Copyright (c) 2016-2017 IRIS, All Rights Reserved.
% $Id: //depot/rfmatlab/anma/mlabScripts/generateLocationMap.m#1 $
% $DateTime: 2017/12/19 16:10:07 $


LMap = cell(options.dim);

measIndex = 1;
if isfield(kalmanOptions,'useProd') && kalmanOptions.useProd
    load('wellLocation', 'wellLocation');

    for K=1:length(wellLocation)
        wellLocationStep = wellLocation{K};
        for I=1:size(wellLocationStep, 1)
            wl = wellLocationStep(I, :);
            LMap{wl(1), wl(2), wl(3)}(end+1) = measIndex;
            measIndex = measIndex + 1;
        end
    end
end

if getOption(kalmanOptions,'useSeismic',0) || getOption(kalmanOptions,'useMRI',0)
    if getOption(kalmanOptions,'useSeismic',0) ~= 0
        load('inputData', 'seismicOptions');
        nStep = length(seismicOptions.seismic_time);
    elseif getOption(kalmanOptions,'useMRI',0) ~= 0 && ...
           isfield(kalmanOptions,'reportTime')
        T = kalmanOptions.reportTime;
        if getOption(kalmanOptions,'thinobs',0) > 0
            T = T(kalmanOptions.measInd); 
        end
        nStep = length(T);
    else
        nStep = 1;
    end
    load('seisLocation', 'seisLocation');

    for K=1:nStep
        seisLocationStep = seisLocation{K}; %#ok<*USENS>
        for I=1:length(seisLocationStep)
            minCoord = seisLocationStep(I).minCoord;
            maxCoord = seisLocationStep(I).maxCoord;

            dataShape = maxCoord - minCoord + [1 1 1];
            dataCount = prod(dataShape);

            c = zeros(1, 3);
            for D=1:dataCount
                [ c(1), c(2), c(3) ] = ind2sub(dataShape, D);
                c = c + minCoord - [1 1 1];

                LMap{c(1), c(2), c(3)}(end+1) = measIndex;
                measIndex = measIndex + 1;
            end
        end
    end
end

if ~hasSize(LMap, options.dim)
    error('LMap changed size to [%d, %d, %d]', size(LMap))
end

save('locationMap', 'LMap');

end
