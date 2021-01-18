function [simulatedEnsemble,filecontentsOut,numberOfFailedSimulations]= ...
    runForwardSimulations(ensemble,~,kalmanOptions,options,~,~)     

% Run the forward simulations.


% loop over ensemble
disp('Running realization: ');
for I = 1:size(ensemble,2)
    
    xValue = ensemble(:,I);
    options = setOptions(options,xValue);
    
    [state]=runFlowSolver3D2QSort(options);
    
    for J = 0:log10(I-1)
       fprintf('\b'); % delete previous counter display
    end
    fprintf('%d', I);
    
    CA = interp1(options.time',state.volTracerArt',kalmanOptions.reportTime);
    CA = options.porosityArt(:).*CA';
    CQ = interp1(options.time',state.volTracerCap',kalmanOptions.reportTime);
    CQ = options.porosityQ(:).*CQ';
    CV = interp1(options.time',state.volTracerVen',kalmanOptions.reportTime);
    CV = options.porosityVen(:).*CV';
    concentration = CA + CQ + CV;
    
    if strcmp(getOption(kalmanOptions,'obsType','concentration'),'concentration')
        measurement = concentration;
        if getOption(kalmanOptions,'thinobs',0) > 0
            % thin out to reduce amount of data
            measurement = measurement(:,kalmanOptions.measInd); 
        end
    else
        concentration = reshape(concentration,options.nx,options.ny,...
                                              options.nz,options.timeEnd);
        [CcrsMax, argMax] = max(concentration,[],4);
        delay = 0;
        tofArt = (argMax-delay) .* options.mask(:,:,:,1);
        if strcmp(kalmanOptions.obsType,'tof')
            measurement = tofArt(:);
        elseif strcmp(kalmanOptions.obsType,'tofAndPeak')
            measurement = [CcrsMax(:),tofArt(:)];
        end
    end
        
    if ~strcmp(kalmanOptions.ES_script,'localMDA')
        % only keep measurements in active cells
        measurement = measurement(options.actnum==1,:);
    end
    
    simulatedEnsemble(:,1,I) = measurement(:);  %#ok<*AGROW>
    filecontentsOut.rateQ(:,I)=state.rateQ(:);
    filecontentsOut.concentration(:,I)=concentration(:);
end

fprintf('\n');
%filecontentsOut = [];
numberOfFailedSimulations = [];
