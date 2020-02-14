function []=localMDA(iter,ensemble,simData,W,Wbase,H,measurement,scale,options,kalmanOptions,dir)
% function []=localMDA(iter,ensemble,simData,W,Wbase,H,measurement,scale,options,kalmanOptions,dir)
%
% ES-MDA filter for local analysis.
%
% Copyright (c) 2016-2017 IRIS, All Rights Reserved.
% $Id: //depot/rfmatlab/anma/mlabScripts/localMDA.m#1 $
% $DateTime: 2017/12/19 16:10:07 $

% This is assumed all over the place
if any(abs(W - ones(size(W, 1), 1)) > 1e-8)
    error('W must be identity')
end

% Wbase is variance, sWbase is std
sWbase = sqrtW(Wbase);

% Truncate debug files only if iter == 0
debugFilesMode = 'a';
if iter == 0
    debugFilesMode = 'w';
end

fidLocalMDA=fopen('debugLocalMDA', debugFilesMode);
tline = 'iteration number ';
fprintf(fidLocalMDA,'%s %d\n',tline,iter);

if (~isfield(kalmanOptions, 'useSeismic') || ~kalmanOptions.useSeismic) &&...
        size(H,1)~=size(measurement,1)
    error('localMDA.m: the dimension of H and measurement is not consistent!')
end

if size(W,1)~=size(measurement,1)
    error('localMDA.m: the dimension of W and measurement is not consistent!')
end

if isfield(kalmanOptions,'useSeismic') && kalmanOptions.useSeismic &&...
        kalmanOptions.numParProcesses == 1
   error('localMDA.m: Seismic only suppoted in parallel runs')
end

ne=size(ensemble,2); % ensemble size
nd=size(W,1); % observation size
numFreeParam = sum(options.freeparam);
numStaticVar = size(options.staticVar, 1);

% generate perturbed data before the start of iteration
perturbedData = repmat(measurement, 1, ne) + randn(length(measurement), ne);

if ~isempty(scale)
    simData=simData.*repmat(scale,1,ne+1);
end
simData=normalizeData(simData, sWbase);

[obj,objStd,objReal]=getDataMismatch(simData,W,measurement);
disp(['obj=',num2str(obj),' objStd=',num2str(objStd)]);
fprintf(fidLocalMDA,'%s \n',' ');
fprintf(fidLocalMDA,'%s %d \t %f \t  %f \n','LocalMDA ', ...
    iter,obj,objStd);
save(strcat(dir,'objRealIter',num2str(iter)),'objReal');

localInfo = getLocalInfo(options, kalmanOptions);

% start the loop for iteration
startIter = iter+1;
for iter=startIter:kalmanOptions.numIter
    % Calculate anomalies
    deltaM = ensemble - repmat(mean(ensemble, 2), 1, ne);
    deltaD = simData - repmat(mean(simData, 2), 1, ne);
    X = perturbedData - simData;
    inc = zeros(size(deltaM));

    if isfield(kalmanOptions, 'localAnalysis') && kalmanOptions.localAnalysis
        % Update free parameters (no location available)
        inc(1:numFreeParam, :) = getIncrement(iter, deltaM(1:numFreeParam, :), deltaD, X, kalmanOptions);

        lastPct = -1;
        for Icell=1:options.numGridBlocks
            progPct = round(Icell / options.numGridBlocks * 100);
            if progPct ~= lastPct
                if lastPct > -1
                    fprintf(repmat('\b', 1, 25));
                end
                fprintf('Performing analysis  %3d%%', progPct);
                lastPct = progPct;
            end

            % Get local data
            areaInfo = getLocalArea(Icell, localInfo);

            deltaDl = getLocalMeasurements(deltaD, areaInfo);
            if size(deltaDl, 1) == 0
                % No data for this cell, skip update
                continue;
            end

            Xl = getLocalMeasurements(X, areaInfo);

            % Get local params
            localParam = 1:options.numGridBlocks:(options.numGridBlocks * numStaticVar);
            localParam = localParam + numFreeParam + Icell - 1;
            deltaMl = deltaM(localParam, :);

            inc(localParam, :) = getIncrement(iter, deltaMl, deltaDl, Xl, kalmanOptions);
        end
        fprintf('\n');
    else
        inc = getIncrement(iter, deltaM, deltaD, X, kalmanOptions);
    end

    % Apply increment
    ensemble = ensemble + inc;

    if isfield(options,'postProcessEnsemble') && ...
            ~strcmp(options.postProcessEnsemble,'')
        [ensemble,~] = feval(options.postProcessEnsemble,ensemble);
    end

    % put the ensemble members within the bounds:
    if ~isfield('options',' writePerm') && (isfield(kalmanOptions,'staticVarLB') || isfield(kalmanOptions,'staticVarUB')) ...
            && (~isempty(kalmanOptions.staticVarLB) || ~isempty(kalmanOptions.staticVarUB))
        [ensemble,ncut]=respectBounds(ensemble,kalmanOptions,options);
        disp(['average number of truncate per realization is ',num2str(ncut)]);
        fprintf(fidLocalMDA,'%s %f \n','average number of truncate per realization is ',ncut);
    end

    % put free parameters within bounds:
    if (isfield(kalmanOptions,'freeparamLB') || isfield(kalmanOptions,'freeparamUB')) ...
            && (~isempty(kalmanOptions.freeparamLB) || ~isempty(kalmanOptions.freeparamUB))
        [ensemble,ncut]=respectBoundsWithReport(ensemble,kalmanOptions,options);
        disp(['average number of truncate per realization is ',num2str(ncut)]);
        fprintf(fidLocalMDA,'%s %f \n','average number of truncate per realization is ',ncut);
    end

    filecontents=[];
    options.existTimeSteps=0;
    
    % simulatedEnsemble will only have predicted data required in SUMMARY section
    numberOfFailedSimulations = [];
    if kalmanOptions.numParProcesses<2
        [simulatedEnsemble,~,numberOfFailedSimulations]= ...
            runForwardSimulations(ensemble,0,kalmanOptions,options,[], ...
            filecontents);
        jtmp=size(simulatedEnsemble,1)*size(simulatedEnsemble,2); % number of total data
        itmp=size(simulatedEnsemble,3); % number of realizations
        simulatedEnsemble=reshape(simulatedEnsemble,jtmp,itmp);
        simData=simulatedEnsemble(H(:,2),:);
    else
        [simData,~]= runMultiCoreSim(0,ensemble,H,kalmanOptions,options);
    end
    
    save(strcat(dir,'simulatedDataIter',num2str(iter)), ...
        'simData','numberOfFailedSimulations')
    clear simulatedEnsemble;
    
    if ~isempty(numberOfFailedSimulations)
        error('Stopped: Some simulations have failed!')
    end
    
    if ~isempty(scale)
     simData=simData.*repmat(scale,1,ne+1);
    end
    simData=normalizeData(simData, sWbase);
    
    [objNew,objStdNew,objReal]=getDataMismatch(simData,W,measurement); %#ok<*ASGLU>
    disp(['objNew=',num2str(objNew),' objStdNew=',num2str(objStdNew)]);
    
    %--
    fprintf(fidLocalMDA,'%s \n',' ');
    fprintf(fidLocalMDA,'%s %d \t %f \t %f \t  %f \n','LocalMDA ', ...
        iter,kalmanOptions.lambda(iter),objNew,objStdNew);
    save(strcat(dir,'objRealIter',num2str(iter)),'objReal');
    
    % save files
    save(strcat(dir,'ensemble',num2str(iter)),'ensemble');
end

fclose(fidLocalMDA);

end

function increment = getIncrement(iter, deltaM, deltaD, X, kalmanOptions)
    ne = size(deltaM, 2);

    [Ud,Wd,Vd]=svd(deltaD,'econ');
    Wd=diag(Wd); % a vector
    fraction = cumsum(Wd)/sum(Wd);
    svdPd = find(fraction > kalmanOptions.tsvdData, 1);
    Vd=Vd(:,1:svdPd);
    Ud=Ud(:,1:svdPd);
    Wd=Wd(1:svdPd);

    Di = Wd./((ne - 1) * kalmanOptions.lambda(iter) + Wd.^2);
    increment = deltaM*(Vd*diag(Di)*(Ud'*X));
end
