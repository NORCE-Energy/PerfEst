function []=AGSMDA(iter,ensemble,simData,W,Wbase,H,measurement,scale,options,kalmanOptions,dir,obsLocation,obsType,weights)

% Iterative ensemble smoother with minimum-average-cost
% function []=RLM_average_cost(iter,ensemble,simData,W,Wbase,H,measurement,options,kalmanOptions,dir,obsLocation,obsType)
%
% INPUTS
% ------
% iter:                 index of the current iteration step (starting from step 0)
% ensemble:             ensemble of parameters (no dynamical states included in the iES)
% simData:              simulated measurements w.r.t "ensemble" and is
%                       normalized by the std of observation errors
% W:                    weight vector for measurements (by default all equal to 1 with the normalization)
% Wbase:                Variance of observation variables, as either a vector or a matrix
% H:                    observation operator w.r.t all available
%                       measurements at different time instants
% measurement:          real measurements at at different time instants
% scale:                scale measurements
% options:              Structure defining the reservoir, see
%                       defaultOptions.m for more info
% kalmanOptions:        Structure defining options for the EnKF/iES parameters, see
%                       defaultKalmanOptions.m for more info.
% dir:                  work directory
% obsLocation:          observation location (for covariance localization)
% obsType:              observation type (for covariance localization)
%
%
% Copyright (c) 2010-2014 IRIS, All Rights Reserved.
% $Id: //depot/rfmatlab/main/Kalman/RLM_average_cost.m#6 $
% $DateTime: 2015/07/30 10:41:33 $

h=kalmanOptions.AGSMDA.h;
lambda=kalmanOptions.AGSMDA.lambda;

% Wbase is variance, sWbase is std
sWbase = sqrtW(Wbase);

% Truncate debug files only if iter == 0
debugFilesMode = 'a';
if iter == 0
    debugFilesMode = 'w';
end

fidAGSMDA=fopen('debugAGSMDA', debugFilesMode);
tline = 'iteration number ';
fprintf(fidAGSMDA,'%s %d\n',tline,iter);


if nargin<13
    if isfield(kalmanOptions,'distanceLoc')
        error('AGSMDA.m: if using distanced-based localizaiton, need to specify obsLocation and obsType')
    else
        obsLocation=[];
        obsType=[];
    end
end
obsLocationEff = obsLocation;
obsTypeEff = obsType;

% If kalmanOptions.ignoreUniformativeMesurements=1 we get erronous output if
% the size of H and measurment is not consistent. This test is only valid if we
% don't use seismic.
if (~isfield(kalmanOptions, 'useSeismic') || ~kalmanOptions.useSeismic) &&...
        size(H,1)~=size(measurement,1)
    error('AGSMDA.m: the dimension of H and measurement is not consistent!')
end

if size(W,1)~=size(measurement,1)
    error('AGSMDA.m: the dimension of W and measurement is not consistent!')
end

if isfield(kalmanOptions,'useSeismic') && kalmanOptions.useSeismic &&...
        kalmanOptions.numParProcesses == 1
    error('AGSMDA.m: Seismic only suppoted in parallel runs')
end

ne=size(ensemble,2); % ensemble size
nd=size(W,1); % observation size

% generate perturbed data before the start of iteration
perturbedData=zeros(nd,ne);
for j=1:ne
    perturbedData(:,j)=measurement+W.*randn(size(measurement));
end

if ~isempty(scale)
    simData=simData.*repmat(scale,1,ne+1);
end
simData=normalizeData(simData, sWbase);

[obj,objStd,objReal]=getDataMismatch(simData,W,measurement);
disp(['obj=',num2str(obj),' objStd=',num2str(objStd)]);
fprintf(fidAGSMDA,'%s \n',' ');
fprintf(fidAGSMDA,'%s %d \t %f \t  %f \n','AGSMDA ', ...
    iter,obj,objStd);
save(strcat(dir,'objRealIter',num2str(iter)),'objReal');

if isfield(kalmanOptions,'distanceLoc') % distance-based localization
    fidLoc=fopen('debugDisLoc', debugFilesMode);
    tline = 'starting iteration';
    fprintf(fidLoc,'%s\n',tline);
end

% start the loop for iteration
iter = iter + 1;
while iter <= kalmanOptions.numIter
    
    if isfield(kalmanOptions,'distanceLoc') % distance-based localization
        tline = 'checking for data with zero variability from ensemble prediction';
        fprintf(fidLoc,'%s\n',tline);
        fprintf(fidLoc,'%s %d\n','number of obs before screening ',length(measurement));
    end
    
    fprintf(fidAGSMDA,'%s %d\n','number of production data is ',size(simData,1));
    disp(['-- Iteration step: ' int2str(iter) '--']);
    disp(['number of production data is ',num2str(size(simData,1))]);
    
    % get the simulated observation of ensemble mean
   
    % normalized by C_{D}^{-1/2}
    simMean=simData*weights';
    deltaD = simData - simMean * ones(1,ne);
    
    % Check the difference between simMean and the sample mean of the simulated observations.
    % This is not essential to the iES, and can be skipped
    %obsDiff = Wbase .*(simMean - mean(simData,2));  %#ok<NASGU>
    %alt_deltaD = simData - mean(simData,2)*ones(1,ne); %#ok<NASGU>
    %save(strcat(dir,'simulatedDataDiff','-',num2str(iter)), 'obsDiff','simMean','simData','deltaD','alt_deltaD');
    
    Witer = W;

    % sWbaseIter is always matrix, even if sWbase isn't
    if(size(sWbase,2) == 1)
        sWbaseIter = diag(sWbase);
    else
        sWbaseIter = sWbase;
    end

    if isfield(kalmanOptions,'ignoreUninformativeMeasurements') && ...
            kalmanOptions.ignoreUninformativeMeasurements==1
        
        varTestMeasurements=var(deltaD,0,2)';
        keepVariables=find(varTestMeasurements~=0);
        
        % adjust W, H and measurement:
        %H=H(keepVariables,:); % this should work ok also if kalmanOptions.Hones==1
        Witer = W(keepVariables);

        sWbaseIter = sWbaseIter(keepVariables, keepVariables);
        %measurement=measurement(keepVariables);
        
        deltaD=deltaD(keepVariables,:);
        simData=simData(keepVariables,:);
        perturbedDataEff=perturbedData(keepVariables,:);
        
        if isfield(kalmanOptions,'distanceLoc') % distance-based localization
            
            % need to exclude statistics obs. (used in facies studies)
            tmp_options = load('inputData.mat','options');
            if isfield(tmp_options.options,'biModel') && (~isempty(tmp_options.options.biModel))
                % For facies studies, additional channel statistics may be
                % used as observations. So 'real_obs_index' below refer to
                % the indices of real observations (e.g., production data),
                % whereas 'obs_type_index' also includes the indices of artifical
                % observations (e.g., channel statistics).
                load('trueSolutionSmoother','real_obs_index','obs_type_index'); %
                
                keepVariables_realObs = intersect(real_obs_index,keepVariables)';
                tmp_index = find(real_obs_index == keepVariables_realObs);
                kept_obs_type_index = obs_type_index(tmp_index); % exclue real observations with zero variations.
                
                obsTypeEff=obsType(kept_obs_type_index,:);
                obsLocationEff=obsLocation(kept_obs_type_index,:);
            else
                load('trueSolutionSmoother','num_prod');
                keepProdVariables = find(varTestMeasurements(1:num_prod)~=0);
                obsTypeEff=obsType(keepProdVariables,:);
                obsLocationEff=obsLocation(keepProdVariables,:);
            end
            clear tmp_options;
            
            fprintf(fidLoc,'%s %d \n','number of obs after screening ',size(obsTypeEff,1));
        end
        disp(['number of data for updating is ',num2str(size(simData,1))])
        
        fprintf(fidAGSMDA,'%s %s\n','number of data for updating is ',num2str(size(simData,1)));
        
        % after this we continues as before, except we fill in zero
        % columns in Ke before reporting, except if all measurements are
        % removed:
        if sum(keepVariables)==0
            warning('rfmatlab:Kalman:AGSMDA','no variability in ensemble prediction');
            return;
        end
    else
        perturbedDataEff = perturbedData;
        obsTypeEff=obsType;
        obsLocationEff=obsLocation;
    end
    
    % svd of deltaD
    [Ud,Wd,Vd]=svd(deltaD,'econ');
    
    val=diag(Wd);
    total=sum(val);
    for j=1:ne
        svdPd=j;
        if (sum(val(1:j))/total > kalmanOptions.tsvdData)
            break
        end
    end
    disp(['svdPd=',num2str(svdPd)]);
    fprintf(fidAGSMDA,'%s %d\n','number of singular value retained ',svdPd);
    Vd=Vd(:,1:svdPd);
    Ud=Ud(:,1:svdPd);
    Wd=val(1:svdPd); % a vector
    clear val;
    
    x1=zeros(ne,svdPd);
    for jVec=1:svdPd
        for j=1:ne
            x1(j,jVec)=Vd(j,jVec)*Wd(jVec)/(h(iter)^(-2)*lambda(iter)+Wd(jVec)^2);
        end
    end
    
    % do the localization
    X = (simData(:,1:ne)-perturbedDataEff);
    I = localization(ensemble,x1,Ud,X,obsTypeEff,obsLocationEff,weights);
    
    % update
    save('TempEns','ensemble')
    ensemble=ensemble - I;
    save('TempEns2','ensemble')
    % weights
    Sigma=h(iter)^2*(sWbaseIter*Ud)*diag(Wd.^2)*(sWbaseIter*Ud)'+lambda(iter)*diag(Witer.^2);
    delta=(simData-perturbedDataEff);
    
    dist = zeros(1,size(ensemble,2));
    for i=1:size(ensemble,2)
        dist(i)=(delta(:,i)'/Sigma)*delta(:,i);
    end
    dist=dist-min(dist);
    weights=weights.*exp(-0.5*dist);
    weights=weights/sum(weights);
    weights=max(weights,0);
    
    if (kalmanOptions.AGSMDA.adaptive==0)
               
        alpha=options.AGSMDA.alpha;
        weights=alpha*weights+(1-alpha)*ones(1,size(ensemble,2));
        weights=weights./sum(weights);
        % Neff=1/sum(weights.^2);
    
    else
        
        Neff=1/sum(weights.^2);
        alpha=Neff/size(ensemble,2);
        weights=alpha*weights+(1-alpha)*ones(1,size(ensemble,2))/size(ensemble,2);
        
    end
    meanstate=ensemble*weights';
    LPae=zeros(size(ensemble));
    for i=1:size(ensemble,2)
        LPae(:,i)=(h(iter)*sqrt(weights(i)/(1-sum(weights.^2))))*(ensemble(:,i)-meanstate);
    end
    
    % Resampling --NEW: inplementing residual resampling.
    if  kalmanOptions.AGSMDA.resample==1
        if (Neff>kalmanOptions.AGSMDA.Nefflim1*size(ensemble,2))&&(Neff<kalmanOptions.AGSMDA.Nefflim2*size(ensemble,2))
            index=[];
            
            Nbar=floor(size(ensemble,2)*weights);
            R=sum(Nbar);
            Nnew=size(ensemble,2)-R;
            newweights=(size(ensemble,2)*weights-Nbar)./Nnew;
            
            % Hindre numeriske feil slik at summen alltid er en
            newweights=newweights./sum(newweights);
            
            
            for i=1:size(ensemble,2)
                index=[index i*ones(1,Nbar(i))]; %#ok<AGROW>
            end
            if (size(index,2)<size(ensemble,2))
                newindex=randsample(1:size(ensemble,2),Nnew,true,newweights);
                finalindex=[index newindex];
            else
                finalindex=index;
            end
            
            ensemble=ensemble(:,finalindex)+h(iter).*LPae*(randn(size(ensemble,2), ...
                size(ensemble,2)));
            
            
            meanstate=ensemble*weights';
            LPae=zeros(size(ensemble));
            for i=1:size(ensemble,2)
                LPae(:,i)=(h(iter)*sqrt(weights(i)/(1-sum(weights.^2))))*(ensemble(:,i)-meanstate);
            end
            
            
            weights=ones(1,size(ensemble,2));
            weights=weights./sum(weights)+(1/size(ensemble,2))*(1/2^53);
            
            
        end
        disp('Resampling has been performed')
    end
    
    % re-initialize the facies field in cases of facies studies
    if isfield(options,'postProcessEnsemble') && ...
            ~strcmp(options.postProcessEnsemble,'')
        [ensemble,~] = feval(options.postProcessEnsemble,ensemble);
    end
    
    
    
    % put the ensemble members within the bounds:
    if ~isfield('options',' writePerm') && (isfield(kalmanOptions,'staticVarLB') || isfield(kalmanOptions,'staticVarUB')) ...
            && (~isempty(kalmanOptions.staticVarLB) || ~isempty(kalmanOptions.staticVarUB))
        [ensemble,ncut]=respectBounds(ensemble,kalmanOptions,options);
        disp(['average number of truncate per realization is ',num2str(ncut)]);
        fprintf(fidAGSMDA,'%s %f \n','average number of truncate per realization is ',ncut);
    end
    
    % put free parameters within bounds:
    if (isfield(kalmanOptions,'freeparamLB') || isfield(kalmanOptions,'freeparamUB')) ...
            && (~isempty(kalmanOptions.freeparamLB) || ~isempty(kalmanOptions.freeparamUB))
        [ensemble,ncut]=respectBoundsWithReport(ensemble,kalmanOptions,options);
        disp(['average number of truncate per realization is ',num2str(ncut)]);
        fprintf(fidAGSMDA,'%s %f \n','average number of truncate per realization is ',ncut);
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
    %simDataOld=simData;
    
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
    fprintf(fidAGSMDA,'%s \n',' ');
    fprintf(fidAGSMDA,'%s %d \t %f \t %f \t  %f \n','AGSMDA ', ...
        iter,lambda(iter),objNew,objStdNew);
    save(strcat(dir,'objRealIter',num2str(iter)),'objReal');
    
    % save files
    % 1) entire ensemble of model varialbes
    save(strcat(dir,'ensemble',num2str(iter)),'ensemble','weights');
    
%     simDataOld=simData;
%     ensembleOld=ensemble;
%     objStd=objStdNew;
%     obj=objNew;
    
%     if ~exist([dir,'ensemble',num2str(iter)],'file')
%         save(strcat(dir,'ensemble',num2str(iter)),'ensemble','lambda','weights');
%     end
%     
%     if ~exist([dir,'simulatedDataIter',num2str(iter)],'file')
%         save(strcat(dir,'simulatedDataIter',num2str(iter)),'simData','numberOfFailedSimulations');
%     end
%     
%     if ~exist([dir,'objRealIter',num2str(iter)],'file')
%         save(strcat(dir,'objRealIter',num2str(iter)),'objReal');
%     end
    
   iter=iter+1;
      
end

% Close debugDisLoc file
if isfield(kalmanOptions,'distanceLoc')
    fclose(fidLoc);
end

fclose(fidAGSMDA);

end
