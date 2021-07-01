function []=RLM_MAC(iter,ensemble,simData,W,Wbase,H,measurement,scale,options,kalmanOptions,dir,obsLocation,obsType)

% RLM_MAC implemented for the case "prodSeis"
% Iterative ensemble smoother with minimum-average-cost
% function []=RLM_average_cost(iter,ensemble,simData,nm,W,Wbase,H,measurement,options,kalmanOptions,dir,obsLocation,obsType)
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
% kalmanOptions:        Structure defining options for the EnKF parameters, see
%                       defaultKalmanOptions.m for more info.
% dir:                  work directory
% obsLocation:          observation location (for covariance localization)
% obsType:              observation type (for covariance localization)
%
%
% Copyright (c) 2010-2016 IRIS, All Rights Reserved.
% $Id: //depot/rfmatlab/main/Kalman/RLM_average_cost.m#18 $
% $DateTime: 2017/05/18 15:13:28 $

% Truncate debug files only if iter == 0
debugFilesMode = 'a';
if iter == 0
    debugFilesMode = 'w';
end

sWbase = sqrtW(Wbase);

fidRLM=fopen('debugRLM_MAC', debugFilesMode);
tline = 'iteration number ';
fprintf(fidRLM,'%s %d\n',tline,iter);


if nargin<12
    if getOption(kalmanOptions,'distanceLoc',false)
        error('RLM_average_cost.m: if using distanced-based localizaiton, need to specify obsLocation and obsType')
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
    error('RLM_avarage_cost.m: the dimension of H and measurement is not consistent!')
end

if size(W,1)~=size(measurement,1)
    error(['RLM_average_cost.m: the dimension of W and measurement is not', ...
        ' consistent!'])
end

if isfield(kalmanOptions,'useSeismic') && kalmanOptions.useSeismic &&...
        kalmanOptions.numParProcesses == 1
    error('RLM_avarage_cost.m: Seismic only suppoted in parallel runs')
end

if isfield(kalmanOptions,'append_mean') && kalmanOptions.append_mean
    ne=size(ensemble,2) - 1; % ensemble size (excluding ensemble mean)
else
    %ne=size(ensemble,2); % ensemble size
    error('Please set kalmanOptions.append_mean = 1 ...')
end
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

iterLambda=0;
if iter > 0
    load(['ensemble',num2str(iter),'.mat'],'lambda');
end
disp(' ');
disp('--------------------------------------------');

if isfield(kalmanOptions,'useSeismic') && kalmanOptions.useSeismic && ...
        isfield(kalmanOptions,'useProd')    && kalmanOptions.useProd
    
    [obj,objStd,objReal,mismatchMtx,individualObjInfo]=getDataMismatch(simData(:,1:ne),W,measurement);
    disp(['obj = ',num2str(obj),' objStd = ',num2str(objStd)]);
    
    disp(['obj_prod =',num2str(individualObjInfo{3}(1)),' objStd_prod = ',num2str(individualObjInfo{3}(2))]);
    disp(['obj_seis =',num2str(individualObjInfo{4}(1)),' objStd_seis = ',num2str(individualObjInfo{4}(2))]);
    disp('--------------------------------------------');
    disp(' ');
    
    fprintf(fidRLM,'%s %s %s \t %s %s %s \t %s %s %s \t  %s %s \n','rlm-mac ', ...
        'iter','iterLambda','obj','(obj_prod,','obj_seis)','objStd',...
        '(objStd_prod','objStd_seis)','lambda','changeM');
    fprintf(fidRLM,'%s \n',' ');
    fprintf(fidRLM,'%s %d %d \t %f %s %f %s %f %s \t %f %s %f %s %f %s \n',...
        'rlm-mac ', iter,iterLambda,obj,'(',individualObjInfo{3}(1),',',...
        individualObjInfo{4}(1),')',objStd,'(',individualObjInfo{3}(2),',',...
        individualObjInfo{4}(2),')');
    fprintf(fidRLM,'%s \n',' ');
    save(['objRealIter',num2str(iter)],'objReal','individualObjInfo','mismatchMtx');
    %load('trueJoindDataSmoother.mat','num_prod','num_seis'); % use num_prod and num_seis to d
    
else
    
    [obj,objStd,objReal]=getDataMismatch(simData(:,1:ne),W,measurement);
    disp(['obj=',num2str(obj),' objStd=',num2str(objStd)]);
    fprintf(fidRLM,'%s %s \t %s \t %s \t %s \t  %s\t %s \n','rlm-mac ', ...
        'iter','iterLambda','objNew','objStdNew','lambda','changeM');
    fprintf(fidRLM,'%s \n',' ');
    fprintf(fidRLM,'%s %d \t %d \t \t %f \t %f \n','rlm-mac ', iter,iterLambda, ...
        obj,objStd);
    fprintf(fidRLM,'%s \n',' ');
    save(['objRealIter',num2str(iter)],'objReal');
    
end


if getOption(kalmanOptions,'distanceLoc',false) % distance-based localization
    fidLoc=fopen('debugDisLoc', debugFilesMode);
    tline = 'starting iteration';
    fprintf(fidLoc,'%s\n',tline);
end

smallReduction=0;
dm_threshold = (kalmanOptions.beta)^2 * length(measurement);
% start the loop for iteration
while (iter < kalmanOptions.maxIter) && (obj > dm_threshold)
    
    if getOption(kalmanOptions,'distanceLoc',false) % distance-based localization
        tline = 'checking for data with zero variability from ensemble prediction';

        fprintf(fidLoc,'%s\n',tline);
        fprintf(fidLoc,'%s %d\n','number of obs before screening ',length(measurement));
    end
    
    fprintf(fidRLM,'%s %d\n','number of data is ',size(simData,1));
    disp(' ');
    disp('--------------------------------------------');
    disp(['-- Iteration step: ' int2str(iter) ' --']);
    disp('--------------------------------------------');
    disp(' ');
    disp(['number of data is ',num2str(size(simData,1))])
    numberOfKeptMeasurements=1:length(measurement);
    
    simMean = simData(:,end); % simMean has already been normalized by C_{D}^{-1/2}
    
    deltaD=zeros(nd,ne);
    for j=1:ne
        deltaD(:,j)=simData(:,j)-simMean;
    end
    
    if isfield(kalmanOptions,'ignoreUninformativeMeasurements') && ...
            kalmanOptions.ignoreUninformativeMeasurements==1
        
        varTestMeasurements=var(deltaD,0,2)';
        %keepVariables=(varTestMeasurements~=0);
        keepVariables=find(varTestMeasurements~=0);
        numberOfKeptMeasurements=numberOfKeptMeasurements(keepVariables);
        deltaD=deltaD(keepVariables,:);
        simData=simData(keepVariables,:);
        perturbedDataEff=perturbedData(keepVariables,:);
        
        %save(['keepVariables_Iter' num2str(iter) '.mat'],'keepVariables');
        
        if getOption(kalmanOptions,'distanceLoc',false) % distance-based localization
           
            % need to exclude statistics obs. (used in facies studies 
            tmp_options = load('inputData.mat','options');
            if isfield(tmp_options.options,'biModel') && (~isempty(tmp_options.options.biModel)) % for facies studies, if additional channel statistics are used as obs
                load('trueSolutionSmoother','real_obs_index','obs_type_index');
                
                keepVariables_realObs = intersect(real_obs_index,keepVariables)';
                tmp_index = find(real_obs_index == keepVariables_realObs);
                kept_obs_type_index = obs_type_index(tmp_index);
                
                obsTypeEff=obsType(kept_obs_type_index,:);
                obsLocationEff=obsLocation(kept_obs_type_index,:);
            else
                load('trueSolutionSmoother','num_prod');
                real_obs_index = (1 : length(H(:,1)))'; % used in "localization.m"
                keepProdVariables = find(varTestMeasurements(1:num_prod)~=0);
                kept_obs_type_index = keepProdVariables; % used in "localization.m"
                tmp_index = keepProdVariables;
                obsTypeEff=obsType(keepProdVariables,:);
                obsLocationEff=obsLocation(keepProdVariables,:);
            end
            clear tmp_options;
            
            fprintf(fidLoc,'%s %d \n','number of obs after screening ',size(obsTypeEff,1));
        end
        disp(['number of data for updating is ',num2str(size(simData,1))])
        
        fprintf(fidRLM,'%s %s\n','number of data for updating is ',num2str(size(simData,1)));
        
        % after this we continues as before, except we fill in zero
        % columns in Ke before reporting, except if all measurements are
        % removed:
        if sum(keepVariables)==0
            warning('rfmatlab:Kalman:RLM_average_cost','no variability in ensemble prediction');
            return;
        end
    else
        perturbedDataEff = perturbedData;
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
    fprintf(fidRLM,'%s %d\n','number of singular value retained ',svdPd);
    Vd=Vd(:,1:svdPd);
    Ud=Ud(:,1:svdPd);
    Wd=val(1:svdPd); % a vector
    clear val;
    
    iterLambda=1;
    
    if iter == 0
        lambda = kalmanOptions.lambda;
    end
    
    while iterLambda < kalmanOptions.maxInnerIter
        % workflow in the RLM_average_cost
        %       1) update estimates
        %       2) runForwardSimulation
        %       3) check if average data mismatch is reduced
        %       4) if yes do next outer iteration unless a certain
        %          stopping condition is met
        %       5) if not, do inner iterations to search for better
        %          estimates
        %       6) if too many inner loops, exit the inner loop
        %          and do next outer iteration
        
        alpha = lambda * sum(Wd.^2) / svdPd;
        %alpha = lambda * sum(Wd) / svdPd;
        %alpha = lambda .* (Wd.^2);
        alpha = max(alpha,1e-3); % in case of ensemble collapge
        
        x1=zeros(ne,svdPd);
        for jVec=1:svdPd
            for j=1:ne
                x1(j,jVec)=Vd(j,jVec)*Wd(jVec)/(alpha+Wd(jVec)^2);
                %x1(j,jVec)=Vd(j,jVec)*Wd(jVec)/(alpha(jVec)+Wd(jVec)^2);
            end
        end
        
        % do the localization
        X = (simData(:,1:ne)-perturbedDataEff);
        I = localization(ensemble,x1,Ud,X,obsTypeEff,obsLocationEff);
        
        %update
        ensembleOld=ensemble; % here ensemble size = ne + 1 (1 for ensemble mean)
        ensemble(:,1:ne)=ensemble(:,1:ne) - I;
        
        % append ensemble mean
        ensemble(:,ne+1)= mean(ensemble(:,1:ne),2);
        
        % re-initialize the facies field in cases of facies studies
        if isfield(options,'postProcessEnsemble') && ...
                ~strcmp(options.postProcessEnsemble,'')
            [ensemble,~] = feval(options.postProcessEnsemble,ensemble);
        end
        
        changeM=sqrt(sum((ensemble(:,end)-ensembleOld(:,end)).^2));
        disp(['average change to model variable ',num2str(changeM)]);
        
        % put the ensemble members within the bounds:
        check_adjustment = (isfield(kalmanOptions,'staticVarLB') && ~isempty(kalmanOptions.staticVarLB)) || ...
            (isfield(kalmanOptions,'staticVarUB') && ~isempty(kalmanOptions.staticVarUB)) || ...
            (isfield(kalmanOptions,'freeparamLB') && ~isempty(kalmanOptions.freeparamLB)) || ...
            (isfield(kalmanOptions,'freeparamUB') && ~isempty(kalmanOptions.freeparamUB));
        if ~isfield('options',' writePerm') && check_adjustment
            [ensemble,ncut]=respectBoundsWithReport(ensemble,kalmanOptions,options);
            disp(['average number of truncate per realization is ',num2str(ncut)]);
            fprintf(fidRLM,'%s %f \n','average number of truncate per realization is ',ncut);
        end
        
        options.existTimeSteps=0;
        simDataOld=simData;
        
        numberOfFailedSimulations = [];
        filecontents = [];
        if kalmanOptions.numParProcesses<2
            [simulatedEnsemble,filecontentsOut,numberOfFailedSimulations]= ...
                runForwardSimulations(ensemble,0,kalmanOptions,options,[], ...
                filecontents);
            jtmp=size(simulatedEnsemble,1)*size(simulatedEnsemble,2); % number of total data
            itmp=size(simulatedEnsemble,3); % number of realizations
            simulatedEnsemble=reshape(simulatedEnsemble,jtmp,itmp);
            simData=simulatedEnsemble(H(:,2),:);
        else
            [simData,~]= runMultiCoreSim(0,ensemble,H,kalmanOptions,options);
        end
        save(strcat(dir,'tmpSimData.mat'),'-v7.3','simData','simulatedEnsemble','filecontentsOut'); % save latest simData before normalizing
        
        if ~isempty(scale)
            simData=simData.*repmat(scale,1,ne+1);
        end
        simData=normalizeData(simData, sWbase);
        
        %%%%%%%%%%%%%%%%%%%%%%%%
        
        disp(' ');
        disp('--------------------------------------------');
        
        if isfield(kalmanOptions,'useSeismic') && kalmanOptions.useSeismic && ...
        isfield(kalmanOptions,'useProd')    && kalmanOptions.useProd
    
            [objNew,objStdNew,objReal,mismatchMtx,individualObjInfo]=getDataMismatch(simData(:,1:ne),W,measurement); %#ok<*ASGLU>
            disp(['obj = ',num2str(objNew),' objStd = ',num2str(objStdNew)]);
            disp(['obj_prod =',num2str(individualObjInfo{3}(1)),' objStd_prod = ',num2str(individualObjInfo{3}(2))]);
            disp(['obj_seis =',num2str(individualObjInfo{4}(1)),' objStd_seis = ',num2str(individualObjInfo{4}(2))]);
            disp('--------------------------------------------');
            disp(' ');
            
            %--
            fprintf(fidRLM,'%s \n',' ');
            fprintf(fidRLM,'%s %d %d \t %f %s %f %s %f %s \t %f %s %f %s %f %s %s %s \n',...
                'rlm-mac ', iter,iterLambda,objNew,'(',individualObjInfo{3}(1),',',...
                individualObjInfo{4}(1),')',objStdNew,'(',individualObjInfo{3}(2),',',...
                individualObjInfo{4}(2),')',lambda,changeM);
            fprintf(fidRLM,'%s \n',' ');
            %save(['objRealIter',num2str(iter),'-',num2str(iterLambda)],'objReal','individualObjInfo','mismatchMtx');
            
        else
            
            [objNew,objStdNew,objReal]=getDataMismatch(simData(:,1:ne),W,measurement);
            disp(' ');
            disp('   --------------------------');
            disp(['   objNew=',num2str(objNew),' objStdNew=',num2str(objStdNew)]);
            disp('   --------------------------');
            disp(' ');
            
            %--
            fprintf(fidRLM,'%s \n',' ');
            fprintf(fidRLM,'%s %d \t %d \t \t %f \t %f \t  %f \t %f \n','rlm-mac ', ...
                iter,iterLambda,objNew,objStdNew,lambda,changeM);
            fprintf(fidRLM,'%s \n',' ');
            %save(['objRealIter',num2str(iter),'-',num2str(iterLambda)],'objReal');
            
        end
        
        %%%%%%%%%%%%%%%%%%%%%%%%
        
        if objNew>obj
            lambda=lambda*kalmanOptions.lambda_increment_factor;
            disp(['increasing Lambda to ',num2str(lambda)]);
            iterLambda=iterLambda+1;
            simData=simDataOld;
            ensemble=ensembleOld;
        else
            changeStd=(objStdNew-objStd)/objStd;
            disp(['changeStd=',num2str(changeStd)]);
            
            lambda=lambda*kalmanOptions.lambda_reduction_factor;%
            
            disp(['reducing Lambda to ',num2str(lambda)]);
            
            if abs(objNew-obj)/abs(obj)*100<kalmanOptions.minReduction
                smallReduction=1;
            end
            
            simDataOld=simData;
            ensembleOld=ensemble;
            objStd=objStdNew;
            obj=objNew;
            break % break the inner loop over lambda
        end
        
    end  % end of inner loop
    
    % if a better update not found
    if iterLambda >= kalmanOptions.maxInnerIter 
        
        lambda = lambda * kalmanOptions.lambda_increment_factor;
        if lambda < kalmanOptions.lambda
            lambda = kalmanOptions.lambda;
        end
        
        tline = 'terminating iterations: iterLambda>=maxInnerIter';
        disp(tline);
        fprintf(fidRLM,'%s \n',tline);
        
    end
    
    % save
    iter=iter+1;
    save(strcat(dir,'ensemble',num2str(iter)),'ensemble','lambda');
    movefile(strcat(dir,'tmpSimData.mat'),strcat(dir,'simulatedDataIter',num2str(iter),'.mat'));
    save(strcat(dir,'objRealIter',num2str(iter)),'objReal');
    if isfield(kalmanOptions,'useSeismic') && kalmanOptions.useSeismic && ...
        isfield(kalmanOptions,'useProd')    && kalmanOptions.useProd
        save(strcat(dir,'objRealIter',num2str(iter)),'individualObjInfo','-append');
    end
        
    if smallReduction==1
        tline = ['terminating iterations: reduction of objective function is less than ', ...
            num2str(kalmanOptions.minReduction),'%'];
        disp(tline);
        fprintf(fidRLM,'%s \n',tline);
        break
    end
    
end  % end of loop for iteration

if iter >= kalmanOptions.maxIter
    tline = 'terminating iterations: iter>=kalmanOptions.maxIter';
    disp(tline);
    fprintf(fidRLM,'%s \n',tline);
end

% Close debugDisLoc file
if getOption(kalmanOptions,'distanceLoc',false)
    fclose(fidLoc);
end

fclose(fidRLM);


