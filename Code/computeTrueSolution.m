function computeTrueSolution()

% function computeTrueSolution()
%
% Compute the "true" solution of a data assimilation problem.
% The simulator is eclipse.
%
% The input data is found from the file inputData.mat (and the
% necessary eclipse files).
%
% kalmanOptions - fields:
%  .timespec          - time for reports.
%  .trueSimName       - name of eclipse file used to generate true
%                       data.
%  .trueScript        - name of eclipse file used to generate true
%                       data.
%
% Copyright(c) International Research Institute of Stavanger (IRIS)
% $Id: //depot/rfmatlab/main/eclipseKalman/computeTrueSolution.m#10 $
% $DateTime: 2018/03/16 17:38:31 $

% set default random number generator
rng('default')

load('inputData','options','kalmanOptions','trueStaticVar')
options=defaultOptions(options); %#ok<*NODEF>
kalmanOptions=defaultKalmanOptions(kalmanOptions);

% help people remember that we now only save static values for the active
% cells
if exist('trueStaticVar','var') &&options.numGridBlocks < options.fieldSize && isequal(size(trueStaticVar,1),options.fieldSize* size(options.staticVar,1)) && ...
        (length(options.staticVar(1,:)) < 6 || ~strcmp(options.staticVar(1,1:6),'FACIES'))
    disp('trueStaticVar should only contain values for active cells')
    disp('this is fixed for you')
    trueStaticVar=reshape(trueStaticVar,options.fieldSize,size(options.staticVar,1),1);
    trueStaticVar=trueStaticVar(options.actnum==1,:);
    trueStaticVar=reshape(trueStaticVar,options.numGridBlocks*size(options.staticVar,1),1);
    save('inputData.mat','-append','trueStaticVar')
end

% this is a prediction, we use the options.prediction to make sure
% that the results obtained from the eclipse run is collected:
options.prediction=1;

% we do not use the append/APPEND options while predicting, this is
% a single run:
if isfield(options,'APPEND')
    options=rmfield(options,'APPEND');
end
if isfield(options,'append')
    options=rmfield(options,'append');
end

% the time steps that we are going to use are defined in
% kalmanOptions.timespec.
options.timespec=kalmanOptions.timespec;


% run with true solution - if it is not computed already
if ~existfile('./trueSolution.mat')
    %disp('Compute true solution')
    [F]=computeTrueSol(kalmanOptions,options,trueStaticVar);
    if F==0
        return
    end
else
    disp(['file trueSolution.mat already exist, so no simulation is', ...
        ' needed!'])
end


% HELP FUNCTIONS - TO SIMPLIFY CODE:
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% The function below is from computePredictedSolution - except that
% the function name is changed.

function [F]=computeTrueSol(kalmanOptions,options,trueStaticVar)

% function [F]=computeTrueSol(kalmanOptions,options,trueStaticVar)
%
%
% Copyright(c) RF-Rogaland Research 2002. All rights reserved.


if ~existfile('./forwSim.mat')
    disp('Compute true solution')
    % the "true" simulation may use a different filename than the
    % simulations in the filter.
    if isfield(kalmanOptions,'trueSimName')
        options.filename=kalmanOptions.trueSimName;
    end
    if isfield(kalmanOptions,'trueScript')
        options.specialScript=kalmanOptions.trueScript;
    end
    
    % initial data:
    numstep=0;
    filecontentsStart=[];
    
    [eclOut,obstimes,filecontents]=extforwecl(trueStaticVar,numstep, ...
        filecontentsStart, options);

     if isfield(kalmanOptions,'iterES') && kalmanOptions.iterES==1
         % eclOut = [ y ]

         state = [];
         y = eclOut;
     else
         % eclOut = [ static; state; y ]

         numStatic = 0;
         if options.returnStatic==1
             numStatic = length(trueStaticVar);
         end

         % Might be unsafe?
         numState = size(options.statevar, 1) * options.numGridBlocks;

         state = eclOut((1 + numStatic):(numStatic + numState), :);
         y = eclOut((1 + numStatic + numState):end, :);
     end

    if isfield(options,'artsimVar')

        %running ARTSim
        [ArtsimState,ArtsimY] = runARTSim(options,filecontents);
        % collect states from both simulators
        state = [state;ArtsimState];
        % collect well measurements from both simulators
        y = [y;ArtsimY];
        
    end
    
    if kalmanOptions.Hones==0
        if kalmanOptions.saveMemory==0
            predictedState=[kron(trueStaticVar,ones(1,size(obstimes,2)));...
                state; y];
        else
            warning('rfmatlab:eclipseKalman',['computeTrueSolution is not yet fully updated to work with ' ...
                'kalmanOptions.saveMemory~=0'])
        end
        save('forwSim','state','filecontents','numstep','y','obstimes', ...
            'filecontentsStart','predictedState')
        disp('saving results of simulation in forwSim.mat')
    else
        predictedState=[kron(trueStaticVar,ones(1,size(obstimes,2)));...
            state; y];
        save('forwSim','state','filecontents','numstep','y','obstimes', ...
            'filecontentsStart','predictedState')
        disp('saving results of simulation in forwSim.mat')
    end
    
end

% since treatment of measurements and uncertainties are very example
% dependent we provide the following possibility for specially designed
% functions:
if isfield(kalmanOptions,'measFunc')
    nin = nargin(kalmanOptions.measFunc);
    nout = nargout(kalmanOptions.measFunc);
    if nin == 1 && nout == 2
        eval(['[W,H] = ',kalmanOptions.measFunc,'(kalmanOptions.keywords);'])
    elseif nin == 1 && nout == 0
        eval([kalmanOptions.measFunc,'(kalmanOptions.keywords)'])
    elseif nin == 0 && nout == 2
        eval(['[W,H] = ',kalmanOptions.measFunc,';'])
    else
        eval(kalmanOptions.measFunc);
    end
end

% Calculate or specify H and W
if isfield(kalmanOptions,'buildMatricesOldWay') && ...
        kalmanOptions.buildMatricesOldWay == 1
    % build measurement matrices
    load('forwSim.mat');
    numVar = size(predictedState,1);
    [W,H]=buildMeasurementMatrices(obstimes, ...
        kalmanOptions,numVar,y);
    save('forwSim','-append','H','W')
else
    load('forwSim')
    if exist('trueSolution.mat','file')
        load trueSolution.mat
    end
    if ~exist('H','var') && ~exist('W','var')
        disp('H and W must be spesified')
        F=0;
        return
    end
end
F=1;

% compute measurements:
if ~exist('trueSolution.mat','file')
    
    load('forwSim','state','filecontents','numstep','y','obstimes','predictedState','filecontentsStart')
    
    %preallocate
    truemeasurement=cell(1,size(predictedState,2));
    measurement=cell(1,size(predictedState,2));
    
    for I=1:size(predictedState,2)
        if isempty(H{I})
            truemeasurement{I} = [];
        else
            if isfield(kalmanOptions,'Hones') && kalmanOptions.Hones==1
                truemeasurement{I} = predictedState(H{I}(:,2),I);
            else
                truemeasurement{I} = H{I}*predictedState(:,I);
            end
        end
    end
    
    numStaticMeas=kalmanOptions.numStaticMeas;
    
    if isfield(kalmanOptions,'historicalData') && ...
            kalmanOptions.historicalData==1
        measurement=truemeasurement; %#ok<NASGU>
    else
        if isfield(kalmanOptions,'state')
            s = 1974+kalmanOptions.state;
            while s > 4.2897e+09, s = round(s/1000); end
            rng(s);
        end
        for I=1:size(truemeasurement,2)
            measurement{I}=addgnoise(truemeasurement{I},W{I});
        end
        if numStaticMeas>0
            staticMeasurement=measurement{1}(1:numStaticMeas);
        end
        % if truemeasurement=0 it is likely that there is no
        % observation, therefore we set measurement=0:
        for I=1:size(truemeasurement,2)
            measurement{I}(truemeasurement{I}==0)=0;
            % non negative dynamic measurements are unphysical in this setting:
            measurement{I}(measurement{I}<0)=0;
        end
        if numStaticMeas>0
            for I=2:size(truemeasurement,2)
                measurement{I}(1:numStaticMeas)=staticMeasurement;
            end
        end
    end
    
    if isfield(kalmanOptions,'trueSimName')
        SMStrue=readXfile(strcat(kalmanOptions.trueSimName,'.SMSPEC'));
    else
        SMStrue=readXfile(strcat(options.filename,'.SMSPEC'));
    end
    
    if exist('measurementIndexTable','var')
        save('trueSolution','predictedState','filecontents','obstimes', ...
            'filecontentsStart','truemeasurement','measurement','SMStrue',...
            'H','W','measurementIndexTable');
    else
        save('trueSolution','predictedState','filecontents','obstimes', ...
            'filecontentsStart','truemeasurement','measurement','SMStrue',...
            'H','W');
        
    end
    
end

if isfield(kalmanOptions,'distanceLoc') && kalmanOptions.distanceLoc == 1
    
    kalmanOptions = Localc_2D();
    %         load('trueSolution')
    %         if ~exist('H','var') || ~exist('obstimes','var')
    %             error('H and obstimes must be specifyed and saved in trueSolution.mat')
    %         end
    
    if ~exist('obsLocation','var') || ~exist('obsType','var')
        
        if ~exist('SMStrue','var')
            if isfield(kalmanOptions,'trueSimName')
                SMStrue=readXfile(strcat(kalmanOptions.trueSimName,'.SMSPEC'));
            else
                SMStrue=readXfile(strcat(options.filename,'.SMSPEC'));
            end
        end
        
        if ~isfield(kalmanOptions,'pickObs')
            numObs = length(obstimes);
        else
            numObs = length(kalmanOptions.pickObs);
        end
        obsLocation=cell(1,numObs);
        obsType=cell(1,numObs);
        
        for I=1:numObs
            measLoc=SMStrue.WGNAMES(H{I}(:,2)-kalmanOptions.index(end)+2,:);
            measType=SMStrue.KEYWORDS(H{I}(:,2)-kalmanOptions.index(end)+2,:);
            obsLocation{I}=measLoc;
            obsType{I}=measType;
        end
        
        save('trueSolution','-append','obsLocation','obsType')
    end
    
end

% compute true solution for seismic
if isfield(kalmanOptions,'useSeismic') && kalmanOptions.useSeismic == 1 
% && (~isfield(kalmanOptions,'historicalData') || kalmanOptions.historicalData == 0)
    disp('Computing seismic data');
    generateSeismicData()
end

% reformating 'measurement','W','H','obsLocation','obsType'
if isfield(kalmanOptions,'iterES') && kalmanOptions.iterES == 1
    
    numGridBlocks = options.numGridBlocks;
    numStaticVar = size(options.staticVar,1);
    numStateVar = size(options.statevar,1);
    numFreeParam = 0;
    if isfield(options,'freeparam')
        numFreeParam = sum(options.freeparam);
    end
    nm = numFreeParam + numStaticVar*numGridBlocks; % number of parameters (including free parameters)
    ns = numStateVar*numGridBlocks; % number of dynamical state variables
    
    if ~isfield(kalmanOptions,'pickObs')
        P = 1:length(obstimes);
    else
        P = kalmanOptions.pickObs;
    end
    
    obsLocation=[];
    obsType=[];
    num_prod = 0;
    num_seis = 0;
    if isfield(kalmanOptions,'useSeismic') && kalmanOptions.useSeismic
        load('trueSeismicData','measurement','truemeasurement','W',...
            'cD_leading_index','cD_leading_number','cA_leading_index','cA_leading_number');
        load('inputData.mat','seismicOptions');
        meas=measurement;
        truemeas=truemeasurement;
        WW=W;
        measurement=[]; truemeasurement=[]; W=[]; H=[];
        for I=1:length(meas)
            measurement=[measurement; meas{I}];
            truemeasurement=[truemeasurement; truemeas{I}];
            W=concatenateW(W, WW{I}); % Even seismic W sould be variance, not std
        end
        num_seis = length(measurement);
        
        clear meas truemeas WW;
        
        if isfield(kalmanOptions,'useProd') && kalmanOptions.useProd
            seis_W = W;
            seis_measurement = measurement;
            seis_truemeasurement = truemeasurement;
            clear W measurement truemeasurement;
        end
    end
    
    if isfield(kalmanOptions,'useProd') && kalmanOptions.useProd
        load('trueSolution','measurement','truemeasurement','W','H','y');
        if isfield(kalmanOptions,'distanceLoc') && kalmanOptions.distanceLoc == 1
            load('trueSolution','obsLocation','obsType');
            measLoc = obsLocation;
            measType = obsType;
            obsLocation = [];
            obsType = [];
        end
        meas=measurement;
        truemeas=truemeasurement;
        WW=W;
        HH=H;
        measurement=[]; truemeasurement=[]; W=[]; H=[];
        for I=1:length(P)
            measurement=[measurement; meas{I}];
            truemeasurement=[truemeasurement; truemeas{I}];
            W=concatenateW(W, WW{I});
            HHtmp=HH{I};
            HHtmp(:,1)=HHtmp(:,1)+size(H,1);
            HHtmp(:,2)=HHtmp(:,2) + size(y,1)*(P(I)-1) - ns -nm;
            H=[H; HHtmp];
            if isfield(kalmanOptions,'distanceLoc') && kalmanOptions.distanceLoc == 1
                obsLocation=[obsLocation; measLoc{I}];
                obsType=[obsType; measType{I}];
            end
        end
        num_prod = length(measurement);
        
        clear meas truemeas WW HH HHtmp y;
        
    end
    
    if  isfield(kalmanOptions,'useSeismic') && kalmanOptions.useSeismic && ...
            isfield(kalmanOptions,'useProd')    && kalmanOptions.useProd
        measurement = [measurement; seis_measurement];
        truemeasurement = [truemeasurement; seis_truemeasurement];
        W=concatenateW(W, seis_W);
        
        clear seis_measurement seis_W seis_truemeasurement;
    end
    
    save('trueSolutionSmoother','measurement','truemeasurement','W','obsLocation','obsType','H','num_prod','num_seis');
    
end


% if isfield(kalmanOptions,'distanceLoc') && kalmanOptions.distanceLoc == 1
%
%     Localc_2D();
%     %    if ~exist('obsLocation','var') || ~exist('obsType','var')
%     obsLocation=cell(size(obstimes));
%     obsType=cell(size(obstimes));
%     for I=1:length(obstimes)
%         measLoc=SMStrue.WGNAMES(H{I}(:,2)-kalmanOptions.index(end)+2,:);
%         measType=SMStrue.KEYWORDS(H{I}(:,2)-kalmanOptions.index(end)+2,:);
%         obsLocation{I}=measLoc;
%         obsType{I}=measType;
%     end
%     %    end
%     save('trueSolution','-append','obsLocation','obsType')
% end
