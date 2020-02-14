function [simData,simulatedEnsemble,filecontentsOut,numberOfFailedSimulations]= ...
    runMultiCoreSim(numStep,ensemble,H,kalmanOptions,options,filecontents,varargin) %#ok<*INUSL>

% Run forecasts on a multi-core machine with linux os.
%
% Input
% -----
% numStep       : step to start from
% ensemble      : ensemble of static variables
% H             : observation operator for iES
% kalmanOptions : options related to the Kalman filter
% options       : options for the reservoir simulator
% filecontents  : contents of restart file
%
% Output
% ------
% simData                   : ensemble of simulated data
% simulatedEnsemble         : ensemble of static and state variables after
%                             simulation
% fileContentsOut           : contents of next restart file
% numberOfFailedSimulations : list with the simulation numbers that fail
%
% Copyright(c) International Research Institute of Stavanger (IRIS). All Rights Reserved.
% $Id: //depot/rfmatlab/main/eclipseKalman/runMultiCoreSim.m#7 $
% $DateTime: 2017/01/12 12:08:59 $
%


numStep = setProperty(varargin,'numStep',numStep);% step to start from
skip_preparation = setProperty(varargin,'skip_preparation',0);
skip_simulation = setProperty(varargin,'skip_simulation',0);
skip_collection = setProperty(varargin,'skip_collection',0);

% write eclipse input for each ensemble member
basefile = options.filename;

% initialize eclipse input
numberOfFailedSimulations = [];
if nargin < 6
    filecontents = [];
end

% **********************for debug only*******************
% disp('Debug: No writing inputs in Sim folder');
% options.existTimeSteps=0;
% **********************for debug only*******************
if ~skip_preparation
    if isfield(kalmanOptions,'useSeismic') && kalmanOptions.useSeismic
        keyword = 'PORO';
        SUCCESS = 0;
        for i = 1 : size(options.staticVar,1)
            if strcmp(keyword,deblank(options.staticVar(i,:)))
                SUCCESS = i;
                break;
            end
        end
        if ~SUCCESS
            error('Keyword is not found')
        end
        numFreeParam = 0;
        if isfield(options,'freeparam')
	        numFreeParam = sum(options.freeparam);
        end
        Poro_index = numFreeParam + (1:options.numGridBlocks) + (i-1)*options.numGridBlocks;
    end

    for I = 1:size(ensemble,2)
        % go to dir
        disp(['Write into Sim',num2str(I)])
        cd(['Sim',num2str(I)]);
        
        
        options.filename = [basefile,'_',num2str(I)];
        if ~isempty(filecontents)
            extforwecl(ensemble(:,I),numStep,filecontents(I),options,'init');
        else
            extforwecl(ensemble(:,I),numStep,[],options,'init');
        end
        if I==1
            options.existTimeSteps=1;
        end
        if isfield(kalmanOptions,'useSeismic') && kalmanOptions.useSeismic
            active_Poro = ensemble(Poro_index,I);
            save active_Poro.mat active_Poro;
        else
            if exist('./active_Poro.mat','file')
                delete ./active_Poro.mat;
            end
        end
        % move up
        cd('..');
    end
end

if ~skip_simulation
    % run eclipse on multiple cores using ppss
    unix('rm ppss*.txt');
    unix('rm ppss.sh_is_running');
    unix('rm -f JOB_LOG/*');
    disp('Running runMultiCore.sh');
    unix('./runMultiCore.sh');
end

% collect results
if ~skip_collection
    disp('Collecting results...');
    simData = [];
    simulatedEnsemble = [];
    for I = 1:size(ensemble,2)
        disp(['Collecting results in Sim',num2str(I)]);
        % go to dir
        cd(['Sim',num2str(I)]);
        
        % read eclipse output
        options.filename = [basefile,'_',num2str(I)];
        if ~isempty(filecontents)
            [simulatedEnsembleTemp,~,filecontentsOut(I)] = extforwecl(ensemble(:,I),numStep,filecontents(I),options,'end');
        else
            [simulatedEnsembleTemp,~,filecontentsOut(I)] = extforwecl(ensemble(:,I),numStep,[],options,'end');
        end
        nt = size(simulatedEnsembleTemp,2);
        jtmp=size(simulatedEnsembleTemp,1)*size(simulatedEnsembleTemp,2); % number of total data
        itmp=size(simulatedEnsembleTemp,3); % number of realizations
        simulatedEnsembleTemp = reshape(simulatedEnsembleTemp,jtmp,itmp);
        
        % production data
        if (isfield(kalmanOptions,'iterES') && kalmanOptions.iterES==1) && ...
                (~isfield(kalmanOptions,'useProd') || kalmanOptions.useProd == 1)
            if isempty(H)
                error('H cannot be empty')
            end
            prodData=simulatedEnsembleTemp(H(:,2),:);
        else
            prodData=[];
        end
        simulatedEnsemble = [simulatedEnsemble simulatedEnsembleTemp];
        clear simulatedEnsembleTemp;
        
        if isfield(kalmanOptions,'useSeismic') && kalmanOptions.useSeismic
            load('./simSeismicData.mat','simSeismicData');
        else
            simSeismicData = [];
        end
        
        % move up
        cd('..');
        simData = [simData [prodData;simSeismicData] ]; %#ok<*AGROW>
    end
end

% reshape
if (length(options.timespec)-numStep > 2) || (isfield(kalmanOptions,'iterES') && kalmanOptions.iterES==1)
    simulatedEnsemble = reshape(simulatedEnsemble,size(simulatedEnsemble,1)/nt,nt,size(ensemble,2));
end
