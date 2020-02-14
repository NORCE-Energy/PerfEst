function setupCase()

disp('Running setupCase.m')

%work_dir = pwd;
%addpath(work_dir);
%if ~exist([work_dir '/localScripts'],'dir')
%mkdir([work_dir '/localScripts'])
%end
%addpath(genpath([work_dir '/localScripts'])); % local MATLAB sripts may be put here
%addpath(genpath([work_dir '/RockphysicsModel'])); % local MATLAB sripts may be put here

warning on
close all;

isReRunTrueSolution = 0; 
if isReRunTrueSolution
    warning off;
    delete forwSim.mat trueSolution.mat trueSolutionSmoother.mat;
    warning on;
end


%%%%%%%%%%%%%%%%%%%%%%%%%%
% Set the options
%%%%%%%%%%%%%%%%%%%%%%%%%%

options=[];

options = getOptions('FULLNORNE.DATA',options);
% options.rftDates = char({'30MAR1998'; '27MAY1998'; '24JUN1998'; '31AUG1998'; '5MAY1999';...
%     '21MAY1999'; '4OCT1999'; '3MAR2000'; '3AUG2000'; '10SEP2000'; '7JUN2001'; '3JUL2001';...
%     '9DEC2003'; '24MAR2005'});
% Dates '2DEC2005'; '1OCT2006' are for Wells B-1AH; K-3H, only exist in
% .DATA file, but no report data exist. Because these 2 wells are drilled
% after the history matching report written in 2005.

% rftStruct.well = [B-4H; D-4H; C-1H; E-3H; C-3H; F-1H; F-2H; E-4H; D-3H; F-3H; B-4AH; F-4H; C-4AH; E-3BH]

% rftDates = getRFTDates('FULLNORNE.DATA',options);



%
% options.statevar=char('PRESSURE','SWAT');
options.statevar=[];
% For seismic calculation
options.dynamicVar=char('PRESSURE','SWAT','SGAS');
options.staticVar = char('PORO','PERMX','NTG');
options.numGridBlocks = sum(options.actnum);
options.historyFile = 'NORPT_NOKH.SCH'; 
options.filename = 'ENKF';
options.specialScript = './runEclipse';
options.useLogPerm = 1;
options.timespec = diff(datenum([options.startDate;options.dates]));
options.unifIn = 0;
options.returnStatic = 1;

% free parameters contain MULTZ (13309), MULTFLT (53), KRW+KRG(8), MULTREGT (27),
% OWC(5); overall = 13383
options.freeparam = [13309, 53, 27, 8, 5]';
options.freeparam_useLog = [1, 1, 1, 0, 0]';
options.freeparamFunc = {'10.^','10.^','10.^','1.*','1.*'};

% options.freeparamFunc = 'norneTransformFreeParam';
%
options.freeparamFileIn = char('multz_template.dat','multflt_template.dat','multregt_template.dat','relperm_template.dat','owc_template.dat');
options.freeparamFileOut = char('multz_simulation.dat','multflt_simulation.dat','multregt_simulation.dat','relperm_simulation.dat','owc_simulation.dat');

% Setup for iES
options.existTimeSteps = 0; % need to set options.existTimeSteps = 0 when using iES

options = defaultOptions(options);

%%%%%%%%%%%%%%%%%%%%%%%%%%
% kalmanOptions
%%%%%%%%%%%%%%%%%%%%%%%%%%

% kalmanOptions.initialEnsemble='initEns'; % containing "ensemble"
% load([kalmanOptions.initialEnsemble,'.mat']);
% kalmanOptions.ensembleSize = size(ensemble,2);

% kalmanOptions.initialEnsemble='initEnsemble300'; % containing "ensemble"
% kalmanOptions.ensembleSize = 300;

% % Only get the first 150 ensembles to runForecast
% 
kalmanOptions.initialEnsemble='initEns'; % containing "ensemble"
kalmanOptions.ensembleSize = 100;
%load initEnsemble300
%ensemble(:,kalmanOptions.ensembleSize+1:end) = [];
%save 'initEns' 'ensemble'
% 
load([kalmanOptions.initialEnsemble,'.mat'],'ensemble');
trueStaticVar=mean(ensemble,2);

kalmanOptions.historicalData=1;

kalmanOptions.staticVarMean=trueStaticVar;
kalmanOptions.ignoreUninformativeMeasurements=0;

% Pick observations
t = options.timespec;
ct = cumsum(t);
dt = 30;
pick = 1;
while pick(end) < length(t)
    v = ct - (ct(pick(end)) + dt);
    [~,ind] = min(abs(v));
    pick = [pick,ind];
end
kalmanOptions.pickObs = pick;
%kalmanOptions.pickObs = 1:length(options.timespec); % Pick all obs for now.

kalmanOptions.timespec=options.timespec;
%kalmanOptions.timespec = 1;
%for numStep = 1:length(kalmanOptions.pickObs)-1
%    kalmanOptions.timespec=[kalmanOptions.timespec...
%        sum(options.timespec(kalmanOptions.pickObs(numStep)+1:kalmanOptions.pickObs(numStep+1)))];
%end

kalmanOptions.state = 2.0427e+05; % Used to control random seed.
kalmanOptions.Hones = 1;
kalmanOptions.trueSimName='FULLNORNE';
kalmanOptions.trueScript='./runTrue';
kalmanOptions = defaultKalmanOptions(kalmanOptions);
kalmanOptions.numParProcesses = 6;
kalmanOptions.saveMemory = 1;
kalmanOptions.keywords = ['WOPRH';'WWPRH';'WGPRH'];
kalmanOptions.measFunc='norneMeasurements';

% Define bounds
dim = ones(options.numGridBlocks,1);
readWoc = [2692.0; 2585.5; 2618.0; 2400.0; 2693.3];
kalmanOptions.freeparamLB = [-4*ones(13309,1); -5*ones(53,1); -5*ones(27,1);...
    0.8*ones(8,1);readWoc-10*ones(size(readWoc))];
kalmanOptions.freeparamUB = [0*ones(13309,1); 2*ones(53,1); 0*ones(27,1);...
    1.5*ones(4,1); 1.0*ones(4,1); readWoc+10*ones(size(readWoc))];
kalmanOptions.staticVarLB = [0.01*dim;0.1*dim;0.01*dim];
kalmanOptions.staticVarUB = [0.4*dim;10*dim;1*dim];

% The bounds of staticVar are calculated based manually HM'ed model
%load staticVarBounds
%kalmanOptions.staticVarLB = staticVarLB;
%kalmanOptions.staticVarUB = staticVarUB; 

% Adjust parameters with bounds
numStaticVar=1:size(options.staticVar,1)*options.numGridBlocks;
if isfield(options,'freeparam')
    numStaticVar=numStaticVar+sum(options.freeparam);
end
% nummeas=options.nummeas;

if isfield(kalmanOptions,'freeparamLB') && isfield(kalmanOptions,'freeparamUB')
    ensemble(1:sum(options.freeparam),:)= ...
      adjustVariableWithInBounds(ensemble(1:sum(options.freeparam),:), ...
				 kalmanOptions.freeparamLB, ...
				 kalmanOptions.freeparamUB);
end

% adjust staticVar within bounds if required:
if isfield(kalmanOptions,'staticVarLB') && isfield(kalmanOptions,'staticVarUB')
  ensemble(numStaticVar,:)= ...
      adjustVariableWithInBounds(ensemble(numStaticVar,:), ...
				 kalmanOptions.staticVarLB, ...
				 kalmanOptions.staticVarUB);
end

disp('Adjusted ensemble and free parameters to be within bounds');
save([kalmanOptions.initialEnsemble,'.mat'],'ensemble');


%%%%%%%%%%%%%%%%%%%%%%%%%%
% Setup for iES
%%%%%%%%%%%%%%%%%%%%%%%%%%
kalmanOptions.iterES = 1;
% iES options can be specified there. So far available options: {'LMEnRML','RLM_average_cost'}
%kalmanOptions.ES_script = 'AGSMDA';
kalmanOptions.ES_script = 'RLM_average_cost';
% for 'RLM_average_cost' only; if isDumpOne = 1, discard 1 model member so that LMEnRML' and 'RLM_average_cost' have equal numbers of forward model runs
kalmanOptions.isDumpOne = 0;

%if ~exist('./reportNum.mat','file')
%    reportNum=resetdate(options,[]);
%    save reportNum.mat reportNum;
%end

if strfind(kalmanOptions.ES_script,'RLM_average_cost')
   kalmanOptions.append_mean = 1; % for RLM-MAC, append ensemble mean to the end of the ensemble, see runIES_multicore.m 
end
%-------------------------------------------------------%
%--Configuration of AGSMDA type iES --%
%-------------------------------------------------------%
kalmanOptions.maxIter = 10; % max number of outer loop iteration
kalmanOptions.maxInnerIter = 5; % max number of inner loop iteration
kalmanOptions.lambda = 1;
%kalmanOptions.numIter = 4;
%kalmanOptions.AGSMDA.adaptive = 1; 
%kalmanOptions.AGSMDA.resample = 1;
%kalmanOptions.AGSMDA.h = [1.0,0.8,0.5,0.2];
%kalmanOptions.AGSMDA.lambda = ones(1,kalmanOptions.numIter)*kalmanOptions.numIter; % initial lambda value
%kalmanOptions.AGSMDA.Nefflim1 = 0;
%kalmanOptions.AGSMDA.Nefflim2 = 1;

kalmanOptions.tsvdData = 0.95; % "energy" threshold that the preserved eigenvalues should sum up to
kalmanOptions.minReduction = 10; % minimum relative change of average data mismatch (in percentage), so 1e-2 actually means .0001
kalmanOptions.retainStaticVarOnly = 1; % only keep the static variables and free parameters

kalmanOptions.lambda_reduction_factor = 0.9; % reduction factor in case to reduce gamma
kalmanOptions.lambda_increment_factor = 2; % increment factor in case to increase gamma 

kalmanOptions.beta = 1; % coeff. in the (squared) data mismatch threshold beta^2 * p (p = observation size) 

kalmanOptions.useProd = 1;
kalmanOptions.useSeismic = 0;

% Adjust measurement noise
kalmanOptions.adjustMeasurementNoise = true;

save('inputData','kalmanOptions','options','trueStaticVar'); 
disp('inputData.mat is created.')

%%%%%%%%%%%%%%%%%%%%%%%%%%
% Localization
%%%%%%%%%%%%%%%%%%%%%%%%%%

kalmanOptions.useLocalization = 1;
kalmanOptions.denoisingLoc = 1;

if isfield(kalmanOptions,'useLocalization') && kalmanOptions.useLocalization
    
    if isfield(kalmanOptions,'denoisingLoc') && kalmanOptions.denoisingLoc == 1
        kalmanOptions.tm = 1; % threshold multiplier for property fields such as PERMX, PORO
        if isfield(options,'freeparam')
            kalmanOptions.freeparam_threshold = 2/sqrt(kalmanOptions.ensembleSize); % (hard-coded) threshold value for free papameters
        end
        kalmanOptions.hybrid_thresholding = 1;
        if kalmanOptions.hybrid_thresholding
            kalmanOptions.hybrid_tm = 2;
        end
    end
    
    if isfield(kalmanOptions,'distanceLoc')  && kalmanOptions.distanceLoc == 1
        
        kalmanOptions.locfunc='GD';
        kalmanOptions.locrange=[1000 1000];
        kalmanOptions.locangle=0;
        kalmanOptions.loctype='WOPR WWPR WGPR';
        
        % table with completion information
        kalmanOptions.wells = {'C-4H ';'B-2H ';'D-1H ';'D-2H ';'B-4H ';'D-4H ';'C-1H ';'E-3H ';'C-2H ';'B-1H ';'C-3H ';'F-1H ';...
            'B-3H ';'E-1H ';'F-2H ';'E-2H ';'E-4H ';'E-4AH';'D-3H ';'D-3AH';'F-3H ';'E-3AH';'B-4AH';'F-4H ';...
            'B-4BH';'D-4AH';'D-1CH';'C-4AH';'B-4DH';'E-3BH';'E-3CH';'E-2AH';'D-3BH';'B-1AH';'B-1BH';'K-3H '};
        
        % zones: Garn, Ile, Tofte, Tilje
        kalmanOptions.comp = [1 0 0 0 0; % C-4H
            0 0 1 0 0; % B-2H
            0 0 1 1 0; % D-1H
            0 0 1 0 0; % D-2H
            1 0 1 1 1; % B-4H
            1 0 1 1 1; % D-4H
            1 0 1 1 1; % C-1H
            1 0 1 1 1; % E-3H
            0 0 0 0 1; % C-2H
            0 0 0 1 0; % B-1H
            1 0 1 1 1; % C-3H
            1 0 1 1 1; % F-1H
            0 0 1 1 0; % B-3H
            0 0 1 0 0; % E-1H
            1 0 1 1 1; % F-2H
            1 0 1 1 1; % E-2H
            1 0 1 1 0; % E-4H
            1 0 0 0 0; % E-4AH
            1 0 1 1 1; % D-3H
            0 0 0 1 0; % D-3AH
            1 0 1 1 1; % F-3AH
            1 0 0 0 0; % E-3H
            1 0 1 1 1; % B-4AH
            1 0 1 1 1; % F-4H
            0 0 1 0 0; % B-4BH
            1 0 0 0 0; % D-4AH
            0 0 1 0 0; % D-1CH
            1 0 1 1 1; % C-4AH
            0 0 1 0 0; % B-4DH
            1 0 1 1 1; % E-3BH
            0 0 1 0 0; % E-3CH
            0 0 1 0 0; % E-2AH
            0 0 1 0 0; % D-3BH
            1 0 0 0 0; % B-1AH
            0 0 0 1 0; % B-1BH
            0 0 1 0 0  % K-3H
            ];
        
        % reservoir layers
        kalmanOptions.layers =  [3 1 7 7 4];
        
    end
    
    save('inputData.mat', 'kalmanOptions','-append');
    disp('inputData.mat with localization is created.')
    
end

%%%%%%%%%%%%%%%%%
% True Solution 
%%%%%%%%%%%%%%%%%

% compute true solution 
if ~existfile('trueSolution.mat')
    computeTrueSolution
end 

disp('Setting up the inputs to runIES are completed')


%-- removing .S* and .X* files --%
%clean_up; 
%close all;


