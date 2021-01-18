function runIES(varargin)
% function runIES(varargin)
%
% Prototype for running iterative Ensemble smoother (IES) with Eclipse.
%
% Input give through file 'inputData.mat':
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Contents of options and kalmanOptions:
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Fields of options (some are optional):
% see defaultOptions for description.
%
% Fields of kalmanOptions (some are optional):
% see defaultKalmanOptions for description.
%
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Input:
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% varargin (integer): variable input to specify the start/restart iteration step
% null input implies to (re-)start from step 0
% example usage: runIES; runIES(); runIES(5)
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Output:
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% At each iteration the following files are written (also see the report
% "Levenberg-Marquardt forms of the iterative ensemble Smoother(LM-EnRML).
% Description of Implementation in IRIS Filter Toolbox and theory" and the
% complementary note "An update note on the implementations of Levenberg-
% Marquardt based iterative ensemble smoothers (iES) in IRIS filter toolbox
% "):
%
% 1. “ensemble#.mat”: files containing updated ensemble of static variables
%    at outer-iteration step #. Specifically, “ensemble0.mat” is the initial
%    ensemble. Variables saved in “ensemble#.mat” include "ensemble"
%    (containing static variables estimated at outer-iteration step #) and
%    "lambda" (the regularization parameter).
%
% 2. “objRealIter#.mat”: files containing data mismatch of each realization
%    at outer-iteration step #. Variable saved in “objRealIter#.mat” is
%    "objReal" (ensemble of data mismatch).
%
% 3. “objRealIter#1-#2.mat”: files containing data mismatch of each
%    realization at the inner-iteration step #2 of the outer-iteration #1.
%    These files are saved only when inner-loop iterations take place
%    (e.g., when average data mismatch is not reduced so that one re-excutes
%    an iteration with a smaller step size). Variable saved in
%    “objRealIter#1-#2.mat” is "objReal" (ensemble of data mismatch).
%
% 4. “simulatedDataIter#.mat”: files containing simulated data at
%    outer-iteration step #. Variables saved in “simulatedDataIter#.mat”
%    include "simData" (simulated data normalized by the STDs of the data)
%    and "numberOfFailedSimulations" (number of failed simulations for the
%    ensemble of reservoir models).
%
% 5. “simulatedDataIter#1-#2.mat” contains simulated data at the inner
%    iteration step #2 of the outer-iteration step #1. These files are saved
%    only when inner-loop iterations take place. Variables saved in
%    “simulatedDataIter#1-#2.mat” include "simData" (simulated data) and
%    "numberOfFailedSimulations" (number of failed simulations for the
%    ensemble of reservoir models).
%
% 6. “simulatedDataDiff-#.mat” (for RLM-MAC only): files containing the
%    differences between the mean of the simulated observations and the
%    simulated observation of the mean of the ensemble reservoir models at
%    outer-iteration step #. These files are used to check the potential
%    different behaviour between aLMEnRML and RLM-MAC. Variables saved in
%    “simulatedDataDiff-#.mat” include "simData" (normalized simulated data
%    of the ensemble reservoir models), "simMean" (normalized simulated data
%    of the mean of the ensemble reservoir models), "deltaD" (deltaD =
%    simData - simMean*ones(1,ensemble_size), "deltaD" ( deltaD = simData -
%    simMean*ones(1,ensemble_size), i.e., the ensemble of deviations of
%    simData to simMean), "alt_deltaD" (deltaD in aLMEnRML, deltaD = simData -
%    <simData>*ones(1,ensemble_size), where <simData> is the ensemble mean
%    of simData), and "obsDiff" (obsDiff = (simMean - <simData>) * STD, i.e.,
%    un-normalized difference between simMean and <simData>).
%
% 7. “debugLMEnRML” (for aLMEnRML) and “debugRLMMAC” (RLM-MAC): ASCII log
%    files containing certain information throughout the iteration, including
%    change of the objective function with iteration, regularization
%    parameter value, number of production data used at each iteration
%    when kalmanOptions.ignoreUninformativeMeasurements is set to 1, number
%    of singular values retained at each iteration, etc.
%
% 8. “debugDisLoc”: ASCII log file generated if covariance localization is
%    conducted. It contains information about the column indices of the
%    Kalman gain that correspond to the observations (well production data)
%    to be localized, the observation sites (in terms of well names), and
%    the localization functions used to taper the Kalman gain.
%
% Copyright (c) 2010-2014 IRIS, All Rights Reserved.
% $Id: //depot/rfmatlab/main/eclipseKalman/runIES.m#26 $
% $DateTime: 2019/01/28 13:50:42 $

% set default random number generator
rng('default')
% restart from a specified iteration step
if isempty(varargin)
    iter = 0;
else
    iter = varargin{1};
end

% load the input:
load('inputData','trueStaticVar','kalmanOptions','Q','options');

if ~exist('Q','var')
    disp('Variable Q is not specified, therefor we are not adding model noise')
    Q=[];
end

% make sure that all default values of kalmanOptions is set:
kalmanOptions=defaultKalmanOptions(kalmanOptions) %#ok<*NODEF,*NOPRT>
if ~isfield(kalmanOptions,'pickObs')
    kalmanOptions.pickObs = 1:length(kalmanOptions.timespec);
end
save('inputData','-append','kalmanOptions')

% make sure that all default values of options is set:
options=defaultOptions(options)

% We use the cholesky factorization of Q:
% (Also assuming that Q does only contains the state variables.)
cQt=chol(Q)';

% check if run is on ram disk (pluto)
if isfield(kalmanOptions,'runWithPluto')
    dir = kalmanOptions.runWithPluto;
else
    dir = './';
end

% compute the true solution - if it is not computed already:
if (isfield(kalmanOptions,'preCompMeas') && ...
        kalmanOptions.preCompMeas==1)
    if ~existfile('./trueSolution.mat') && (~isfield(kalmanOptions,'useProd') || kalmanOptions.useProd == 1)
        computeTrueSolution;
    end
else
    truefilecontents=[];
end

%if isfield(kalmanOptions,'distanceLoc')
%    checkForDisLoc;
%end

% generate or read initial ensemble:
if iter == 0
    ensemble=generateInitialEnsemble(kalmanOptions,options);
    if ~isfield(kalmanOptions,'initialEnsemble')
        save('initial_ensemble','ensemble');
    end
    if isfield(kalmanOptions,'append_mean') && kalmanOptions.append_mean
        ensemble = [ensemble, mean(ensemble,2)];
    end
    save('ensemble0.mat','ensemble');
    weights=1/size(ensemble,2).*ones(1,size(ensemble,2));
else
    load([dir,'ensemble',num2str(iter)]); % must contain "ensemble"
end

if isfield(options,'biModel') && (~isempty(options.biModel)) && isfield(kalmanOptions,'distanceLoc') % for facies studies, if additional channel statistics are used as obs
    save('trueSolutionSmoother','-append','real_obs_index','obs_type_index');
end

filecontents=[];
options.existTimeSteps=0;
load trueSolutionSmoother.mat;
if iter  == 0 && ~exist('simulatedDataIter0.mat','file')
    numberOfFailedSimulations = [];
    if kalmanOptions.numParProcesses<2
        [simulatedEnsemble,filecontentsOut,numberOfFailedSimulations]= runForwardSimulations(ensemble,0,...
            kalmanOptions,options,cQt,filecontents); %#ok<*ASGLU>
        jtmp=size(simulatedEnsemble,1)*size(simulatedEnsemble,2); % number of total data
        itmp=size(simulatedEnsemble,3); % number of realizationspw
        simulatedEnsemble=reshape(simulatedEnsemble,jtmp,itmp);
        simData=simulatedEnsemble(H(:,2),:);
    else
        [simData,simulatedEnsemble] = runMultiCoreSim(0,ensemble,H,kalmanOptions,options);
    end
    save(strcat(dir,'simulatedDataIter',num2str(iter)),'simData','simulatedEnsemble',...
        'numberOfFailedSimulations','filecontentsOut');
else
    load([dir,'simulatedDataIter',num2str(iter)],'simData'); % must contain "simData"
end

if ~exist('reportNum.mat','file') || ~exist('trueSolution.mat','file') || ~exist('forwSim.mat','file')
    %error('Recompute true solution');
end

% normalize
if getOption(kalmanOptions,'adjustMeasurementNoise',false)
    W = adjustMeasurementNoise;
end
Wbase=W;
W=ones(size(Wbase, 1), 1); % vector with unit weights
scale = [];
if ~isnan(getOption(kalmanOptions,'scaling'))
    scaling = kalmanOptions.scaling;
    if size(scaling,1) == 1
        scaling = scaling';
    end
    if length(scaling) == 2
        scale = [scaling(1)*ones(num_prod,1);scaling(2)*ones(num_seis,1)];
    elseif length(scaling) == length(measurement)
        scale = scaling;
    else
         error('scaling has wrong dimension');
    end
    measurement=measurement.*scale;
end
measurement= normalizeData(measurement, sqrtW(Wbase));

if strcmp(kalmanOptions.ES_script,'AGSMDA')
    disp('Running AGSMDA')
    AGSMDA(iter,ensemble,simData,W,Wbase,H,measurement,scale,options,kalmanOptions,dir,obsLocation,obsType,weights);
elseif strcmp(kalmanOptions.ES_script,'RLM_MAC')
    disp('Running RLM_MAC') 
    RLM_MAC(iter,ensemble,simData,W,Wbase,H,measurement,scale,options,kalmanOptions,dir,obsLocation,obsType);
elseif strcmp(kalmanOptions.ES_script,'localMDA')
    disp('Running localMDA') 
    localMDA(iter,ensemble,simData,W,Wbase,H,measurement,scale,options,kalmanOptions,dir);
else
    error('Smoother type not defined!');
end

