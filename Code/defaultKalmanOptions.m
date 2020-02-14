function kalmanOptions=defaultKalmanOptions(kalmanOptions)

% function kalmanOptions=defaultKalmanOptions(kalmanOptions)
%
% if no input argument is given, kalmanOptions is read from inputData.mat.
%
% Set default values of KalmanOptions.
%
% Note that there are some fields that do not have default values and
% must be specified manually to be able to run the programs, these are
% marked NEEDED.    
%  
% -.analyticCovMat         - Flag. If set to 1, use analytic covariance
%                            matrix for measurement noise
%                            (default=1).
% -.beta (iES)             - Scalar. Coeffient in the (squared) data mismatch 
%                            threshold beta^2 * p (p = observation size) 
%                            that determines one of the stopping conditions
% -.buildMatricesOldWay    - Flag. If set to 1, the "old way" of
%                            specifying H and W can be used. 
% -.checkCorrelation       - Flag. If set to 1, only free
%                            parameters with significant
%                            correlation to the measurements are 
%                            updated. 
% -.cleanup                - Text string. Name of a file where a script
%                            used for clean up after running parallel
%                            simulations can be included. See
%                            runParSim.m. 
% -.computePredictions     - Vector. Compute predictions using the analyzed 
%                            parameters at the timesteps defined in the 
%                            vector. 
% -.cutVariables           - Scalor. Experimental variable. See
%                            addModelNoise.m.
% -.deleteOldEnsemble      - Flag. If set to 1, delete old ensemble file
%                            (to avoid filling the hard disc'). 
% -.ensembleSize           - Scalar. Number of ensemble members.
% -.EnSRF                  - Flag. If set to 1, the EnSRF filter is used
%                            instead of the default EnKF filter.
% -.ES_script (iES)		    - Text string. iES algorithm to use. So far 
%                             available options: {'LMEnRML','RLM_average_cost'}	
% -.filterMeasurement      - Flag. If set to 1, the measurement will be
%                            filtered in runEnKF according to their
%                            relative distance to the mean of the
%                            ensemble members results in terms of
%                            std. (default = 0, no filtering)
% -.forgettingFactor:      - Flag. If set, add model noise according to the
%                            aposterior standard deviation of the
%                            ensemble (apostStdDev) (implemented in
%                            addModelNoise.m). 
% -.freeparamLB            - Scalar. Lowerbound on free parameters (see
%                            defaultOptions). 
% -.freeparamUB            - Scalar. Upperbound on free parameters (see
%                            defaultOptions).
% -.GMEnKF                 - Structure. If exists use AGM instead of EnKF.
% -.GMEnKF.adaptive        - Flag. If set to 1 adaptive alpha is
%                            calculated, if set to zero user must 
%                            specify -.GMEnKF.adaptive.alpha
% -.GMEnKF.resample        - Scalar. 0 is ordinary resampling and 1 residual 
%                            resampling (recomended), i.e, use resampling or
%                            not.
% - GMEnKF.h               - Scalar between zero and one. Dampening factor 
%                            of the update step. 
% -.GMEnKF.Nefflim         - Scalar. If >=0.8 no resampling. If between 0.8
%                            and 1 do resampling, higher values leads to 
%                            more often resampling. 
% -.HEnKFalpha             - Matrix. Dampening factor alpha (matrix) used with the
%                            Hierarchical EnKF. 
% -.HEnKFNse               - Scalar. Number of sub-ensembles when using the
%                            Hierarchical EnKF. Existent if runEnKF
%                            should call HEnKF.m instead of EnKF.
% -.historicalData         - Flag. If set to 1, the data is historical
%                            data and therefore no measurement noise is
%                            added. 
% -.Hones                  - Flag. Should be set to 0 or 1. Tells which
%                            formate that is used when specifying
%                            the H matrix. See EnKF.m. (default=1)   
% -.ignoreUninformativMeasurments - Flag. If set to 1, if there is no
%                            variance in the values of some measurements
%                            for the ensemble members, these measurements
%                            are not used in computing Ke and updating
%                            the ensemble. Default is 1. 
% -.indepHorLayers         - Flag. The horizontal layers are generated
%                            independently if this field exist.  
% -.initialEnsmble         - Text string. Filename of the .mat file where
%                            the initial ensemble is placed. The variable
%                            ensemble should be stored in this .mat file,
%                            and each ensemble member should specify the
%                            variables to be estimated. If this field
%                            exists, generateInitialEnsemble will use this
%                            initial ensemble.  NEEDED if no
%                            .meanCorrLength, .staticVarMean,
%                            .staticVarStdDev, and .stdCorrLength.
% -.isDumpOne (iES)         - Flag. (for 'RLM_average_cost' only) if 
%                             isDumpOne = 1, discard 1 model member. 
% -.iterate                - Structure or set to 1. If exist, use iterated
%                            EnKF. If set to 1, default values are used
%                            for the structure fields below.    
% -.iterate.epochs         - Scalar. Maximum allowed epochs used for
%                            training of the neural network.
% -.iterate.goal           - Scalar. Convergence criteria used for 
%                            training of the neural network.
% -.iterate.itFile         - Text string. If a small timestep is used, this
%                            option gives the name of the file
%                            containing the TSTEP keyword. This file
%                            is added to the file NEXTTIME.SCH. 
% -.iterate.maxIt          - Scalar. Maximum allowed iterations at each
%                            assimilation step.
% -.iterate.updMeas        - Scalar. 1: use nural networks, 2: use rerun, 
%                            3: use small timestep.
% -.iterES (iES)           - Flag. When set to 1, run iterative
%                             ensemble smoother.
% -.lambda (iES)           - Positive scalar. initial value of the damping 
%                            coeffient in an iES.
% -.lambda_reduction_factor -Scalar. Reduction factor in case to reduce gamma in an
%                            iES.
% -.lambda_increment_factor -Scalar. Increment factor in case to increase gamma in 
%                            an iES.
% -.maxInnerIter (iES)     - Scalar. Maximum number of inner
%                            loop iterations in an iES
% -.maxIter (iES)          - Scalar. Maximum number of outer
%                            loop iterations in an iES
% -.maxUncertainty         - Vector. Do not add more model noise to the a
%                            static variable when it has passed
%                            this value. The field should be
%                            specified for a vector of
%                            lenght(trueStaticVar). It is used in
%                            addModelNoise. 
% -.meanCorrLength         - Vector with the same length as number of
%                            variables in options.staticVar (one value
%                            for each static variable). Mean of
%                            correlation length used while generating
%                            fields of static variables. NEEDED if no
%                            .initialEnsemble.  
% -.measFunc               - Text string. Name of function used to generate
%                            measurement. See computeTrueSolution.m. 
% -.measTrust              - Scalar. See removeMeasurements.m. (default=3
%                            if .filterMeasurement exists) 
% -.minMeasUnc             - Scalar. Lower limit for the measurement
%                            uncertainty used in calculation of
%                            history match in
%                            'testHistMatchNew.m'. This value is
%                            not used in the actual EnKF
%                            simulation.  
% -.minReduction (iES)     - Scalar. minimum relative change of average data 
%                            mismatch (in percentage) in iES
% -.noisePriorToSimulation - Flag. If the variable is set to 0 we run
%                            simulation -> model noise -> filter,
%                            otherwise (default is 1)
%                            model noise -> simulation -> filter.
% -.numFields              - Scalar. Number of seismic fields (used by LEnKF).
% -.numParProcesses        - Scalar. Number of processing units to use
%                            with the parallelScript. (default=1)
% -.numStaticMeas          - Scalar. Number of measured static variables
%                            (default 0).
% -.parallab               - Flag. Set to 1 if you are running parallel
%                            processes on parallab.
% -.permvar                - Vector. Variance in model noise of
%                            staticVariables (i.e., permeability (at
%                            each step)) (default is 0). 
% -.postProcessEnsemble    - String. If set, use this function after each EnKF
%                            update (e.g. used in relation to facies
%                            updating).
% -.predTrust              - Scalar. See removeMeasurements.m. (default=3
%                            if .filterMeasurement exists) 
% -.preCompMeas            - Flag. If set to 1, the true solution is
%                            computed in advance and the data is saved in
%                            trueSolution.mat. (default=1)   
% -.relative               - Vector containing the indexes of which
%                            measurement types of the measurement matrix that
%                            are measured with relative accuracy.  
%                            The accuracy is given in field
%                            -.relativeWeight.
% -.relativeWeight         - Vector, the same size as .relative. Each
%                            element gives a standard deviation for the
%                            measurement type measured with relative
%                            uncertainty, see -.relative. The formula
%                            used are:  index=kalmanOptions.relative;
%                            W(index,index)=diag(1e-6+(abs(truemeasurement(index)).*...
%                            kalmanOptions.relativeWeight).^2); 
% -.retainStaticVarOnly    - Flag. In an iES one only keeps the static variables 
%                            and free parameters  
% -.runWithPluto           - String. Save results and ensembles at this
%                            location. This string is used when the case is run 
%                            on a ram disk (on Pluto).  
% -.saveEnsemble           - Flag. Do not save the ensemble if this is
%                            set to 0 (default=1, i.e., save the
%                            ensemble).
% -.saveMemory             - Flag. To avoid generation and storage of large
%                            matrices(=1). Default is off (=0) for
%                            backward compatibility. 
% -.staticVarLB            - Vector, the same length as
%                            trueStaticVar. Lowerbound for the
%                            estimated static variables.  
% -.staticVarMean          - Vector, the same length as
%                            trueStaticVar. Mean used to generate 
%                            initial ensemble for the static
%                            variables. NEEDED if no .initialEnsemble. 
% -.staticVarStdDev        - Vector, the same length as
%                            trueStaticVar. Standard deviation used to
%                            generate the initial fields of the
%                            static variables. NEEDED if no .initialEnsemble.
% -.staticVarUB            - Vector, the same length as
%                            trueStaticVar. Upperbound for the estimated
%                            static variables. 
% -.stencil                - Vector with three numbers that specifies the area 
%                            that should be used in the localized EnKF (LEnKF).
% -.stdCorrLength          - Vector with the same length as number of
%                            variables in options.staticVar (one value
%                            for each static variable). Standard
%                            deviation in correlation length used to
%                            generate initial ensemble for the static
%                            variables. NEEDED if no .initialEnsemble.  
% -.syntheticGenerated     - Flag. Set to 1 if the .historicalData
%                            options is used, but the data are
%                            generated from synthetic data. Used in
%                            oppsum.m. 
% -.sqrtfilter               Scalar. Existent if runEnKF should call
%                            sqrtFilter.m instead of EnKF.m 
% -.timespec               - Vector. Equal to options.timespec, see
%                            defaultOptions.m. NEEDED.
% -.trueSimName            - Text string. Name of eclipse data file,
%                            <>.DATA, used to generate true data. 
% -.trueScript             - Text string. options.specialScript is the
%                            name of a special script for running
%                            eclipse when generating true data. 
% -.tsvdData (iES)         - Scalar in [0, 1]. threshold for truncated SVD 
%                            that the preserved eigenvalues should sum up to
% -.zeroMeanNoiseStaticVar - Flag. If set to 1, the noise in the static
%                            variables are adjusted such that it always
%                            has zero mean.   
% -.state                  - Scalar. If set, use the value to control the
%                            random seed. Use this option to reproduce a
%                            case.   
% -.truncMeas              - Scalar. Value used for truncating measurements
%                            used in calculation of history match
%                            in 'testHistMatchNew.m'. (default=-1)
% -.useLocalization        - Flag. If set to 1, use localization.
% -.useProd                - Flag. It is possible to run only with seismic
%                            data, and useProd should then be 0. Defalut is 1.
%
% Copyright(c) International Research Institute of Stavanger (IRIS). All Right Reserved.
% $Id: //depot/rfmatlab/main/eclipseKalman/defaultKalmanOptions.m#12 $
% $DateTime: 2017/05/18 15:13:28 $

if nargin==0
  load('inputData','kalmanOptions','options')
end
if nargin==1
  load('inputData','options')
end

% default is to use analytic covariance matrix in the Kalman
% filter.
if ~isfield(kalmanOptions,'analyticCovMat')
  kalmanOptions.analyticCovMat=1;
end

% translate options to internal variables:
if ~isfield(kalmanOptions,'ensembleSize')
    if isfield(kalmanOptions,'initialEnsemble')
        if exist(kalmanOptions.initialEnsemble,'file')
            load(kalmanOptions.initialEnsemble);
            kalmanOptions.ensembleSize = size(ensemble,2);
        end
    else
        kalmanOptions.ensembleSize = 100;
    end
end
if ~isfield(kalmanOptions,'initialEnsemble')
  if ~isfield(kalmanOptions,'staticVarMean')
    error('Missing kalmanOptions.staticVarMean')
  end
  if ~isfield(kalmanOptions,'staticVarStdDev')
    error('Missing kalmanOptions.staticVarStdDev')
  end
  if ~isfield(kalmanOptions,'meanCorrLength')
    error('Missing kalmanOptions.meanCorrLength')
  end
  if ~isfield(kalmanOptions,'stdCorrLength')
    error('Missing kalmanOptions.stdCorrLength')
  end
end 

if ~isfield(kalmanOptions,'timespec')
  error('Missing kalmanOptions.timespec')
end


if isfield(kalmanOptions,'numperm')
  disp('kalmanOptions.numperm  is obsolete, use')
  disp('kalmanOptions.numStaticMeas instead!')
  error('See above!')
end

if ~isfield(kalmanOptions,'permvar') % variance in static variables.
  kalmanOptions.permvar=0;
end
if ~isfield(kalmanOptions,'noisePriorToSimulation')
  kalmanOptions.noisePriorToSimulation=1;
end

% The default is to save the ensemble, check if we are required
% not to do it: (The field is not required.)
if ~isfield(kalmanOptions,'saveEnsemble')
  kalmanOptions.saveEnsemble=1;
end

% if either of the fields kalmanOptions.staticVarUB or
% kalmanOptions.staticVarLB is defined we need both:
if isfield(kalmanOptions,'staticVarLB') && ~isfield(kalmanOptions, ...
						   'staticVarUB') 
  kalmanOptions.staticVarUB=[]; % no upper bound.
end
if ~isfield(kalmanOptions,'staticVarLB') && isfield(kalmanOptions, ...
						   'staticVarUB') 
  kalmanOptions.staticVarLB=[]; % no lower bound.
end

% if either of the fields kalmanOptions.freeparamUB or
% kalmanOptions.freeparamLB is defined we need both:
if isfield(kalmanOptions,'freeparamLB') && ~isfield(kalmanOptions, ...
						   'freeparamUB') 
  if isfield(options,'freeparam')
    if sum(options.freeparam)~=length(kalmanOptions.freeparamLB)
      display(['Inconsistency between sum(options.freeparam) and' ...
	       ' length of freeparam lower/upper bounds']);
      exit
    end
  end
  kalmanOptions.freeparamUB=[]; % no upper bound.
end
if ~isfield(kalmanOptions,'freeparamLB') && isfield(kalmanOptions, ...
						   'freeparamUB') 
  if isfield(options,'freeparam') 
    if sum(options.freeparam)~=length(kalmanOptions.freeparamUB)
      display(['Inconsistency between sum(options.freeparam) and' ...
	       ' length of freeparam lower/upper bounds']);
      exit
    end
  end  
  kalmanOptions.freeparamLB=[]; % no lower bound.
end


if ~isfield(kalmanOptions,'numStaticMeas')
  kalmanOptions.numStaticMeas=0;
end

if ~isfield(kalmanOptions,'preCompMeas')
  kalmanOptions.preCompMeas=1;
end

if ~isfield(kalmanOptions,'filterMeasurement')
  kalmanOptions.filterMeasurement=0;
end
if kalmanOptions.filterMeasurement==1 
  if ~isfield(kalmanOptions,'predTrust')
    kalmanOptions.predTrust=3;
  end
  if ~isfield(kalmanOptions,'measTrust')
    kalmanOptions.measTrust=3;
  end
end

% to restrict the use of memory the following option is included:
if ~isfield(kalmanOptions,'saveMemory')
  kalmanOptions.saveMemory=0;
end

if ~isfield(kalmanOptions,'numParProcesses')
  kalmanOptions.numParProcesses=1;
end

if ~isfield(kalmanOptions,'Hones')
  kalmanOptions.Hones=0;
end

if ~isfield(kalmanOptions,'truncMeas')
  kalmanOptions.truncMeas=-1;
end

if ~isfield(kalmanOptions,'ignoreUninformativeMeasurements')
  kalmanOptions.ignoreUninformativeMeasurements=1;
end

% are we running on parallab?
if ~isfield(kalmanOptions,'parallab')
  if isunix
    [~,w]=unix('uname -a|grep -c fimm');
    if str2double(w)==1
      warning('rfmatlab:eclipseKalman',...
      'I think you are running on parallab, therfore we set')
      disp('kalmanOptions.parallab=1')
      kalmanOptions.parallab=1;
    end
  end
end
if ~isfield(kalmanOptions,'parallab')
    if kalmanOptions.saveMemory==1
        kalmanOptions.parallab=0;
    end
end

if ~isfield(kalmanOptions,'useProd')
    kalmanOptions.useProd = 1;
end

% variables needed if use distance-based localization
% if isfield(kalmanOptions,'distanceLoc')
%     % can only have one variable in the file containing distance-based
%     % localization
%     test=load([kalmanOptions.distanceLoc,'.mat']);
%     test=fieldnames(test);
%     load([kalmanOptions.distanceLoc,'.mat']);
%     kalmanOptions.alphaDB=eval(test{1});  
%     clear test;
%     index=variableIndex(options);
%     kalmanOptions.variable = index.name;
%     kalmanOptions.index = [index.first; index.last]';
%     kalmanOptions.staticVar=options.staticVar;
% end

% save('inputData','-append','kalmanOptions')

