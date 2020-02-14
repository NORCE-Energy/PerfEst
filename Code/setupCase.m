function setupCase(model, varargin) 
useProd = setProperty(varargin,'useProd', true);
useSeismic = setProperty(varargin,'useSeismic', false);
useOrigEns = setProperty(varargin,'useOrigEns', false);
useExisting = setProperty(varargin,'useExisting', false);
useLoc = setProperty(varargin,'useLoc', false);
ensSize = setProperty(varargin,'ensSize', 100);

load('pathInfo.mat');

trueDir = 'SimTrue';

options.modelName = model;
options.paths = pathInfo;

modelCheckout(options, trueDir, true);
load('modelInfo.mat');

% Get ensemble generation function
oldDir = cd(modelDir);
eval(sprintf('ensGenFunc = @%s;', modelInfo.ensGenFunc));
cd(oldDir);

options.eclVer = '2016.1';

options.modelInfo = modelInfo;
options.dim = modelInfo.dim;
options.filename = modelInfo.filename;
options.startDate = modelInfo.startDate;
options.initFile = modelInfo.initFile;
options.paths.modelDir = modelDir;

if isfield(modelInfo, 'actnum')
    options.actnum = modelInfo.actnum;
else
    options.actnum = readEclipseField(...
        fullfile(modelDir, modelInfo.actnumFile), 'ACTNUM');
end

if isfield(modelInfo, 'freeparam')
    options.freeparam = modelInfo.freeparam;
    options.freeparamFunc = modelInfo.freeparamFunc;
    options.freeparamFileIn = modelInfo.freeparamFileIn;
    options.freeparamFileOut = modelInfo.freeparamFileOut;

    kalmanOptions.freeparamLB = modelInfo.freeparamLB;
    kalmanOptions.freeparamUB = modelInfo.freeparamUB;
end

options.statevar= [];
options.staticVar = char(modelInfo.staticVar);

if useOrigEns
    % Tried 20 (too easy)
    % Tried 1  (unstable)
    options.trueRealization = 20;
else
    sdevScale = 1;
    options.ensParams = {
    %   Name     Sdev  Range  WSdevMult
        'PORO',   sdevScale,  4,    0.1;
        'PERMX',  sdevScale,  4,    0.2;
        'PERMY',  sdevScale,  4,    0.2;
        'PERMZ',  sdevScale,  4,    0.1;
    };

    options.genEnsRangeRatio = 0.7;
    options.genEnsSdevRange = 5;
    options.genEnsKriging = true;
end

options.numGridBlocks = sum(options.actnum);
options.fieldSize = prod(options.dim);
options.binary = 2;
options.timespec = generateTimespec(25, 1);
options.historyFile = 'HISTORY.INC';
options.existTimeSteps = 0;
options.returnStatic = 0;

options = defaultOptions(options);

useExistingTruth = useExisting &&...
    exist('trueVar.mat', 'file');

useExistingEns = useExistingTruth &&...
    exist('initEns.mat', 'file');

if ~useExistingEns
    numEns = ensSize;

    [ ensemble, trueStaticVar ] = ensGenFunc(numEns, options);

    if ~useExistingTruth
        save('trueVar.mat', 'trueStaticVar');
    else
        load('trueVar.mat', 'trueStaticVar');
    end
else
    load('initEns.mat', 'ensemble')
    load('trueVar.mat', 'trueStaticVar');

    % Truncate input ensemble to requested size
    ensemble = ensemble(:, 1:ensSize);
end

kalmanOptions.numParProcesses = 6;

kalmanOptions.ensembleSize = size(ensemble, 2);
kalmanOptions.historicalData=0;
dim = ones(options.numGridBlocks,1);
kalmanOptions.staticVarLB = kron(modelInfo.staticVarLB, dim);
kalmanOptions.staticVarUB = kron(modelInfo.staticVarUB, dim);
kalmanOptions.ignoreUninformativeMeasurements=1;
kalmanOptions.timespec=options.timespec(1:40);
kalmanOptions.state = 3.0427e+05;
kalmanOptions.Hones = 1;
kalmanOptions.trueSimName = modelInfo.filename;
kalmanOptions.trueScript = './trueScript';

kalmanOptions.useProd = useProd;
kalmanOptions.useSeismic = useSeismic;

if (~kalmanOptions.useProd) && (~kalmanOptions.useSeismic);
   error('At least we have to use either production or seismic data, if not both ...');
end

seismicOptions = [];

if kalmanOptions.useSeismic
    options.seismicScript = 'generateSeismicData';
    options.matlabRuntime = '/opt/MATLAB';
    %options.seismicRunScript = [ 'run_' options.seismicScript '.sh' ];

    if isfield(options, 'seismicRunScript')
        copyfile(fullfile(options.paths.scriptDir, options.seismicScript), options.seismicScript);
        copyfile(fullfile(options.paths.scriptDir, options.seismicRunScript), options.seismicRunScript);
    end

    seismicOptions = modelInfo.seismicOptions;

    seismicOptions.dataTypes = { 'top_surf', 'b_bottom_surf', 'c_res_tt' };
    seismicOptions.seismic_time = [ 1 10 35 ];
    seismicOptions.snr = 4;

    load(fullfile(modelDir, seismicOptions.seisAttrFile), 'DZ1', 'Top');
    seismicOptions.Top = Top;
    seismicOptions.DZ = DZ1;

    seismicOptions = defaultSeismicOptions(seismicOptions);
end

if useLoc
    kalmanOptions.useLocalization = 1;
    kalmanOptions.denoisingLoc = 1;
    kalmanOptions.tm = 4;
    kalmanOptions.freeparam_threshold = 0.2;
end

% For iES
kalmanOptions.ignoreUninformativeMeasurements = 0;
kalmanOptions.iterES = 1;

%kalmanOptions.ES_script = 'AGSMDA';
%kalmanOptions.ES_script = 'RLM_average_cost';
kalmanOptions.ES_script = 'localMDA';

if strcmp(kalmanOptions.ES_script, 'RLM_average_cost')

    kalmanOptions.append_mean = 1;

    kalmanOptions.maxIter = 6; % max number of outer loop iteration
    kalmanOptions.maxInnerIter = 5; % max number of inner loop iteration
    kalmanOptions.lambda = 1; % initial lambda value
    kalmanOptions.minReduction = 1e-2; % minimum relative change of average data mismatch (in percentage), so 1e-2 actually means .0001
    kalmanOptions.retainStaticVarOnly = 0; % only keep the static variables and free parameters

    kalmanOptions.lambda_reduction_factor = 0.9; % reduction factor in case to reduce gamma
    kalmanOptions.lambda_increment_factor = 2; % increment factor in case to increase gamma 

    kalmanOptions.beta = 1; % coeff. in the (squared) data mismatch threshold beta^2 * p (p = observation size) 

elseif strcmp(kalmanOptions.ES_script, 'AGSMDA')

    kalmanOptions.numIter = 4;
    kalmanOptions.AGSMDA.h = repmat(0.3, 1, kalmanOptions.numIter);
    kalmanOptions.AGSMDA.lambda = repmat(kalmanOptions.numIter, 1, kalmanOptions.numIter);
    kalmanOptions.AGSMDA.resample = 0;
    kalmanOptions.AGSMDA.adaptive = 1;
    kalmanOptions.AGSMDA.Nefflim1 = 0;
    kalmanOptions.AGSMDA.Nefflim2 = 1;

elseif strcmp(kalmanOptions.ES_script, 'localMDA')

    kalmanOptions.numIter = 4;
    kalmanOptions.lambda = repmat(kalmanOptions.numIter, 1, kalmanOptions.numIter);
    kalmanOptions.localRadius = 16;
    kalmanOptions.LMapFile = 'locationMap.mat';
    kalmanOptions.localTaperFunc = 'sine';
    kalmanOptions.localAnalysis = true

end

kalmanOptions.tsvdData = 0.9; % "energy" threshold that the preserved eigenvalues should sum up to

kalmanOptions.initialEnsemble = 'initEns.mat';
kalmanOptions.measFunc = 'generateH';
kalmanOptions.keywords = { 'WBHP', '-[PI]-', 50; 'WOPR', '-P-', 100; 'WWCT', '-P-' 0.02 };

kalmanOptions = defaultKalmanOptions(kalmanOptions);

ensemble = respectBounds(ensemble, kalmanOptions, options);
trueStaticVar = respectBounds(trueStaticVar, kalmanOptions, options);

save(kalmanOptions.initialEnsemble, 'ensemble');

numIter=length(kalmanOptions.timespec);

actnum = options.actnum;
save('inputData','numIter','kalmanOptions','options','seismicOptions','actnum',...
    'trueStaticVar');

% Generate true history
generateTrueHistory(trueStaticVar, options, {
    '^''BR-P-9''',  50;
    '^''BR-P-10''', 50;
    '^''BR-P-11''', 50;
    '^''BR-P-12''', 50;
    '^''BR-P-13''', 50;
    '^''BR-I-9''',  50;
    '^''BR-I-10''', 50});

% Make true run dir
copyfile('inputData.mat', fullfile(trueDir, 'inputData.mat'));
copyfile('HISTORY_TRUE.INC', fullfile(trueDir, options.historyFile));

useExistingProd = useExisting &&...
    exist('trueSolution.mat', 'file') &&...
    exist('trueSolutionSmoother.mat', 'file');

useExistingSeis = useExisting && exist('trueSeismicData.mat', 'file');

if ~useExistingProd && exist('trueSolution.mat', 'file')
    delete('trueSolution.mat')
end

if ~useExistingProd && exist('trueSolutionSmoother.mat', 'file')
    delete('trueSolutionSmoother.mat')
end

if ~useExistingSeis && exist('trueSeismicData.mat', 'file')
    delete('trueSeismicData.mat')
end

% Copy seismic scripts to true dir
if useSeismic && isfield(options, 'seismicRunScript')
    copyfile(options.seismicScript, fullfile(trueDir, options.seismicScript));
    copyfile(options.seismicRunScript, fullfile(trueDir, options.seismicRunScript));
end

root = cd(trueDir);
sfid = fopen(kalmanOptions.trueScript, 'w');
%fprintf(sfid, '#!/bin/sh\n');
%fprintf(sfid, 'unset LD_LIBRARY_PATH\n');
%fprintf(sfid, 'unset XFILESEARCHPATH\n');
%fprintf(sfid, 'unset OSG_LD_LIBRARY_PATH\n');
%fprintf(sfid, 'flow use_TUNING=true %s > ecl.out 2>&1\n', kalmanOptions.trueSimName);
fprintf(sfid, '@eclipse -ver %s %s > ecl.out 2>&1\n', options.eclVer, kalmanOptions.trueSimName);
fclose(sfid);

fileattrib(kalmanOptions.trueScript, '+x');

if ~useExistingProd || ~useExistingSeis
    computeTrueSolution();
    generateLocationMap();
end

cd(root);

if ~useExistingProd
    copyfile(fullfile(trueDir, 'trueSolution.mat'), 'trueSolution.mat');
    copyfile(fullfile(trueDir, 'trueSolutionSmoother.mat'), 'trueSolutionSmoother.mat');
    copyfile(fullfile(trueDir, 'forwSim.mat'), 'forwSim.mat');
    copyfile(fullfile(trueDir, 'reportNum.mat'), 'reportNum.mat');
    copyfile(fullfile(trueDir, 'locationMap.mat'), 'locationMap.mat');
    copyfile(fullfile(trueDir, 'NEXTTIME.SCH'), 'NEXTTIME.SCH');
    copyfile(fullfile(trueDir, sprintf('%s.SMSPEC', modelInfo.filename)),...
        sprintf('%s.SMSPEC', modelInfo.filename));
else
    % Adapt exsisting solution
    tss = load('trueSolutionSmoother.mat');
    lm = load('locationMap.mat');

    if useProd && tss.num_prod == 0
        error('Existing truth does not have production measurements');
    end

    if useSeismic && tss.num_seis == 0
        error('Existing truth does not have seismic measurements');
    end

    if ~useProd
        display('Chopping production data');

        tss.W = tss.W(tss.num_prod + 1:end);
        tss.measurement = tss.measurement(tss.num_prod + 1:end);
        tss.truemeasurement = tss.truemeasurement(tss.num_prod + 1:end);

        lm.LMap = chopLMap(lm.LMap, tss.num_prod + 1, []);
        tss.num_prod = 0;
    end

    if ~useSeismic
        display('Chopping seismic data');
        tss.W = tss.W(1:end-tss.num_seis);
        tss.measurement = tss.measurement(1:end-tss.num_seis);
        tss.truemeasurement = tss.truemeasurement(1:end-tss.num_seis);

        lm.LMap = chopLMap(lm.LMap, [], tss.num_prod);
        tss.num_seis = 0;
    end

    save('trueSolutionSmoother.mat', '-struct', 'tss');
    save('locationMap.mat', '-struct', 'lm');
end

if useSeismic && ~useExistingSeis
    copyfile(fullfile(trueDir, 'trueSeismicData.mat'), 'trueSeismicData.mat');
end

if ~useExistingProd
    display('Write history file');
    writeHistory(fullfile(trueDir, kalmanOptions.trueSimName), options,...
        fullfile(modelDir, modelInfo.compdatTemplate), 'saveConnectionStates', true);
end

setupMultiCore(options, kalmanOptions);

end
