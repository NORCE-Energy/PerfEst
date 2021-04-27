function initialize

if ~existfile('Cfinered.mat')
    makeTrimmedData
    disp('load trimmed data')


    load('trimmedData','Cfine','perfObs','prm')




    totConc=sum(abs(Cfine),4);
    [numIremoved,numIendRemoved,numJremoved,numJendRemoved,numKremoved,numKendRemoved,comprTotConc]=removeUnactiveLayers(totConc);
    origGrid=size(Cfine);
    Cfine=Cfine(numIremoved+1:end-numIendRemoved,numJremoved+1:end-numJendRemoved,numKremoved+1:end-numKendRemoved,:);
    perfObs=perfObs(numIremoved+1:end-numIendRemoved,numJremoved+1:end-numJendRemoved,numKremoved+1:end-numKendRemoved,:);
    redGrid=size(Cfine);

    save('Cfinered','-v7.3','Cfine','perfObs','numIremoved','numIendRemoved','numJremoved','numJendRemoved','numKremoved','numKendRemoved','redGrid','origGrid')
else
    load('Cfinered')
    load('trimmedData','prm')
    disp('load tree')
end
% $$$ load('simfullbrainindicator-634x515x1_filament.mat'); %#ok<*LOAD>
% $$$ load('simfullbraindisc-634x515x1-tree_filament.mat','data');
% $$$ 
% $$$ Cart = results.Cmat.arterial.im;
[nxF,nyF,nzF,nt] = size(Cfine); %#ok<*ASGLU>
% $$$ Cven = results.Cmat.venous.im;
% $$$ Cfine = Cart + Cven; 
% $$$ perfObs = results.perf;
[~,argMaxFine] = max(Cfine,[],4);
options.fineMask = NaN(nxF,nyF,nzF);
options.fineMask(argMaxFine>1) = 1;
% $$$ Cfine = options.fineMask.*Cfine;
% $$$ perfObs = options.fineMask.*perfObs;

% boundary conditions   FIX/CHECK ( the values found from prm are not coinciding with the paper's values)
prsArt = prm.tree.arterial.bndpress;  % kPa % 10.6; % mmHg
prsVen = prm.tree.venous.bndpress; % kPa  %1.6; % mmHg

%Domain - scale for part cutted away
redFakt=redGrid(1:3)./origGrid(1:3);
xL = redFakt(1)*prm.fov(1)/1000; 
yL = redFakt(2)*prm.fov(2)/1000; 
zL = redFakt(3)*prm.fov(3)/1000;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% determine the upscale level
ny = 4 %128;
nx = round(xL/(yL/4))% 158;
nz = round(zL/(yL/4))
options.L = [nx,ny,nz];
nn=nx*ny*nz;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Ccrs = upscale(Cfine, [nx,ny,nz,nt]);
Pcrs = upscale(perfObs, [nx,ny,nz]);
Ccrs(isnan(Ccrs)) = 0;
Pcrs(isnan(Pcrs)) = 0;

[CcrsMax, argMax] = max(Ccrs,[],4);
maxIndicator = (Ccrs(:,:,:,2:nt-1) > Ccrs(:,:,:,3:nt)) .* ...
               (Ccrs(:,:,:,2:nt-1) > Ccrs(:,:,:,1:nt-2));
accIndicator = sum(maxIndicator,4);

% create options for runFlowSolver3D2Q

% time
options.time=0:1:nt;

%Viscosity CHECK (consistent with paper)
options.visc0 = 3e-3;
visc = options.visc0*ones(nx,ny,nz,2); 

% Load prior information CHECK (values refer to Frog)
options.pA = 0.1; % porosity arteries (0.05 in true) (previous value 0.04)
options.pV = 0.1; % porosity venous (0.1 in true) (previous value 0.08)
options.pQ = 1e-4; % porosity capillaries (not given in true)
options.kA = 1e-13; % permeability arteries (1e-13 in true) (previous value 1e-12)
options.kV = 1e-13; % permeability venous (5e-13 in true) (previous value 2e-12)
%kQ = 3e-9; % permeability capillaries (computed from true alpha 1e-3)
options.kQ = 1e-7; % permeability capillaries 
options.pVeins = 0.1; % porosity in veins (same as true)


%%% NB NB linjene under m? korrigeres for det som er fjernet

% BigBrain Vener: dim(346 x 448 x 319)
% % %    211   173    75
% % %    189   184   126
% % %    330   192   133
% % %    176   218   181
termSink.x=([211 189 330 176]-numIremoved)/redGrid(1)*xL;
termSink.y=([173 184 192 218]-numJremoved)/redGrid(2)*yL;
termSink.z=([ 75 126 133 181]-numKremoved)/redGrid(3)*zL;
termSink.p0=prsVen*1000*[1 1 1 1];%133*10*[1 1 1 1];
termSink.perm=1e-6*[1 1 1 1];
    
% BigBrain Arterier:  dim(346 x 448 x 319)
% % %    185   180    98
% % %    171   184   144
termSrc.x=([185 171]-numIremoved)/redGrid(1)*xL;
termSrc.y=([180 184]-numJremoved)/redGrid(2)*yL;
termSrc.z=([ 98 144]-numKremoved)/redGrid(3)*zL;
termSrc.p0=prsArt*1000*[1 1];%133*85*[1 1];
termSrc.perm=1e-6*[1 1];


% Terminal sources:
options.termSink.i=1+floor(nx*termSink.x/xL);
options.termSink.j=1+floor(ny*termSink.y/yL);
options.termSink.k=1+floor(nz*termSink.z/zL);
options.termSink.p0=termSink.p0;
options.termSink.perm=termSink.perm/options.visc0;
clear termSink;
options.termSrc.i=1+floor(nx*termSrc.x/xL);
options.termSrc.j=1+floor(ny*termSrc.y/yL);
options.termSrc.k=1+floor(nz*termSrc.z/zL);
options.termSrc.p0=termSrc.p0;
options.termSrc.perm=termSrc.perm/options.visc0;
clear termSrc;

%%% NB options delen m? kanskje flyttes til senere i koden
load('/home/gen/Data/simfullbraindisc-384x384x256-tree','data')
% the lines below remove some empty parts of the domain:
data.arterial.tree.bw=data.arterial.tree.bw(numIremoved+1:end-numIendRemoved,numJremoved+1:end-numJendRemoved,numKremoved+1:end-numKendRemoved);
data.venous.tree.bw=data.venous.tree.bw(numIremoved+1:end-numIendRemoved,numJremoved+1:end-numJendRemoved,numKremoved+1:end-numKendRemoved);
prm.dim=size(data.venous.tree.bw);
%if ~exist('priorTissue.mat','file')
%    getPrior(prm,data,options);
%end
%load('priorTissue.mat')
permArtCrs=8e-10
permQCrs=8e-10
permVenCrs=8e-10
permArtCrs=2e-9
permQCrs=2e-9
permVenCrs=2e-9
% correction factor 0.4 based on first initial ensemble
% mean(sum(fullmeasurement)./sum(simData))
poroArtCrs= 0.075 * (13/84)
poroQCrs=0.075 * (7/84)
poroVenCrs = 0.075 * (64/84)
% if ~exist('priorFrog.mat','file')
%     getPriorFrog(prm,data,options);
% end
% load('priorFrog.mat')

% find boundary cells -- must be tuned to example at hand
% FIX/CHECK - do not include for now
% $$$ boundaryCells = []; 
% $$$ N = floor(nyF / ny);
% $$$ M = options.fineMask;
% $$$ M(isnan(M)) = 0;
% $$$ for j=1:N:nyF-N+1
% $$$     if sum(M(end,j:j+N-1)>0) > 0 
% $$$             boundaryCells = [ boundaryCells, floor(1+j/N) ]; %#ok<*AGROW>
% $$$     end
% $$$ end
% $$$ options.boundaryCells = boundaryCells;

% inflow profile   CHECK (currnet values are from initilizeBrain.m)
cValBnd = 5;%25e6;
timeStart=0;
cValBndVen=0;%25e6;
timeEnd = nt;
% $$$ 
% $$$ cValBnd = 5;%25e6;
% $$$ timeStart=0;
% $$$ cValBndVen=0;%25e6;
% $$$ timeStartVen=6;
% $$$ timeEnd = nt;



termSrc.p0=prsArt*1000*[1 1];%133*85*[1 1];
termSrc.perm=1e-6*[1 1];

prsValueWestQ1=1000*prsArt;
prsValueEastQ1=1000*prsArt;
prsValueSouthQ1=1000*prsArt;
prsValueNorthQ1=1000*prsArt;
prsValueDownQ1=1000*prsArt;
prsValueUpQ1=1000*prsArt;
prsValueWestQ2=1000*prsVen;
prsValueEastQ2=1000*prsVen;
prsValueSouthQ2=1000*prsVen;
prsValueNorthQ2=1000*prsVen;
prsValueDownQ2=1000*prsVen;
prsValueUpQ2=1000*prsVen;

% no-flow up and down - note that this overides any corresponding prs cnd
options.noFlowWestQ1 = 1;
options.noFlowEastQ1 = 1;
options.noFlowSouthQ1 = 1;
options.noFlowNorthQ1 = 1;
options.noFlowDownQ1 = 1;
options.noFlowUpQ1 = 1;
options.noFlowDownQ2 = 1;
options.noFlowUpQ2 = 1;
options.noFlowEastQ2=1;
options.noFlowWestQ2=1;
options.noFlowSouthQ2 = 1;
options.noFlowNorthQ2 = 1;

options.xL=xL;
options.yL=yL;
options.zL=zL;
options.nx=nx;
options.ny=ny;
options.nz=nz;
options.visc=visc;

options.prsValueWestQ1 = prsValueWestQ1;
options.prsValueEastQ1 = prsValueEastQ1;
options.prsValueSouthQ1 = prsValueSouthQ1;
options.prsValueNorthQ1 = prsValueNorthQ1;
options.prsValueDownQ1 = prsValueDownQ1 ;
options.prsValueUpQ1 = prsValueUpQ1;
options.prsValueWestQ2 = prsValueWestQ2;
options.prsValueEastQ2 = prsValueEastQ2;
options.prsValueSouthQ2 = prsValueSouthQ2;
options.prsValueNorthQ2 = prsValueNorthQ2;
options.prsValueDownQ2 = prsValueDownQ2;
options.prsValueUpQ2 = prsValueUpQ2;
options.cValBnd=cValBnd;
options.timeStart=timeStart;
%options.gammaShape=1.1;
%options.gammaScale=1.5;
scalingaif = 1e6;  % CHECK (as in initializeBrain.m)
options.empiricalAIF = [prm.reporttimeline,prm.reportaifval*scalingaif];
options.timeEnd=timeEnd;

%chop background
options.mask = zeros(nx,ny,nz,2);
options.maskQ = zeros(nx,ny,nz);
if size(argMax,3)==1
    options.mask(:,:,1,1) = (argMax>1);
    options.mask(:,:,1,2) = (argMax>1);
    options.maskQ(:,:,1) = (argMax>1);
else
    options.mask(:,:,:,1) = (argMax>1);
    options.mask(:,:,:,2) = (argMax>1);
    options.maskQ = (argMax>1);
end



delay = 0; % was 13 before CHECK/FIX
tofArtObserved=(argMax-delay) .* options.mask(:,:,1,1);
if size(argMax,3)==1
    actnum = options.mask(:,:,1,1) + options.mask(:,:,1,2);
else
    actnum = options.mask(:,:,:,1) + options.mask(:,:,:,2);
end
actnum = min(actnum,1);
options.actnum = reshape(actnum,nn,1);

% smoother options
options.dim = [nx,ny,nz];
options.statevar=[];
options.dynamicVar=char('PRESSURE','VEL');
%options.staticVar = char('TRANXART','TRANXVEN','TRANYART','TRANYVEN','TRANZART','TRANZVEN','TRANQ','POROART','POROVEN','POROQ');
options.staticVar = char('PERMXART','PERMXVEN','PERMYART','PERMYVEN','PERMZART','PERMZVEN','TRANQ','POROART','POROVEN','POROQ');
options.fieldSize = nn;
options.useLogPerm = 1;
options.returnStatic = 1;
%options.freeparam = [];
options.numGridBlocks = sum(options.actnum);
options.timespec = options.time;
options.nummeas = sum(options.actnum);

% Setup for iES
options.existTimeSteps = 0; % need to set options.existTimeSteps = 0 when using iES

options = defaultOptions(options);

% initial guess (in 3D)
x0=[log(permArtCrs); log(permVenCrs); log(permArtCrs); log(permVenCrs); ...
    log(permArtCrs); log(permVenCrs); ...
    log(permQCrs); poroArtCrs; poroVenCrs;  poroQCrs]

% kalmanOptions.initialEnsemble='initEns'; % generate this 'on-the-fly'
kalmanOptions.ensembleSize = 200;
trueStaticVar = [];
kalmanOptions.historicalData=1;
tmp=repmat(x0,1,nn)';
kalmanOptions.staticVarMean = tmp(:);
kalmanOptions.ignoreUninformativeMeasurements=0;
kalmanOptions.timespec=options.timespec;
kalmanOptions.reportTime = prm.reporttimeline;

kalmanOptions.state = 2.0427e+05; % Used to control random seed.
kalmanOptions.Hones = 1;
kalmanOptions.numParProcesses = 1;
kalmanOptions.saveMemory = 1;
kalmanOptions.thinobs = 10;
obsType = 'concentration';
if strcmp(obsType,'concentration') % CHECK/FIX
    measurement = reshape(Ccrs,size(Ccrs,1)*size(Ccrs,2)*size(Ccrs,3),size(Ccrs,4));
    measurement = measurement * scalingaif;
    fullmeasurement=measurement;
    save('fullMeas','fullmeasurement')
    if getOption(kalmanOptions,'thinobs',0) > 0 
        % thin out to reduce amount of data
        Nr = floor(length(options.time)/kalmanOptions.thinobs);
        if  Nr > 1 && size(measurement,1) > 1
            
            CV=cov(measurement);
            [a,b]=max(diag(CV));
            measInd=b;
            indCand=1:size(CV,1);
            indCand(b)=[];
            while length(measInd) < Nr % we use Nr measurement points
                detM=zeros(length(indCand),1);
                for I=1:length(indCand)
                    detM(I)=det(CV([measInd indCand(I)],[measInd indCand(I)]));
                end
                [a,b]=max(detM);
                measInd=sort([measInd indCand(b)]);
                indCand(b)=[];
            end
            
        else
            measInd = 1:size(measurement,2);
            [a,b]=max(diff(fullmeasurement'));
            measInd=1:max(b)
        end
        measInd=1:150;
        measInd = unique(measInd);
        kalmanOptions.measInd = measInd;
        measurement = measurement(:,measInd); 
    end
    
elseif strcmp(obsType,'tof')
    measurement = tofArtObserved(:);
elseif strcmp(obsType,'tofAndPeak')
    measurement = [CcrsMax(:) * scalingaif,tofArtObserved(:)];
end
kalmanOptions.obsType = obsType;

% Define bounds     CHECK/FIX
dim = ones(options.fieldSize,1);
%kalmanOptions.staticVarStdDev = [1*dim;1*dim;1*dim;1*dim;1*dim;1*dim;1*dim;0.1*dim;0.1*dim;1e-5*dim];

kalmanOptions.meanCorrLength = 0.01; %floor(options.L/2.1333);
kalmanOptions.stdCorrLength = 0.01;

kalmanOptions.meanCorrLength = 1; %floor(options.L/2.1333);
kalmanOptions.stdCorrLength = 0; % floor(options.L/2.1333)/5;

% Compute initial ensemble for porosity and transmissibility
na = ones(options.numGridBlocks,1);
permLB = -32;
permUB = -14;
permQLB = -32;
permQUB = -14;
poroLB = 0.001;
poroUB = 0.999;
poroQLB = 1e-6;
poroQUB = 0.1;
if strcmp(options.staticVar(1,1:4),'TRAN')
    permLB = ceil(permLB - log(options.visc0));
    permUB = ceil(permUB - log(options.visc0));
    permQLB = ceil(permQLB - log(options.visc0));
    permQUB = ceil(permQUB - log(options.visc0));
end
staticVarLB = [permLB*na;permLB*na;permLB*na;permLB*na;permLB*na;permLB*na;permQLB*na;poroLB*na;poroLB*na;poroQLB*na];
staticVarUB = [permUB*na;permUB*na;permUB*na;permUB*na;permUB*na;permUB*na;permQUB*na;poroUB*na;poroUB*na;poroQUB*na];
kalmanOptions.threshold = 0;
B=load('../RunPerf_1/fullMeas.mat');
for I=1:nn
    corrFact=sum(fullmeasurement(I,:))/sum(B.fullmeasurement);
    kalmanOptions.staticVarMean(7*nn+I:nn:end)=kalmanOptions.staticVarMean(7*nn+I:nn:end)*corrFact;
end
C=load('../RunPerf_1/finalState.mat')
kalmanOptions.staticVarMean(1:7*nn)=log(C.options.permQ);
% increase perm for arteries and veins
kalmanOptions.staticVarMean(1:6*nn)=log(C.options.permQ)+3;
% increas perm for arteries further
%kalmanOptions.staticVarMean(1:nn)=log(C.options.permQ)+5;
%kalmanOptions.staticVarMean(2*nn+1:3*nn)=log(C.options.permQ)+5;
%kalmanOptions.staticVarMean(4*nn+1:5*nn)=log(C.options.permQ)+5;
kalmanOptions.staticVarStdDev = 0.1*abs(kalmanOptions.staticVarMean);
kalmanOptions.staticVarStdDev(1:6*nn)=1;
kalmanOptions.staticVarStdDev(6*nn+1:7*nn)=2;
newEns=0;

load('resPostPros','P_ms','uc')
rateQ=P_ms(:)./uc;
rateQ(isnan(rateQ))=eps;
% formula based on simulation:
kalmanOptions.staticVarMean(6*nn+(1:nn))=log(rateQ(:))-15.2528;
for I=1;6; %[1 3 5]
    kalmanOptions.staticVarMean(nn*(I-1)+(1:nn))=log(rateQ(:))-15.2528+3;
end
kalmanOptions.staticVarStdDev(6*nn+1:7*nn)=0.05;

if 0 %existfile('./simulatedDataIter0.mat') && existfile('./resPostPros.mat')
    newEns=1;
    SD=load('simulatedDataIter0','filecontentsOut');
    E0=load('ensemble0','ensemble');
    load('resPostPros','P_ms','uc')
    permQens=E0.ensemble(6*options.numGridBlocks+1:7*options.numGridBlocks,:);
    rateQ=SD.filecontentsOut.rateQ(options.actnum==1,:);
    for I = 1:nn
        if options.actnum(I)==1
            [a,b]=min(abs(rateQ(:)*uc-P_ms(I)));
            kalmanOptions.staticVarMean(I:2*nn:5*nn)=permQens(b)+3;
            kalmanOptions.staticVarMean(6*nn+I)=permQens(b);
        else
            kalmanOptions.staticVarMean(6*nn+I)=permQLB;
        end
    end
    kalmanOptions.staticVarStdDev(6*nn+1:7*nn)=0.05;
    %clear SD E0
end
    

if ~exist('initial_ensemble.mat','file')
    % ensemble = getFrogEnsemble(kalmanOptions,options);
    !rm simulatedDataIter*
    !rm initial_ensemble.mat
    !rm ensemble*.mat
    ensemble = generateInitialEnsemble(kalmanOptions,options);
    ensemble = adjustVariableWithInBounds(ensemble,staticVarLB,staticVarUB);
    save('initial_ensemble','ensemble');
    disp('initial_ensemble.mat is created.')
else
    load('initial_ensemble','ensemble');
end
kalmanOptions.initialEnsemble = 'initial_ensemble';
kalmanOptions.staticVarLB = staticVarLB;
kalmanOptions.staticVarUB = staticVarUB;
meanEnsemble = mean(ensemble,2);
options = setOptions(options,meanEnsemble);
%options.permZ = 2e-10*ones(nx,ny,nz+1,2)/options.visc0; % FIX not used but must be defined

%Pre-process ...
initialState=runFlowSolver3D2QSort(options);

CA = interp1(options.time',initialState.volTracerArt',prm.reporttimeline);
CA = options.porosityArt(:).*CA';
CQ = interp1(options.time',initialState.volTracerCap',prm.reporttimeline);
CQ = options.porosityQ(:).*CQ';
CV = interp1(options.time',initialState.volTracerVen',prm.reporttimeline);
CV = options.porosityVen(:).*CV';
ref = CA + CQ + CV; %#ok<*NASGU>


% Setup for iES
kalmanOptions.iterES = 1;
% iES options can be specified there.
kalmanOptions.ES_script = 'RLM_MAC';
kalmanOptions.useProd = 0;
kalmanOptions.useSeismic = 0;
kalmanOptions.useMRI = 1;

if strfind(kalmanOptions.ES_script,'RLM_MAC')
    
    % for 'RLM_average_cost' only; if isDumpOne = 1, discard 1 model member so that LMEnRML' and 'RLM_average_cost' have equal numbers of forward model runs
    kalmanOptions.isDumpOne = 0;
    kalmanOptions.append_mean = 1; % for RLM-MAC, append ensemble mean to the end of the ensemble, see runIES_multicore.m
    kalmanOptions.maxIter = 30; % max number of outer loop iteration
    kalmanOptions.maxInnerIter = 5; % max number of inner loop iteration
    kalmanOptions.lambda = 1;
    kalmanOptions.minReduction = 10; % minimum relative change of average data mismatch (in percentage), so 1e-2 actually means .0001
    kalmanOptions.retainStaticVarOnly = 1; % only keep the static variables and free parameters
    
    kalmanOptions.lambda_reduction_factor = 0.9; % reduction factor in case to reduce gamma
    kalmanOptions.lambda_increment_factor = 2; % increment factor in case to increase gamma
    measurement = measurement(options.actnum==1,:); % only keep measurements for active cells
    ref = ref(options.actnum==1,:);
    
elseif strfind(kalmanOptions.ES_script,'AGSMDA')
    %-------------------------------------------------------%
    %--Configuration of AGSMDA type iES --%
    %-------------------------------------------------------%
    kalmanOptions.numIter = 4;
    kalmanOptions.AGSMDA.adaptive = 1;
    kalmanOptions.AGSMDA.resample = 1;
    kalmanOptions.AGSMDA.h = [1.0,0.8,0.5,0.2];
    kalmanOptions.AGSMDA.lambda = ones(1,kalmanOptions.numIter)*kalmanOptions.numIter; % initial lambda value
    kalmanOptions.AGSMDA.Nefflim1 = 0;
    kalmanOptions.AGSMDA.Nefflim2 = 1;
    measurement = measurement(options.actnum==1,:); % only keep measurements for active cells
    ref = ref(options.actnum==1,:);
    
elseif strfind(kalmanOptions.ES_script,'localMDA')
    
    kalmanOptions.numIter = 5;
    kalmanOptions.lambda = repmat(kalmanOptions.numIter, 1, kalmanOptions.numIter);
    kalmanOptions.localRadius = 4;
    kalmanOptions.LMapFile = 'locationMap.mat';
    kalmanOptions.localTaperFunc = 'linear';
    kalmanOptions.localAnalysis = true;
    Nr = 1;
    if isfield(kalmanOptions,'reportTime')
        T = kalmanOptions.reportTime;
        if getOption(kalmanOptions,'thinobs',0) > 0
            T = T(kalmanOptions.measInd); 
        end
        Nr = length(T);
    end
    for I = 1:Nr
        seisLocation{I}.minCoord = [ 1 1 1 ]; %#ok<*AGROW,*SAGROW>
        seisLocation{I}.maxCoord = options.dim;
    end
    options.freeparam = [];
    save('../Data/seisLocation.mat', 'seisLocation');
    generateLocationMap(options,kalmanOptions);
   
else
    error(['Unknown script: ', kalmanOptions.ES_script]);
end

kalmanOptions.tsvdData = 0.999; % "energy" threshold that the preserved eigenvalues should sum up to
kalmanOptions.beta = 1e-2; % coeff. in the (squared) data mismatch threshold beta^2 * p (p = observation size)
Q = [];

measurement = measurement(:);
kalmanOptions.measurement = measurement;
nm = length(measurement);
W = (max(measurement)*0.1)^2 * ones(nm,1); %max(1e-12, (1e-1*measurement).^2);
H(:,1) = 1:nm;
H(:,2) = 1:nm;
obsLocation = []; obsType = [];
save('trueSolutionSmoother.mat','H','W','measurement','obsLocation','obsType','tofArtObserved','perfObs','Pcrs');
disp('trueSolutionSmoother.mat is created.')

% Localization
kalmanOptions.useLocalization = 1;
kalmanOptions.useLocalization = 0;
if getOption(kalmanOptions,'useLocalization',false)
    %kalmanOptions.autoAdaLoc = 1;
    kalmanOptions.denoisingLoc = 1;
    if getOption(kalmanOptions,'autoAdaLoc',false)
        kalmanOptions.nstd = 1;
        for I = 1:size(options.staticVar,1)
            kalmanOptions.thresholding_mode{I} = 'soft';
        end
    %kalmanOptions.locFilter = 3;
    %kalmanOptions.tm = 1; % threshold multiplier for property fields such as PERMX, PORO
    elseif getOption(kalmanOptions,'denoisingLoc',false)
        if isfield(options,'freeparam')
            kalmanOptions.freeparam_threshold = 2/sqrt(kalmanOptions.ensembleSize); % (hard-coded) threshold value for free papameters
        end
        kalmanOptions.hybrid_thresholding = 1;
        if kalmanOptions.hybrid_thresholding
            kalmanOptions.hybrid_tm = 2; %#ok<*STRNU>
        end
    end
end

save('inputData','kalmanOptions','options','trueStaticVar','Q');
disp('inputData.mat is created.')

clear Cfine % stored elsewhere, to save time
disp('start saving initialState')
save('initialState.mat','-v7.3');
disp('initialState.mat is created.')
%end


end

