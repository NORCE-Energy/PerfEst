%function initializeBrain
%% set up for contrast data from Erlend's simulation
if ~exist('trimmedData.mat','file')
    disp('start loading data')
    tic
    warning off;
    load('/home/gen/Data/simfullbrainindicator-384x384x256.mat','results','prm');
    load('/home/gen/Data/simfullbrainsolve-384x384x256','disc')
    warning on;
    toc
    disp('data is loaded')
    
    Cart = results.Cmat.arterial.im;
    [nxF,nyF,nzF,nt] = size(Cart); %#ok<*ASGLU>
    Cven = results.Cmat.venous.im;
    % ioVeins = logical((Cart(:,:,7) > 7e-11) + (Cven(:,:,100) > 1e-10));
    % art = data.arterial.tree.bw;
    % art(ioVeins) = 0;
    % art = repmat(art,1,1,nt);
    % ven = data.venous.tree.bw;
    % ven(ioVeins) = 0;
    % ven = repmat(ven,1,1,nt);
    %Cart(art == 1) = Cart(art==1)/0.05;
    %Cven(ven == 1) = Cven(ven == 1)/0.1;
    %Cven(art == 1) = 0;
    Cfine = Cart + Cven;
    perfObs = results.perf;
    [maxCfine,argMaxFine] = max(Cfine,[],4);
    options.fineMask = NaN(nxF,nyF,nzF);
    options.fineMask(argMaxFine>1) = 1;
    for I=1:size(Cfine,4)
        Cfine(:,:,:,I) = options.fineMask.*Cfine(:,:,:,I);
    end
    perfObs = options.fineMask.*perfObs;
%     totConc=sum(abs(Cfine),4);
%     [numIremoved,numIendRemoved,numJremoved,numJendRemoved,numKremoved,numKendRemoved,comprTotConc]=removeUnactiveLayers(totConc);
%     Cfine=Cfine(numIremoved+1:end-numIendRemoved,numJremoved+1:end-numJendRemoved,numKremoved+1:end-numKendRemoved,:);
%     perfObs=perfObs(numIremoved+1:end-numIendRemoved,numJremoved+1:end-numJendRemoved,numKremoved+1:end-numKendRemoved);
    save('trimmedData','-v7.3','Cfine','perfObs','Cart','Cven','prm','disc') 
end
%% starting from trimmedData
if ~exist('contrastData.mat','file')
    load('trimmedData','Cfine','perfObs','prm')
    totConc=sum(abs(Cfine),4);
    [numIremoved,numIendRemoved,numJremoved,numJendRemoved,numKremoved,numKendRemoved,comprTotConc]=removeUnactiveLayers(totConc);
    Cfine=Cfine(numIremoved+1:end-numIendRemoved,numJremoved+1:end-numJendRemoved,numKremoved+1:end-numKendRemoved,:);
    xL = size(Cfine,1)*prm.fov(1)/1000;
    yL = size(Cfine,2)*prm.fov(2)/1000;
    zL = size(Cfine,3)*prm.fov(3)/1000;
    [Ccrs,N_el_Ccrs]=upscale(Cfine,[2 2 1 150]);
    nt=size(Ccrs,4);
    save contrastData Ccrs N_el_Ccrs  nt xL yL zL prm numIremoved numIendRemoved numJremoved numJendRemoved numKremoved numKendRemoved comprTotConc
end
% clear all
% if ~exist('contrastData.mat','file')
%     load('trimmedData','Cfine','perfObs')
%     totConc=sum(abs(Cfine),4);
%     [numIremoved,numIendRemoved,numJremoved,numJendRemoved,numKremoved,numKendRemoved,comprTotConc]=removeUnactiveLayers(totConc);
%     Cfine=Cfine(numIremoved+1:end-numIendRemoved,numJremoved+1:end-numJendRemoved,numKremoved+1:end-numKendRemoved,:);
%     perfObs=perfObs(numIremoved+1:end-numIendRemoved,numJremoved+1:end-numJendRemoved,numKremoved+1:end-numKendRemoved);
%     
%     [nxF,nyF,nzF,nt] = size(Cfine);
%     
%     nx=floor(nxF/nxF);
%     ny=floor(nyF/nyF);
%     nz=floor(nzF/nzF);
%     Ccrs=zeros(nx,ny,nz,nt);
%     Pcrs=zeros(nx,ny,nz);
%     stepZ=round(nzF/nz);stepY=round(nyF/ny);stepX=round(nyX/nx);
%     
%     for k=1:stepZ:nzF
%         for j=1:stepY:nyF
%             for i=1:stepX:nxF
%                 Ccrs(floor(1+i/4),floor(1+j/4),floor(1+k/4),:) = mean(reshape(Cfine(i:i+3,j:j+3,k:k+3,:),[16*4,nt]),'omitnan');
%                 Pcrs(floor(1+i/4),floor(1+j/4),floor(1+k/4)) = mean(reshape(perfObs(i:i+3,j:j+3,k:k+3),[16*4,1]),'omitnan');
%             end
%         end
%     end
%     Ccrs(isnan(Ccrs)) = 0;
%     Pcrs(isnan(Pcrs)) = 0;
%     % clear Cfine;
%     save contrastData Ccrs Pcrs nx ny nz nxF nyF nzF  nt
% end
%% upscaled model is prepared
load contrastData
load trimmedData prm
nt=size(Ccrs,4);
% skip void part

%trim some boundary cells in order to comply with bc:
disp('trimming')
%Ccrs(158,115:117,:)=0;
totConc=sum(abs(Ccrs),4);
[numIremoved,numIendRemoved,numJremoved,numJendRemoved,numKremoved,numKendRemoved,comprTotConc]=removeUnactiveLayers(totConc);
if numJendRemoved>0
    numJendRemoved=numJendRemoved-1; % we keep one additional layer to prperare for inflow through tubes
end
Ccrs=Ccrs(numIremoved+1:end-numIendRemoved,numJremoved+1:end-numJendRemoved,numKremoved+1:end-numKendRemoved,:);
%Pcrs=Pcrs(numIremoved+1:end-numIendRemoved,numJremoved+1:end-numJendRemoved,numKremoved+1:end-numKendRemoved);
nx=size(Ccrs,1);
ny=size(Ccrs,2);
nz=size(Ccrs,3);

[CcrsMax, argMax] = max(Ccrs,[],4);

maxIndicator = (Ccrs(:,:,:,2:nt-1) > Ccrs(:,:,:,3:nt)) .*  (Ccrs(:,:,:,2:nt-1) > Ccrs(:,:,:,1:nt-2));
accIndicator = sum(maxIndicator,4);

%clear options;

% create options for runFlowSolver3D2Q

% time
options.time=0:1:nt;

% %Discretization
% nx = 158;
% ny = 128;
% nz = 1;

%Domain
%load('trimmedData','prm')

nn=nx*ny*nz;

%Viscosity
options.visc0 = 3e-3;
visc = options.visc0*ones(nx,ny,nz,2); 

% Load prior information
options.pA = 0.1; % porosity arteries (0.05 in true)
options.pV = 0.1; % porosity venous (0.1 in true)
options.pQ = 1e-4; % porosity capillaries (not given in true)
options.kA = 1e-12; % permeability arteries (1e-12 in true)
options.kV = 1e-12; % permeability venous (5e-12 in true)
%options.kQ = 1.5e-9; % permeability capillaries (computed from true alpha 5e-4)
options.kQ = 5e-8; % permeability capillaries 
options.pVeins = 0.5; % porosity in veins (same as true)
% if ~exist('priorBrain.mat','file')
%     load('simfullbraindisc-384x384x256-tree.mat','data');
%     getPriorBrain(prm,data,options);
% end
% load('priorBrain.mat')

% inflow profile
cValBnd = 5;%25e6;
timeStart=0;
cValBndVen=0;%25e6;
timeEnd = nt;

% boundary conditions
prsArt = 10.6; % mmHg
prsVen = 1.6; % mmHg
prsValueWestQ1=133*prsArt;
prsValueEastQ1=133*prsArt;
prsValueSouthQ1=133*prsArt;
prsValueNorthQ1=133*prsArt;
prsValueDownQ1=133*prsArt;
prsValueUpQ1=133*prsArt;
prsValueWestQ2=133*prsVen;
prsValueEastQ2=133*prsVen;
prsValueSouthQ2=133*prsVen;
prsValueNorthQ2=133*prsVen;
prsValueDownQ2=133*prsVen;
prsValueUpQ2=133*prsVen;

% no-flow up and down - note that this overides any corresponding prs cnd
options.noFlowWestQ1 = 1;
options.noFlowEastQ1 = 1;
options.noFlowSouthQ1 = 1;
options.noFlowNorthQ1 = 0;
options.noFlowDownQ1 = 1;
options.noFlowUpQ1 = 1;
options.noFlowDownQ2 = 1;
options.noFlowUpQ2 = 1;
options.noFlowEastQ2=1;
options.noFlowWestQ2=1;
options.noFlowSouthQ2 = 1;
options.noFlowNorthQ2 = 0;

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
scalingaif = 1e6;
options.empiricalAIF = [prm.reporttimeline,prm.reportaifval*scalingaif];
options.timeEnd=timeEnd;

%chop background
options.mask = zeros(nx,ny,nz,2);
options.mask(:,:,:,1) = (argMax>1);
options.mask(:,:,:,2) = (argMax>1);
options.maskQ = zeros(nx,ny,nz);
options.maskQ(:,:,:) = (argMax>1);

% separate venes and arteries:
% options.mask(150:158,50:51,1,1) = 0;
% options.mask(150:158,52:53,1,2) = 0;
% options.mask(150:158,79,1,1) = 0;
% options.mask(150:158,77:78,1,2) = 0;
% options.maskQ(150:158,50:53,1) = 0;
% options.maskQ(150:158,77:79,1) = 0;

delay = 0; % was 13 before
tofArtObserved=(argMax-delay) .* options.mask(:,:,:,1);
actnum = options.mask(:,:,:,1) + options.mask(:,:,:,2);
actnum = min(actnum,1);
options.actnum = reshape(actnum,nn,1);

% smoother options
options.dim = [nx,ny,nz];
options.statevar=[];
options.dynamicVar=char('PRESSURE','VEL');
options.staticVar = char('PERMXART','PERMXVEN','PERMYART','PERMYVEN','PERMZART','PERMZVEN','PERMQ','POROART','POROVEN','POROQ');
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

% x0(1:nn) = log(permArtCrs(numIremoved+1:end-numIendRemoved,numJremoved+1:end-numJendRemoved,numKremoved+1:end-numKendRemoved));
% x0(nn+1:2*nn) = log(permVenCrs(numIremoved+1:end-numIendRemoved,numJremoved+1:end-numJendRemoved,numKremoved+1:end-numKendRemoved));
% x0(2*nn+1:3*nn) = log(permArtCrs(numIremoved+1:end-numIendRemoved,numJremoved+1:end-numJendRemoved,numKremoved+1:end-numKendRemoved));
% x0(3*nn+1:4*nn) = log(permVenCrs(numIremoved+1:end-numIendRemoved,numJremoved+1:end-numJendRemoved,numKremoved+1:end-numKendRemoved));
% x0(4*nn+1:5*nn) = log(permArtCrs(numIremoved+1:end-numIendRemoved,numJremoved+1:end-numJendRemoved,numKremoved+1:end-numKendRemoved));
% x0(5*nn+1:6*nn) = log(permVenCrs(numIremoved+1:end-numIendRemoved,numJremoved+1:end-numJendRemoved,numKremoved+1:end-numKendRemoved));
% x0(6*nn+1:7*nn) = log(permQCrs(numIremoved+1:end-numIendRemoved,numJremoved+1:end-numJendRemoved,numKremoved+1:end-numKendRemoved));
% x0(7*nn+1:8*nn) = poroArtCrs(numIremoved+1:end-numIendRemoved,numJremoved+1:end-numJendRemoved,numKremoved+1:end-numKendRemoved);
% x0(8*nn+1:9*nn) = poroVenCrs(numIremoved+1:end-numIendRemoved,numJremoved+1:end-numJendRemoved,numKremoved+1:end-numKendRemoved);
% x0(9*nn+1:10*nn) = poroQCrs(numIremoved+1:end-numIendRemoved,numJremoved+1:end-numJendRemoved,numKremoved+1:end-numKendRemoved);

% kalmanOptions.initialEnsemble='initEns'; % generate this 'on-the-fly'
kalmanOptions.ensembleSize = 200;
trueStaticVar = [];
kalmanOptions.historicalData=1;
%kalmanOptions.staticVarMean = x0';
kalmanOptions.ignoreUninformativeMeasurements=0;
kalmanOptions.timespec=options.timespec;
kalmanOptions.reportTime = prm.reporttimeline;

kalmanOptions.state = 2.0427e+05; % Used to control random seed.
kalmanOptions.Hones = 1;
%kalmanOptions = defaultKalmanOptions(kalmanOptions);
kalmanOptions.numParProcesses = 1;
kalmanOptions.saveMemory = 1;

% find measurement times

%obsType = 'concentration';
obsType='concentration'
if strcmp(obsType,'concentration')
    %findMeasurementTimes
    %   scalingaif = 1e6;
    measurement = reshape(Ccrs*scalingaif,size(Ccrs,1)*size(Ccrs,2)*size(Ccrs,3),size(Ccrs,4));
    A=load('../meanvalsSubjo2.mat')
    measurement=A.meanvalbrainandvessels';
    if size(measurement,1)>1
        CV=cov(measurement);
        [a,b]=max(diag(CV));
        ind=b;
        indCand=1:size(CV,1);
        indCand(b)=[];
        while length(ind)<15 % we use 15 measurement points
            detM=zeros(length(indCand),1);
            for I=1:length(indCand)
                detM(I)=det(CV([ind indCand(I)],[ind indCand(I)]));
            end
            [a,b]=max(detM);
            ind=sort([ind indCand(b)]);
            indCand(b)=[];
        end
        kalmanOptions.thinobs = ind;
        %if getOption(kalmanOptions,'thinobs',0) > 0
        % thin out to reduce amount of data
        measurement = measurement(:,kalmanOptions.thinobs);
    end
    %end
    %measurement = measurement * scalingaif;
elseif strcmp(obsType,'tof')
    measurement = tofArtObserved(:);
elseif strcmp(obsType,'tofAndPeak')
    measurement = [CcrsMax(:) * scalingaif,tofArtObserved(:)];
elseif strcmp(obsType,'handTailed')
    measurement = reshape(Ccrs,size(Ccrs,1)*size(Ccrs,2)*size(Ccrs,3),size(Ccrs,4))* scalingaif;
    [maxmeas,timeOfMaxMeas]=max(measurement,[],2,'omitnan');
    measurement=sum(measurement,2,'omitnan');
    measurement=[measurement maxmeas timeOfMaxMeas];
    noiseLevels=0.1*max(measurement,[],'omitnan');noiseLevels(3)=0.5;noiseLevels(1)=0.7
end
kalmanOptions.obsType = obsType;

% Define bounds
dim = ones(options.fieldSize,1);
kalmanOptions.staticVarStdDev = [1*dim;1*dim;1*dim;1*dim;1*dim;1*dim;1*dim;0.1*dim;0.1*dim;1e-5*dim];
kalmanOptions.meanCorrLength = 60;
kalmanOptions.stdCorrLength = 1;

% Compute initial ensemble for porosity and transmissibility
na = ones(options.numGridBlocks,1);
permLB = -32;
permUB = -6;
permQLB = -32;
permQUB = -6;
poroLB = 0.001;
poroUB = 0.25;%0.999;
poroQLB = 0.001;%1e-6;
poroQUB = 0.25;
if strcmp(options.staticVar(1,1:4),'TRAN')
    permLB = ceil(permLB - log(options.visc0));
    permUB = ceil(permUB - log(options.visc0));
end
staticVarLB = [permLB*na;permLB*na;permLB*na;permLB*na;permLB*na;permLB*na;permQLB*na;poroLB*na;poroLB*na;poroQLB*na];
staticVarUB = [permUB*na;permUB*na;permUB*na;permUB*na;permUB*na;permUB*na;permQUB*na;poroUB*na;poroUB*na;poroQUB*na];

%load('simfullbrainsolve-384x384x256','disc')
if na>1
    load('discmask','discmask') % contains the  disc.mask in discmask
    compressVen=zeros(nx,ny,nz);
    for k=1:4:nzF-3
        for j=1:4:nyF-3
            for i=1:4:nxF-3
                compressVen(floor(1+i/4),floor(1+j/4),floor(1+k/4))=max(reshape(discmask.venous.tree.root.bw(i:i+3,j:j+3,k:k+3),[16*4,1]));
            end
        end
    end
    compressVen=compressVen(numIremoved+1:end-numIendRemoved,numJremoved+1:end-numJendRemoved,numKremoved+1:end-numKendRemoved);
    
    compressArt=zeros(nx,ny,nz);
    for k=1:4:nzF-3
        for j=1:4:nyF-3
            for i=1:4:nxF-3
                compressArt(floor(1+i/4),floor(1+j/4),floor(1+k/4))=max(reshape(discmask.arterial.tree.root.bw(i:i+3,j:j+3,k:k+3),[16*4,1]));
            end
        end
    end
    compressArt=compressArt(numIremoved+1:end-numIendRemoved,numJremoved+1:end-numJendRemoved,numKremoved+1:end-numKendRemoved);
end
save dataForNewFunction kalmanOptions options measurement numIremoved numIendRemoved numJendRemoved numKendRemoved
%compressArt compressVen
%% prepare for tubes to/from arteiral/veneous roots
if na>1
    options=prepareTubes(options,compressArt,compressVen);
end
%% prepare initial ensemble
disp('generating initial ensemble')
tic
if ~exist('initial_ensemble.mat','file')
    if na>1
        ensemble = ensembleFromBrainSimulation(kalmanOptions,options,measurement);
        ensemble = adjustVariableWithInBounds(ensemble,staticVarLB,staticVarUB);
        save('initial_ensemble','ensemble');
        disp('initial_ensemble.mat is created.')
    else
        %         for I=1:200
        %             ensemble(:,I)=[staticVarUB(1:7)+0.5*(staticVarLB(1:7)-staticVarUB(1:7)).*rand(7,1);...
        %                 staticVarLB(8:10)+0.25*(staticVarUB(8:10)-staticVarLB(8:10)).*rand(3,1)];
        %         end
        ensemble=repmat([log(0.3e-7); log(0.6e-7);log(0.3e-7); log(0.6e-7);log(0.3e-7); log(0.6e-7);log(1.5e-6);0.025;0.015;0.01],1,200)+...
            diag([1; 1; 1; 1; 1; 1; 1; 0.003;0.003;0.003])*randn(10,200);
        %        ensemble=repmat([log(0.3e-7); log(0.6e-7);log(0.3e-7); log(0.6e-7);log(0.3e-7); log(0.6e-7);log(1.5e-6);0.05;0.05;0.05],1,200)+...
        %              diag([1; 1; 1; 1; 1; 1; 1; 0.003;0.003;0.003])*randn(10,200);
        ensemble(1,:)=ensemble(1,:)+5;
        for I=2:7
            ensemble(I,:)=ensemble(1,:);
        end
        ensemble = adjustVariableWithInBounds(ensemble,staticVarLB,staticVarUB);
        save('initial_ensemble','ensemble');
        disp('initial_ensemble.mat is created.')
    end
else
    load('initial_ensemble','ensemble');
end
toc
disp('initial ensemble is generated')
kalmanOptions.initialEnsemble = 'initial_ensemble';
kalmanOptions.staticVarLB = staticVarLB;
kalmanOptions.staticVarUB = staticVarUB;

meanEnsemble = mean(ensemble,2);
options = setOptions(options,meanEnsemble);
kalmanOptions.staticVarMean = meanEnsemble;

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
    kalmanOptions.maxIter = 10; % max number of outer loop iteration
    kalmanOptions.maxInnerIter = 5; % max number of inner loop iteration
    kalmanOptions.lambda = 1;
    kalmanOptions.minReduction = 1; % minimum relative change of average data mismatch (in percentage), so 1e-2 actually means .0001
    kalmanOptions.retainStaticVarOnly = 1; % only keep the static variables and free parameters
    
    kalmanOptions.lambda_reduction_factor = 0.9; % reduction factor in case to reduce gamma
    kalmanOptions.lambda_increment_factor = 2; % increment factor in case to increase gamma
    measurement = measurement(options.actnum==1,:); % only keep measurements for active cells
    ref = ref(options.actnum==1,:);
    if na==1
        kalmanOptions.maxIter = 100;
    end
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
            T = T(1:kalmanOptions.thinobs:end); 
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
%W = max(1e-12, (1e-1*measurement).^2);
%W = (0.1^2)*ones(size(measurement)); %max(1e-6, (1e-1*measurement).^2); % data is order 1000 higher than for the frog
W = ((max(measurement(:))/10)^2)*ones(size(measurement)); 
if strcmp(obsType,'handTailed')
    W(1:length(measurement)/3)=noiseLevels(1);
    W(1+length(measurement)/3:2*(length(measurement)/3))=noiseLevels(2);
    W(1+2*(length(measurement)/3):end)=noiseLevels(3);
end
    
H(:,1) = 1:nm;
H(:,2) = 1:nm;
obsLocation = []; obsType = [];
save('trueSolutionSmoother.mat','H','W','measurement','obsLocation','obsType','tofArtObserved')%,'perfObs','Pcrs');
disp('trueSolutionSmoother.mat is created.')


% Localization
kalmanOptions.useLocalization = 1;
if getOption(kalmanOptions,'useLocalization',false)
    kalmanOptions.autoAdaLoc = 1;
    %kalmanOptions.denoisingLoc = 1;
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

load('initial_ensemble');
load('inputData','options')
xValue = mean(ensemble,2);%mean(ensemble,2);
options = setOptions(options,xValue); %#ok<*NODEF>

%Post-process ...
load trimmedData prm
initialState=runFlowSolver3D2QSort(options); %#ok<*NASGU>
CA = interp1(options.time',initialState.volTracerArt',prm.reporttimeline);
CA = options.porosityArt(:).*CA';
CQ = interp1(options.time',initialState.volTracerCap',prm.reporttimeline);
CQ = options.porosityQ(:).*CQ';
CV = interp1(options.time',initialState.volTracerVen',prm.reporttimeline);
CV = options.porosityVen(:).*CV';
ref = CA + CQ + CV;

save('initialState.mat','initialState','options','ref');    
%save('initialState.mat','-v7.3');
%disp('initialState.mat is created.')
%% prepare for run
clear ystat
figure
plot(measurement,'r'),hold on
for I=1:size(ensemble,2)
    Y=bloodFlow(ensemble(:,I),measurement,W,options,prm);
    ystat(I)=norm(Y);
    plot(sqrt(W).*Y+measurement,'b')
    plot(measurement,'r')
end
[a,b]=min(ystat)
%%
if 1 % set to 0 to skip this part
    %     Cstart(1:7,1)=staticVarUB(1:7)-1;
    %     Cstart(3:4)=-15;
    %     Cstart(7)=-15;
    %     Cstart(8:10,1)=0.062*rand(3,1);Cstart(8:10)=Cstart(8:10)*(0.062/sum(Cstart(8:10)));
    Cstart=ensemble(:,b);
    %Cstart=Copt;
    iderord=2;
      %denom=1e-7*ones(size(Cstart));denom(1:7)=10*denom(1:7)
    denom=1e-4*ones(10,1);    
    ConMat=[eye(length(Cstart));-eye(length(Cstart))];
    bound=[staticVarLB;-staticVarUB];
    Dstart=100;
    Dstop=1e-7;%0.1*min(denom)
    D=eye(10);%D(8,8)=1000;D(9,9)=1000;D(10,10)=1000;
    options.maxit=1000;
    prefunc=[];
    func='bloodFlow';
    numprepar=0;
    varargin={measurement,W,options,prm};
    [Copt,Y,varargout] = ...
        LM ( Cstart,1,iderord,denom,ConMat,bound,Dstart,Dstop,...
        D,options,prefunc,func,numprepar,measurement,W,options,prm);
    [J,iderordout,chkder,B]=derivation(Copt,Y,denom,iderord,ConMat,bound,func,measurement,W,options,prm);
    plot(sqrt(W).*Y+measurement,'go')
end
%% run and post process
runSmoothFrog
%perfTest=perfObs;
%perfTest(maxCfine>0)=perfTest(maxCfine>0)+eps;
load perfWithArtVen
Cfine=kalmanOptions.measurement';
Cfine=reshape(Cfine,1,1,1,150);
[P,N]=postProcess([1 1 1],perfTest,Cfine*1e-6,prm,'regions')
