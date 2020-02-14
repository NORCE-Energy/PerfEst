function generateSeismicData(sim_index)

% function generateSeismicData(sim_index)
%
% generate observed AVA traces or Impedance 
%
% Copyright (c) 2010-2016 IRIS, All Rights Reserved.
% $Id: //depot/rfmatlab/main/seismicHM/generateSeismicData.m#6 $
% $DateTime: 2018/12/05 14:14:11 $


% check input
if nargin < 1
    sim_index = 0;
end

% check input : sim_index is string if function is compiled
if ischar(sim_index)
    sim_index = str2double(sim_index);
end

% load input
load('inputData.mat');
active_index = find(options.actnum==1);
inactive_index = find(options.actnum==0);
    
% do nothing if seismic is not set
if ~isfield(kalmanOptions,'useSeismic') || kalmanOptions.useSeismic == 0
    return;
end

% set default random number generator
rng(sim_index);

% initialize 
if sim_index < 1 % true measurements

    cA_leading_index = []; % approximation coeff.
    cA_leading_number = []; % approximation coeff.
    cD_leading_index = []; % detail coeff.
    cD_leading_number = []; % detail coeff.
    FD_WDEC_rec = []; 
    
    measurement = [];
    truemeasurement = [];
    W = [];
    
    dataFile = kalmanOptions.trueSimName;
    
else % simulated measurements
    
    disp(['Generating seismic for realization ',num2str(sim_index)]);

    load('../trueSeismicData.mat','W','cD_leading_index',...
         'cD_leading_number','cA_leading_index','cA_leading_number');
    simSeismicData = [];
    dataFile = sprintf('%s_%d', options.filename, sim_index); %#ok<*NODEF>
    
end

if sim_index >= 1 || (~isfield(kalmanOptions,'historicalData') || ...
                     kalmanOptions.historicalData == 0)
    % find poro and ntg
    filename = [dataFile,'.INIT'];
    if ~exist(['./' filename],'file')
        error([filename ' does not exist']);
    end
    filecontents=readXfile(filename,2);
    if ~isfield(filecontents,'NTG')
        active_NTG = [];
    else
        active_NTG = filecontents.NTG;
    end
    if ~isfield(filecontents,'PORO')
        error('PORO not specified in INIT file');
    else
        active_Poro = filecontents.PORO;
    end
    
    PORO = NaN .* ones(options.fieldSize,1);
    PORO(active_index) = active_Poro;
    PORO = reshape(PORO,options.dim);
    if ~isempty(active_NTG)
        NTG = NaN .* ones(options.fieldSize,1);
        NTG(active_index) = active_NTG;
        NTG = reshape(NTG,options.dim);
        VCLAY = 1 - NTG;
    else
        VCLAY = [];
    end
end

seismic_time = seismicOptions.seismic_time;
data_in_use = seismicOptions.data_in_use;
threshold = seismicOptions.threshold_multiplier;
dataType = seismicOptions.dataType; % AVO or Impedance
colorRange = seismicOptions.colorRange; % colored noise
colored_noise = false;
if ~isempty(colorRange)
    colored_noise = true;
end
showFigure =  seismicOptions.showFigure;

% These are the values for over, under, and side burden
Vp_shale = seismicOptions.Vp_shale;
Vs_shale = seismicOptions.Vs_shale;
Den_shale = seismicOptions.Den_shale;

NT = length(data_in_use);
FD_noisy = cell(1,NT); % FD = 4D
FD_true = cell(1,NT);
FD_rec = cell(1,NT);

options.statevar=char('PRESSURE','SWAT','SGAS','WAT_DEN','OIL_DEN','GAS_DEN','RS');

if exist('preProcessedSeismicData.mat','file') && sim_index == 0
    load('preProcessedSeismicData.mat');
else
    
    for seismic_time_index = 1 : length(seismic_time)
        
        close all;
        
        % compute simulated measurements
        if sim_index >= 1 || (~isfield(kalmanOptions,'historicalData') || ...
                kalmanOptions.historicalData == 0)
            
            filename=strcat(dataFile,sprintf('.X%04i',seismic_time(seismic_time_index))); % make sure there are restart files like *.X0001
            if ~exist(['./' filename],'file')
                error([filename ' does not exist']);
            end
            
            filecontents=readXfile(filename,2);
            
            for i=1:size(options.statevar,1)
                state_name = options.statevar(i,:);
                if ~isfield(filecontents,deblank(state_name))
                    eval(['filecontents.',deblank(state_name),' = zeros(options.numGridBlocks,1);']);
                end
                
                % obtain 'PRESSURE','SWAT','SGAS','OIL_DEN','WAT_DEN','GAS_DEN'
                eval([state_name '= NaN .* ones(options.fieldSize,1);'])
                eval([state_name '(active_index,1) = filecontents.' state_name ';' ]) %
                
                eval([state_name ' = reshape(' state_name ',options.dim);' ]) %
            end
            
            % unit conversion
            if isfield(seismicOptions,'fieldUnits') && (seismicOptions.fieldUnits == 1)
                PRESSURE = PRESSURE * 0.0689476; % Psi to Bar
                WAT_DEN = WAT_DEN * 16.01846; % Pound/ft3 to Kg/m3
                GAS_DEN = GAS_DEN * 16.01846; % Pound/ft3 to Kg/m3
                OIL_DEN = OIL_DEN * 16.01846; % Pound/ft3 to Kg/m3
            end
            
            % ensure physical range
            PRESSURE(~isnan(PRESSURE)) = max(1,PRESSURE(~isnan(PRESSURE)));
            PORO(~isnan(PORO)) = max(1e-3,PORO(~isnan(PORO)));
            PORO(~isnan(PORO)) = min(1,PORO(~isnan(PORO)));
            VCLAY(~isnan(VCLAY)) = max(1e-3,VCLAY(~isnan(VCLAY)));
            VCLAY(~isnan(VCLAY)) = min(1,VCLAY(~isnan(VCLAY)));
            RS(~isnan(RS)) = max(1e-3,RS(~isnan(RS)));
            WAT_DEN(~isnan(WAT_DEN)) = max(1e-3,WAT_DEN(~isnan(WAT_DEN)));
            GAS_DEN(~isnan(GAS_DEN)) = max(1e-3,GAS_DEN(~isnan(GAS_DEN)));
            OIL_DEN(~isnan(OIL_DEN)) = max(1e-3,OIL_DEN(~isnan(OIL_DEN)));
            SWAT(~isnan(SWAT)) = max(0,SWAT(~isnan(SWAT)));
            SWAT(~isnan(SWAT)) = min(1,SWAT(~isnan(SWAT)));
            SGAS(~isnan(SGAS)) = max(0,SGAS(~isnan(SGAS)));
            SGAS(~isnan(SGAS)) = min(1,SGAS(~isnan(SGAS)));
            
            [Vp_sat,Vs_sat,Den_sat,TS] = RockPhysicsModel(SWAT,SGAS,WAT_DEN,OIL_DEN,GAS_DEN,PRESSURE,PORO,VCLAY,RS);
            
            if strcmp(dataType,'AVO')
                [data_true,R0,G, ~,~,Time,Ricker,Tricker] = AVO(Vp_sat,Vs_sat,Den_sat);
                dim = [options.dim(1:2),size(Time,3)];
            elseif strcmp(dataType,'Impedance')
                data_true{1} = Vp_sat.*Den_sat; % p-impedance
                %data_true{1}(inactive_index) = mean(data_true{1}(active_index)); %Vp_shale .* Den_shale;
                data_true{2} = Vs_sat.*Den_sat; % s-impedance
                %data_true{2}(inactive_index) = mean(data_true{2}(active_index)); %Vs_shale .* Den_shale;
                dim = options.dim;
            elseif strcmp(dataType,'Timeshift')
                data_true{1} = TS;
                dim = options.dim;
            end
            
            % stds
            NT = length(data_true);
            std_data = cell(1,NT);
            noise = zeros(1,NT);
            data = cell(1,NT);
            for I = 1:NT
                
                if strcmp(dataType,'AVO')
                    
                    snr = seismicOptions.snr;
                    if snr > 1e-9
                        std_data{I} = std(data_true{I}(:));
                        noise(I) = std_data{I} * sqrt(snr);
                        data{I} = data_true{I} + noise(I) * randn(size(data_true{I}));
                    else
                        data{I} = data_true{I};
                    end
                    
                else
                    
                    noise_std = seismicOptions.noise_std;
                    if noise_std > 1.0e-6
                        if ~colored_noise
                            data{I} = data_true{I} + noise_std * randn(size(data_true{I}));
                        else
                            if length(dim) ~= length(colorRange)
                                error('color range must match the dimention');
                            end
                            currentNoise = fastGaussian3d(dim,noise_std,colorRange);
                            data{I} = data_true{I} + reshape(currentNoise,dim(1),dim(2),dim(3));
                        end
                    else
                        data{I} = data_true{I};
                    end
                    
                end
                
            end
            
        else % collect real data
            load(['Reservoir_gridded_data_',num2str(seismic_time_index),'.mat'],'Data_res');
            Data_res = permute(Data_res,[2,3,1]);
            Data_res = Data_res(:);
            Data_res_active_mean = mean(Data_res(active_index));
            Data_res(inactive_index) = Data_res_active_mean;
            data_true{1} = reshape(Data_res,options.dim);
            data{1} = data_true{1};
        end
        
        for I = 1:length(data_in_use)
            FD_noisy{I} = cat(4,FD_noisy{data_in_use(I)},data{data_in_use(I)}); % FD = 4D
            FD_true{I} = cat(4,FD_true{data_in_use(I)},data_true{data_in_use(I)});
        end
        
    end
    
    % process FD data (difference and / or average)
    if ~isnan(getOption(seismicOptions,'use_diff_data'))
        FD_noisy_tmp = FD_noisy;
        FD_true_tmp = FD_true;
        clear FD_noisy FD_true;
        for I = 1:length(data_in_use)
            for J = 1:size(FD_noisy_tmp{I},4)-1
                FD_noisy{I}(:,:,:,J) = FD_noisy_tmp{I}(:,:,:,J+1) - FD_noisy_tmp{I}(:,:,:,1);
                FD_true{I}(:,:,:,J) = FD_true_tmp{I}(:,:,:,J+1) - FD_true_tmp{I}(:,:,:,1);
            end
        end
    end
    if ~isnan(getOption(seismicOptions,'avg_data'))
        FD_noisy_tmp = FD_noisy;
        FD_true_tmp = FD_true;
        clear FD_noisy FD_true;
        avgMask = seismicOptions.avg_data;
        for I = 1:length(data_in_use)
            for K = 1:size(avgMask,2)                
                FD_noisy{I}(:,:,K,:) = mean(FD_noisy_tmp{I}(:,:,avgMask(1,K):avgMask(2,K),:),3,'omitnan');
                FD_true{I}(:,:,K,:) = mean(FD_true_tmp{I}(:,:,avgMask(1,K):avgMask(2,K),:),3,'omitnan');
            end
        end
    end
    
    % replace NaN with zeros
    if ~isnan(getOption(seismicOptions,'repzero'))
        for I = 1:length(data_in_use)
            FD_noisy{I}(isnan(FD_noisy{I})) = 0;
            FD_true{I}(isnan(FD_true{I})) = 0;
        end
    else % replace NaN with mean
        for I = 1:length(data_in_use)
            FD_noisy{I}(isnan(FD_noisy{I})) = mean(FD_noisy{I}(~isnan(FD_noisy{I})));
            FD_true{I}(isnan(FD_true{I})) = mean(FD_true{I}(~isnan(FD_true{I})));
        end
    end
    
    % use cumulative timeshift
    if getOption(seismicOptions,'cumsum','false')
        for I = 1:length(data_in_use)
            FD_noisy{I} = cumsum(FD_noisy{I},3);
            FD_true{I} = cumsum(FD_true{I},3);
        end
    end
    
end % if exist('preProcessedSeismicData.mat','file')

% wavelet transform on data
for seismic_time_index = 1 : size(FD_true{1},4)
    
    tmp_measurement = [];
    if sim_index < 1
        tmp_truemeasurement = [];
        tmp_W = [];
    end
    if isfield(seismicOptions,'waveletTransform') && (seismicOptions.waveletTransform ==1) % if do wavelets transform
              
        cD_end_index = 0;
        cA_end_index = 0;
        if sim_index < 1
            tmp_cD_leading_index = [];
            tmp_cD_leading_number = [];
            tmp_cA_leading_index = [];
            tmp_cA_leading_number = [];
            tmp_WDEC_rec = [];
        end
        
        for data_type_index = 1 : length(data_in_use) 
            
            if sim_index < 1 % create leading coefficients
                
                % noisy trace
                current_data = FD_noisy{data_in_use(data_type_index)}(:,:,:,seismic_time_index);
                output = wavelet3D_sparse_representation(current_data,threshold);
                
                tmp_WDEC_rec = [tmp_WDEC_rec; output.WDEC_rec]; 
                
                tmp_cD_leading_index = [tmp_cD_leading_index; output.leading_cD_coeff_index]; %#ok<*AGROW>
                tmp_cD_leading_number = [tmp_cD_leading_number; length(output.leading_cD_coeff_index)];
                current_cD_coeff_index = output.leading_cD_coeff_index; 
                
                tmp_cA_leading_index = [tmp_cA_leading_index; output.leading_cA_coeff_index];
                tmp_cA_leading_number = [tmp_cA_leading_number; length(output.leading_cA_coeff_index)];                
                current_cA_coeff_index = output.leading_cA_coeff_index; 
                
                tmp_measurement = [tmp_measurement; output.leading_coeff(:)]; %#ok<*FNDSB>
                tmp_W = [tmp_W; output.est_noise_level']; % .* ones(length(output.leading_coeff),1)];
                
                % reconstructed data
                FD_rec{data_type_index} = cat(4,FD_rec{data_type_index},waverec3(output.WDEC_rec));
                
            else
                
                current_cD_leading_index = cD_leading_index{seismic_time_index}; %#ok<*USENS>
                current_cD_leading_number = cD_leading_number{seismic_time_index};
                current_cA_leading_index = cA_leading_index{seismic_time_index};
                current_cA_leading_number = cA_leading_number{seismic_time_index};
                cD_start_index = cD_end_index + 1;
                cD_end_index = cD_end_index + current_cD_leading_number(data_type_index);
                current_cD_coeff_index = current_cD_leading_index(cD_start_index:cD_end_index);
                cA_start_index = cA_end_index + 1;
                cA_end_index = cA_end_index + current_cA_leading_number(data_type_index);
                current_cA_coeff_index = current_cA_leading_index(cA_start_index:cA_end_index); 
                
            end
            
            % simulated seismic measurements
            current_data = FD_true{data_in_use(data_type_index)}(:,:,:,seismic_time_index);
            output = wavelet3D_sparse_representation(current_data,0);
            cD_leading_coeff = output.cD_in_vec(current_cD_coeff_index);
            cA_leading_coeff = output.cA_in_vec(current_cA_coeff_index);
            if sim_index < 1  % reference data
                tmp_truemeasurement = [tmp_truemeasurement;  [cA_leading_coeff;cD_leading_coeff]];
            else
                tmp_measurement = [tmp_measurement;  [cA_leading_coeff;cD_leading_coeff]];
            end
            
        end

    else
        
        R0 = sum(R0.^2,3);
        G = sum(G.^2,3);
        R0 = real(R0);
        G = real(G);
        if sim_index < 1
            tmp_truemeasurement = [tmp_truemeasurement; R0(:);G(:)];
            R0 = R0 + randn(size(R0));
            G = G + randn(size(G));
            tmp_measurement = [tmp_measurement;R0(:);G(:)];
            tmp_W = [tmp_W; ones(numel(R0) + numel(G),1)];
        else
            tmp_measurement = [tmp_measurement;R0(:);G(:)];
        end
        
    end

    if sim_index < 1
        % tmp_W is std, while W should always be saved as variance
        W = [W {tmp_W.^2}];
        truemeasurement = [truemeasurement {tmp_truemeasurement}];

        cD_leading_index = [cD_leading_index {tmp_cD_leading_index}];
        cD_leading_number = [cD_leading_number {tmp_cD_leading_number}];

        cA_leading_index = [cA_leading_index {tmp_cA_leading_index}];
        cA_leading_number = [cA_leading_number {tmp_cA_leading_number}];

        FD_WDEC_rec = [FD_WDEC_rec {tmp_WDEC_rec}];
        measurement = [measurement {tmp_measurement}];
    else
        simSeismicData = [simSeismicData; tmp_measurement];
    end
    
end

% plot seismic data
if showFigure
    if strcmp(dataType,'AVO')
        if isfield(seismicOptions,'waveletTransform') && ...
                (seismicOptions.waveletTransform ==1) && sim_index < 1
            plotSeismic(FD_noisy,FD_true,FD_rec,Time,Ricker,Tricker);
        else
            plotSeismic(FD_noisy,FD_true,[],Time,Ricker,Tricker);
        end
    else
        if isfield(seismicOptions,'waveletTransform') && ...
                (seismicOptions.waveletTransform ==1) && sim_index < 1
            plotSeismic(FD_noisy,FD_true,FD_rec);
        else
            plotSeismic(FD_noisy,FD_true);
        end
    end
end

if sim_index < 1
    save trueSeismicData.mat truemeasurement measurement W ...
        cD_leading_index cD_leading_number cA_leading_index cA_leading_number ...
        FD_noisy FD_true FD_rec FD_WDEC_rec;
    if strcmp(dataType,'AVO')
        save('trueSeismicData.mat','Time','Ricker','Tricker','-append');
    end
else
    save simSeismicData.mat simSeismicData;
    save_data = seismicOptions.save_data;
    if save_data
        save('simSeismicData.mat', 'FD_true','FD_noisy','-append');
    end
end

close all

