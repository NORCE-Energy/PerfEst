function I = localization(ensemble,x1,Ud,X,obsTypeEff,obsLocationEff,weights)

% function I = localization(ensemble,x1,Ud,X,obsTypeEff,obsLocationEff,weights)
% 
% Main function for computing localization. Currently distance based and
% correlation based localization is computed. It is possible to apply
% distance based on production data and correlation based on seismic data.
% 
% Input:
% ------
% ensemble:              Ensemble of state vectors
% xl:                    Singular values
% Ud:                    From SVD 
% X:                     Simulated minus real observations
% obsTypeEff:            Production data types 
% obsLocationEff:        Production data location
% weights:               Weights from AGSMDA
%
% Output: 
% -------
% I:                     KalmanGain * X
%
% Copyright (c) 2010-2014 IRIS, All Rights Reserved.
% $Id: //depot/rfmatlab/main/preProcess/localization.m#8 $
% $DateTime: 2019/02/22 14:05:44 $


load trueSolutionSmoother;
load inputData;

if ~isfield(kalmanOptions,'ensembleSize')
    N = size(ensemble,2);
    if isfield(kalmanOptions, 'append_mean') && kalmanOptions.append_mean
        N = N - 1;
    end
else
    N = kalmanOptions.ensembleSize;
end
N_M = size(ensemble,1);
deltaM=zeros(N_M,N);
if nargin > 6
    meanEn = ensemble*weights';
else
    meanEn = mean(ensemble(:,1:N),2);
end
for j=1:N
    deltaM(:,j)=ensemble(:,j)-meanEn;
end
% I = zeros(N_M,1);
% K = deltaM*x1*Ud';

if isfield(kalmanOptions,'distanceLoc') && kalmanOptions.distanceLoc == 1 % distance-based localization
    
    coeff = deltaM*x1*Ud';
    coeff = distance_Loc(coeff,obsTypeEff,obsLocationEff);
    I = coeff*X;
    disp('-----------------------------------')
    disp('Distance-based localization')
    disp('-----------------------------------')
    
elseif isfield(kalmanOptions,'denoisingLoc') && kalmanOptions.denoisingLoc == 1
    
    disp('-----------------------------------')
    disp('Calculating tapering matrix')
    disp('-----------------------------------')
    
    coeff=deltaM*x1;
    proj_data = Ud' * X;
    corr_mtx = getCorrMtx(ensemble(:,1:N),proj_data);
    target_domain = 'wavelet';
    if isfield(kalmanOptions,'target_domain')
        target_domain = kalmanOptions.target_domain;
    end
    tm = 1;
    if isfield(kalmanOptions,'tm')
        tm = kalmanOptions.tm;
    end
    TPTR = 'sqtwolog';
    if isfield(kalmanOptions,'TPTR')
        TPTR = kalmanOptions.TPTR;
    end
    freeparam_threshold = 0.1;
    if isfield(kalmanOptions,'freeparam_threshold')
        freeparam_threshold = kalmanOptions.freeparam_threshold;
    end
    hybrid_thresholding = 0;
    if isfield(kalmanOptions,'hybrid_thresholding')
        hybrid_thresholding = kalmanOptions.hybrid_thresholding;
    end
    [tapering_matrix,~]= adaptive_Loc(corr_mtx,options,'tm',tm,...
        'target_domain',target_domain,'TPTR',TPTR,...
        'freeparam_threshold',freeparam_threshold,...
        'hybrid_thresholding',hybrid_thresholding);
    
    %save(['./correlation_mtx_' num2str(iter) '.mat'],'corr_mtx','tapering_matrix','output');
    
    disp('-----------------------------------')
    disp('Denoising-based localization')
    disp('-----------------------------------')
    
    coeff = tapering_matrix .* coeff;
    I = coeff*(Ud'*X);
    
elseif isfield(kalmanOptions,'autoAdaLoc') && kalmanOptions.autoAdaLoc == 1
    
    disp('-----------------------------------')
    disp('Calculating tapering matrix')
    disp('-----------------------------------')
    
    coeff=deltaM*x1;
    proj_data = Ud' * X;
    
    % we need to always use initial ensemble to avoid member correlations
    %load('ensemble0.mat');
    tapering_matrix = autoAdaLoc(ensemble(:,1:N),proj_data,options,kalmanOptions);
    
    disp('------------------------------------')
    disp('Denoising-based autoAda localization')
    disp('------------------------------------')
    
    coeff = tapering_matrix .* coeff;
    I = coeff*(Ud'*X);
    
else
    
    I = (deltaM*x1) * (Ud'*X); 
    
end


