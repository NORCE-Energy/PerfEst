function [corr_mtx,model_zero_spread_index,data_zero_spread_index] = getCorrMtx(ensemble,simData,varargin)
% [corr_mtx,model_zero_spread_index,data_zero_spread_index] = getCorrMtx(ensemble,simData,varargin)
% Get correlation matrix between model variables and simulated data
%
% -- INPUTS --
% ensemble: ensemble of reservoir models
% simData: corresponding simulated observations
%
% -- OUTPUTS --
% corr_mtx: correlation matrix between ensemble and simData
% % Copyright (c) 2010-2014 IRIS, All Rights Reserved.
% $Id: //depot/rfmatlab/main/preProcess/getCorrMtx.m#3 $
% $DateTime: 2018/04/25 15:43:43 $

ne = size(simData,2);
if size(ensemble,2) ~= ne
   error('Sizes of ensemble and simData are not equal!')
end

% corr_mtx = zeros(size(ensemble,1),size(simData,1));
% 
% for i = 1: size(ensemble,1)
%     for j = 1: size(simData,1)
%         tmp = corrcoef(ensemble(i,:),simData(j,:)); % tmp is a 2x2 matrix
%         corr_mtx(i,j) = tmp(1,2);
%     end
% end


% normalization_file = setProperty(varargin,'normalization_file','./trueSeismicDataSmoother.mat');
% 
% if ~isempty(normalization_file)
%     load(normalization_file,'W');
%     simData = simData ./ repmat(W,1,ne);
% end

std_model = std(ensemble,[],2);
std_data = std(simData,[],2);

model_zero_spread_index = find(std_model<1e-6);% zero or tiny model ensemble spread
data_zero_spread_index = find(std_data<1e-6);% zero or tiny model ensemble spread

normalized_ensemble = (ensemble - mean(ensemble,2)*ones(1,ne)) ./ (std_model * ones(1,ne));
normalized_simData = (simData - mean(simData,2) * ones(1,ne)) ./ (std_data * ones(1,ne));

corr_mtx = ( normalized_ensemble * normalized_simData' ) ./ ne;


if ~isempty(model_zero_spread_index)
    corr_mtx(model_zero_spread_index,:) = 0;
end

if ~isempty(data_zero_spread_index)
    corr_mtx(:,data_zero_spread_index) = 0;
end
