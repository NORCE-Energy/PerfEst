function [tapering_matrix,output,corr_mtx,corr_mtx_shuffled]= autoAdaLoc(ensemble,simData,options,kalmanOptions,varargin)
%%
% [tapering_matrix,output]= autoAdaLoc(ensemble,simData,varargin)
% Compute a tapering matrix based on sample correlation fields
% 
% -- INPUTS --
% ensemble:             ensemble of model variables
% simData:              simulated data corresponding to ensemble
% options:              options from setupCase.m (see Matlab toolbox manual)
% varargin:             reserved for future development
%
% -- OUTPUTS --
% tapering_matrix:      constructed tapering matrix
% output:               a structure variable used to record other information
%
% Copyright (c) 2010-2018 IRIS, All Rights Reserved.
% $Id:  $
% $DateTime:  $

rp_index = generate_randperm_index('repeat',1,'ensize',size(ensemble,2));
%shuffled_ensemble = ensemble(:,rp_index);
corr_mtx = getCorrMtx(ensemble,simData);
%corr_mtx_shuffled = getCorrMtx(shuffled_ensemble,simData);
corr_mtx_shuffled = getCorrMtx(randn(size(ensemble)),simData);

%
tapering_matrix = zeros(size(corr_mtx));
output = [];

%
num_active = numel(find(options.actnum));
ensize = kalmanOptions.ensembleSize;
nstd = kalmanOptions.nstd;
% if exist free param
starting_index = 0;
if isfield(options,'freeparam')
    disp(' Tapering freeparam ');
    prop_index = 1 : sum(options.freeparam);
    
    thresholding_mode = 'hard';
    [current_tapering,current_output] = tapering_function(corr_mtx(prop_index,:),corr_mtx_shuffled(prop_index,:),ensize,nstd,thresholding_mode,'freepara',1);
    %[current_tapering,current_output] = tapering_function_theoretical(corr_mtx(prop_index,:),ensize);
    output = current_output;
    
    tapering_matrix(prop_index,:) = current_tapering;
    starting_index = starting_index + sum(options.freeparam);
end

%
for para_index = 1 : size(options.staticVar,1) % number of parameter fields
    
    prop_index = starting_index + (1:num_active) + (para_index-1)*num_active;
    
    thresholding_mode = kalmanOptions.thresholding_mode{para_index};
    [current_tapering,current_output] = tapering_function(corr_mtx(prop_index,:),corr_mtx_shuffled(prop_index,:),ensize,nstd,thresholding_mode,'freepara',0);
    
    % if 'est_noise_std' exists, so are 'cutoff_point_lower' and 'cutoff_point_upper'
    if isfield(output,'est_noise_std')
        output.est_noise_std = [output.est_noise_std;current_output.est_noise_std];
        if isfield(output,'cutoff_point_lower')
            output.cutoff_point_lower = [output.cutoff_point_lower;current_output.cutoff_point_lower];
            output.cutoff_point_upper = [output.cutoff_point_upper;current_output.cutoff_point_upper];
        end
    else
        output = current_output;
    end
    
    tapering_matrix(prop_index,:) = current_tapering;
end

output.rp_index = rp_index;

%

