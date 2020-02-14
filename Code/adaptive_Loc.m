function [tapering_matrix,output]= adaptive_Loc(corr_mtx,options,varargin)
%%
% tapering_matrix = adaptive_Loc(corr_mtx,options,varargin)
% Compute a tapering matrix by denoising correlation fields
% 
%%
% Copyright (c) 2010-2016 IRIS, All Rights Reserved.
% $Id: //depot/rfmatlab/main/preProcess/adaptive_Loc.m#6 $
% $DateTime: 2017/12/20 12:37:17 $

% wavelet settings
wavelet_name = setProperty(varargin,'wavelet_name','db3'); % family of wavelet
wavelet_scale = setProperty(varargin,'wavelet_scale',3); % level of decomposition

% num_ones: check dimension by examining how many 1s in options.dim
num_ones = numel(find(options.dim ==1));

if num_ones < 1 % 3D case
    dim = 3;
elseif num_ones < 2 % 2D case
    dim = 2;
else % 1D case
    dim = 1;
end

% post_processing: do post-processing to the tapering coefficients through
% image processing
post_processing = setProperty(varargin,'post_processing',0);

% hybrid thresholding: using threhold values computed both by universal
% rule and correlation std

hybrid_thresholding = setProperty(varargin,'hybrid_thresholding',0);

% tm: threshold multiplier
tm = setProperty(varargin,'tm',1); 
output.tm = tm;

disp(['tm = ' num2str(tm)]);

% dist_constraint: use distance-based constraint to reduce the number of tapering regions
dist_constraint = setProperty(varargin,'dist_constraint',0); 
output.dist_constraint = dist_constraint;
if dist_constraint
    CL = setProperty(varargin,'CL',[1901,1800,13]);% critical length in meters
    output.CL = CL;
end

% target_domain: estimate noise level in wavelet or the original data domain.
% available options: {'wavelet','space'}
target_domain = setProperty(varargin,'target_domain','wavelet');

% TPTR: threshold rule
% available options = {'rigrsure','heursure','sqtwolog','minimaxi','bayesian','identity'}
TPTR = setProperty(varargin,'TPTR','sqtwolog'); % 

%
tapering_matrix = zeros(size(corr_mtx));

active_index = find(options.actnum);
inactive_index = setdiff(1:prod(options.dim),active_index);
num_active = numel(active_index);
curr_corr_field = zeros(options.dim);

load('./inputData.mat','kalmanOptions');

% if exist free param
starting_index = 0;
if isfield(options,'freeparam')
    freeparam_threshold = setProperty(varargin,'freeparam_threshold',kalmanOptions.freeparam_threshold);
    disp(['freeparam threshold = ' num2str(freeparam_threshold)]);
    prop_index = 1 : sum(options.freeparam);
    current_corr_field =  corr_mtx(prop_index,:);
    current_tapering = ones(size(current_corr_field));
    zero_indices = find(abs(current_corr_field(:)) < freeparam_threshold);
    current_tapering(zero_indices) = 0;
    tapering_matrix(prop_index,:) = current_tapering;

    starting_index = starting_index + sum(options.freeparam);
end

%
for obs_index = 1 : size(tapering_matrix,2) % number of observations
    for para_index = 1 : size(options.staticVar,1) % number of parameter fields
        prop_index = starting_index + (1:num_active) + (para_index-1)*num_active;
        
        current_tapering = zeros(prod(options.dim),1);
        
        active_corr_field =  corr_mtx(prop_index,obs_index);
        curr_corr_field(active_index) = active_corr_field; 
        
        %
        prop_mean = mean(curr_corr_field(:));
        curr_corr_field(inactive_index) = prop_mean;
        curr_corr_field = squeeze(curr_corr_field); 
        
        switch lower(target_domain)
            
            case 'space'            
                est_noise_std = median(abs(curr_corr_field(:)))/0.6745; % MAD estimator
                %th = sqrt(2*log(numel(prop_dev))) * est_noise_std;
                
            case 'wavelet'
                % Estimate noise level in the wavelet domain
                switch dim
                    case 1
                        [C,S] = wavedec(curr_corr_field,wavelet_scale,wavelet_name);
                        h_band_index = (numel(C) - S(end-1) + 1):numel(C);
                        h_subband = C(h_band_index);
                    case 2
                        [C,S] = wavedec2(curr_corr_field,wavelet_scale,wavelet_name);
                        hh_band_num = prod(S(end-1,:));
                        h_subband = C((end - hh_band_num + 1):end);
                    case 3
                        WDEC = wavedec3(curr_corr_field,wavelet_scale,wavelet_name);
                        h_subband = WDEC.dec{end}(:);
                    otherwise
                        error('Property must have a dimension!')
                end
                
                %est_noise_std = median(abs(h_subband(:)))/0.6745; % MAD estimator
                est_noise_std = median( abs( h_subband(:) - mean(h_subband(:)) ) )/0.6745; % MAD estimator
            otherwise
                error('Choose either wavelet or space domain')
        end
        
        output.est_noise_std(obs_index,para_index) = est_noise_std;
        
        switch lower(TPTR)
            case 'bayesian'
                std_data = std(curr_corr_field(:));
                threshold = est_noise_std^2 / sqrt(abs(std_data^2 - est_noise_std^2));% threshold                
            case {'rigrsure','heursure','sqtwolog','minimaxi'}
                threshold = thselect(curr_corr_field(:),TPTR) * est_noise_std;
            case 'identity'
                threshold = est_noise_std;
            otherwise
                error('Thresholding rule not supported!')
        end
        
        threshold = threshold*tm; % adjust the threshold value by multiplying it by a factor
        
        if isfield(kalmanOptions,'hybrid_thresholding') && kalmanOptions.hybrid_thresholding
            if isfield(kalmanOptions,'ensembleSize')
                threshold = max(threshold,kalmanOptions.hybrid_tm/sqrt(kalmanOptions.ensembleSize));
            else
                threshold = max(threshold,0.3);
            end
        end

        output.threshold(obs_index,para_index) = threshold;
        
        denoised_corr = curr_corr_field;
        denoised_corr( abs(curr_corr_field(:)) < threshold) = 0;
        current_tapering(find(denoised_corr)) = 1; %#ok<*FNDSB>
        
        current_tapering = reshape(current_tapering, size(curr_corr_field));
        
        % post processing the tapering coeffients; right now is for 2D maps
        if post_processing
            
            % remove objects containing fewer than $num_pixel$ pixels
            num_pixel = setProperty(varargin,'num_pixel',30);
            bw = bwareaopen(current_tapering,num_pixel);
            
            % filling interior gaps
            bw = imfill(bw,'holes');
            
            % smoothing objects (erode twice)
            seD = strel('diamond',1);
            bw = imerode(bw,seD);
            bw = imerode(bw,seD);
            
            % try removing small objects and filling holes again
            current_tapering = bwareaopen(bw,num_pixel);
            current_tapering = imfill(current_tapering,'holes');
        end
        
        % Do distance-based tapering w.r.t the location of the maximum absolute correlation.
        if dist_constraint
            curr_corr_field = reshape(curr_corr_field,options.dim); % revert the effect of "squeeze"
            c = find(abs(curr_corr_field)==max(abs(curr_corr_field(:))));
            [max_i,max_j,max_k] = ind2sub(options.dim,c(1)); % use index "1" in case of multiple max values
            tapering_coeff = getDistLoc_3D('GC',CL,[0 0 0],[max_i,max_j,max_k]); 
            current_tapering = current_tapering .* squeeze(tapering_coeff); 
        end
        
        %
        current_tapering = reshape(current_tapering,numel(current_tapering),1);
        
        tapering_matrix(prop_index,obs_index) = current_tapering(active_index);
    end
end
