function [tapering_coefficient,output] = tapering_function(corr_field,corr_field_shuffled,ensize,nstd,thresholding_mode,varargin)
% function tapering_coefficient = tapering_function(corr_field,,corr_field_shuffled,varargin)
%
% -- INPUTS --
% corr_field:               correlation fields
% corr_field_shuffled:      shuffled correlation fields for estimation of noise levels 
% ensize:                   ensemble size
% nstd:                     # of STDs
% thresholding_mode:        available options {'hard','soft'}
% varargin:                 reserved for future development
%
% -- OUTPUTS --
% tapering_coefficient:     tapering coefficients
% output:                   other information
%
% Made by XILU, Jan 2018
% Copyright (c) 2010-2018 IRIS, All Rights Reserved.
% $Id:  $
% $DateTime:  $

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% options for varargin
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% "tf" stands for "tapering function", and should point to an already
% constructed function handle (fh) in the form of
% fh = @(corr_field,corr_field_shuffled,kalmanOptioins,freepara) ......
% If tf is empty, then the default tapering function will be called
%tf = setProperty(varargin,'tf',[]);
tf = setProperty(varargin,'tf','rationalFunc_thresholding');

% "freepara" stands for "free parameter". If it's a free parameter, the the 
% noise std is estimated as 1/sqrt(ensize); otherwise a MAD estimator is
% adopted
freepara = setProperty(varargin,'freepara',0);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if isempty(tf)    
    [tapering_coefficient,output] = three_pieces_thresholding(corr_field,corr_field_shuffled,freepara,ensize,nstd,thresholding_mode);
else
    [tapering_coefficient,output] = feval(tf,corr_field,corr_field_shuffled,freepara,ensize,nstd,thresholding_mode);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% embedded function
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    function [tc, output] = three_pieces_thresholding(cf,cf_shuffled,fp,ensize,nstd,thresholding_mode,varargin)
        % cf = corr_field, nl = noise_level, tc = tapering_coefficient
        % fp = freepara;
       
        
        tc = zeros(size(cf));
        
        for i = 1 : size(cf,2) % size(cf,2) corresponds to observation size
            current_cf = cf(:,i);
            
            if fp
                est_noise_std = 1/sqrt(ensize);
            else
                est_noise_std = median( abs(cf_shuffled(:,i) ) ) / 0.6745; % MAD estimator
            end
            
            %
            cutoff_point_upper = sqrt(2*log(numel(current_cf))) * est_noise_std; % universal threshold rule
            
            if strcmpi(thresholding_mode,'hard')
                cutoff_point_lower = cutoff_point_upper; % hard thresholding
            else %
                % soft thresholding. "nstd" stands for number of STDs (e.g., ko.nstd = 1)
                cutoff_point_lower = nstd * est_noise_std; 
                cutoff_point_upper = 1;
                %cutoff_point_upper = max(abs(current_cf(:)));
            end
            
            if (cutoff_point_upper < cutoff_point_lower)
               error('Lower truncation point is larger than upper one'); 
            end
            
            output.est_noise_std(i) = est_noise_std;
            output.cutoff_point_lower(i) = cutoff_point_lower;
            output.cutoff_point_upper(i) = cutoff_point_upper;
            
            %
            current_tc = zeros(numel(current_cf),1); 
            
            % if (cutoff_point_lower == cutoff_point_upper) => hard thresholding
            % else => soft thresholding
            set_upper = find( abs(current_cf(:)) >  cutoff_point_upper);
            current_tc(set_upper) = 1;
            
            if cutoff_point_lower ~= cutoff_point_upper 
                set_lower = find( abs(current_cf(:)) >= cutoff_point_lower);
                
                set_middle = setdiff(set_lower,set_upper);
                
                %current_tc(set_middle) = ( exp( current_cf(set_middle).^2 ./ (2*est_noise_std) ) - exp( cutoff_point_lower^2 / (2*est_noise_std) ) ) / ...
                %    ( exp( cutoff_point_upper^2 / (2*est_noise_std) ) - exp( cutoff_point_lower^2 / (2*est_noise_std) ) );
                
                current_tc(set_middle) = polynomial_tapering(abs(current_cf(set_middle)), cutoff_point_lower,cutoff_point_upper);                
            end
            
            tc(:,i) = current_tc(:);                  
        end
        
    end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% embeded function
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

end