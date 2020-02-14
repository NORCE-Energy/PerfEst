function [tc, output] = rationalFunc_thresholding(cf,cf_shuffled,fp,ensize,nstd,thresholding_mode,varargin)
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
            output.est_noise_std(i) = est_noise_std;
            
            %
            current_tc = zeros(numel(current_cf),1); 
            
            %
            cutoff_point = sqrt(2*log(numel(current_cf))) * est_noise_std; % universal threshold rule
            output.threshold(i) = cutoff_point;
            
            if strcmpi(thresholding_mode,'hard')
                set_upper = find( abs(current_cf(:)) >  cutoff_point);
                current_tc(set_upper) = 1; %#ok<*FNDSB>
            else %
                current_tc = rationalFunction(1-abs(current_cf(:)),abs(1-cutoff_point));
            end
            
            tc(:,i) = current_tc(:);                  
        end
        
    end