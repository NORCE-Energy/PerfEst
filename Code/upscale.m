function [F_coarse, N_el] = upscale(F_fine, D, nanflag)

% Upscale F_fine to a grid with dimension D, where D(1) is in the
% row-direction, D(2) is in the column direction, and so on. Excess cells
% are added to the last upscaled cell. The script omits NaN if nanflag is
% 'omitnan' (default), and includes NaN if nanflag is 'includenan'. Note
% that 'omitnan' divides by the total number of cells in the block,
% including the cells with NaN value. 


D_f = size(F_fine); % the fine grid dimension
if length(D) < length(D_f)
    D = [D,D_f(length(D)+1:end)]; % pad
elseif length(D) > length(D_f)
    D = D(1:length(D_f)); % chop
end

N_c = floor(D_f ./ D); % number of gridcells used for upscaling
if min(N_c) < 1
    error('Upscaled dimension must be smaller than the original dimension!');
end

% initialize
A = F_fine;
N_A = ones(size(A)); 
N_nan = isnan(A); % nan mask
if nargin < 3
    nanflag = 'omitnan';
end

% loop over dimension
for J = 1:length(D)
    
    D_A = size(A);
    B = zeros([D(J),D_A(2:end)]);
    N_B = zeros([D(J),D_A(2:end)]);
    N_C = zeros([D(J),D_A(2:end)]);
    
    % loop over cells
    for I = 1:D(J)
        
        % calculate start and stop index for this dimension
        idx_start = N_c(J)*(I-1)+1;
        idx_stop = N_c(J)*I;
        
        % handle case with excess cells
        if I == D(J) && idx_stop < D_f(J) 
            idx_stop = D_f(J);
        end
        
        % create index cell array
        nd = ndims(A);
        C = repmat({':'},1,nd-1);
        
        % calculate temporary matrix sum       
        V = A(idx_start:idx_stop,C{:});
        W = N_A(idx_start:idx_stop,C{:});
        U = N_nan(idx_start:idx_stop,C{:});
        B(I,C{:}) = sum(V,1, nanflag);
        N_B(I,C{:}) = sum(W,1);
        N_C(I,C{:}) = sum(U,1);
        
    end
    
    % permute and clear for next loop
    A = permute(B,[2:nd,1]);
    N_A = permute(N_B,[2:nd,1]);
    N_nan = permute(N_C,[2:nd,1]);
    clear B N_B N_C;
    
end

% divide by number of elements
N_el = N_A;
N_el(N_el==N_nan) = 0;
F_coarse = A ./ N_el;

% make sure dimensions are correct (preserve unit dimensions)
F_coarse = reshape(F_coarse,D);
N_el = reshape(N_el,D);

