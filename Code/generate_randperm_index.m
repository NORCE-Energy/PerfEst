function rp_index = generate_randperm_index(varargin)
%function rp_index = generate_randperm_index(varargin)

% number of repeat experiments (default 20)
num_repeat = setProperty(varargin,'repeat',20);
ensize = setProperty(varargin,'ensize',100);
rp_index = [];

count = 0;

while count < num_repeat
    
    current_perm_index = randperm(ensize);
    
    unpermuted_index = find(current_perm_index == (1:ensize));
    
    if length(unpermuted_index) ==1 % only one index not changed
        while true
            rand_exchange = randi([1, ensize]);
            if rand_exchange ~= unpermuted_index
                current_perm_index(unpermuted_index) = current_perm_index(rand_exchange);
                current_perm_index(rand_exchange) = unpermuted_index;
                break
            end
        end
    elseif length(unpermuted_index) > 1 % more than one index not changed
        shift_index = circshift(unpermuted_index(:),1);
        current_perm_index(unpermuted_index) = current_perm_index(shift_index);
    else % none is not changed
        disp('randperm ok')
    end
    
    rp_index = [rp_index,current_perm_index(:)]; %#ok<*AGROW>
    count = count + 1;
end