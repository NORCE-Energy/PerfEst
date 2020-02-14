function y=rationalFunction(dist,lc,varargin) 
% the 5th-order rational function, cf. Eq (4.10) of DAO office note
% 96-03R1 (by G. Gaspari and E. Cohn)

% a, b \geq 0

        z = abs(dist(:)) ./ lc;
        
        index_1 = find(z <= 1);
        index_2 = find(z <= 2);
        index_12 = setdiff(index_2,index_1); 
        
        y = zeros(length(z),1);
        
        y(index_1) = 1 - ( (z(index_1).^5) ./ 4 )  ...
                       + ( (z(index_1).^4) ./ 2 ) ...
                       + ( 5 .* (z(index_1).^3) ./ 8 ) ... 
                       - ( 5 .* (z(index_1).^2) ./ 3 );
        
        y(index_12) =   ( (z(index_12).^5) ./ 12 ) ...
                      - ( (z(index_12).^4) ./ 2 ) ...
                      + (5 .* (z(index_12).^3) ./ 8)  ...
                      + (5 .* (z(index_12)) .^2 ./ 3) ... 
                      - (5 .* z(index_12) ) ...
                      - (2 ./ (3 .* z(index_12)) ) + 4;
        
     