function P_svd = perfusion_svd(ca,Ccrs,dt)

% find svd of A
r = zeros(size(ca));
r(1) = ca(1);
A = toeplitz(ca,r);
[U,L,V] = svd(A);
threshold = 0.2*max(diag(L));
keepValues = diag(L) >= threshold;
L = L(keepValues,keepValues);
U = U(:,keepValues);
V = V(:,keepValues);
        
% compute perfusion
nxComp = size(Ccrs,1);
nyComp = size(Ccrs,2);
nzComp = size(Ccrs,3);
P_svd = zeros(nxComp,nyComp,nzComp);
for I = 1:nxComp
    for J = 1:nyComp
        for K = 1:nzComp
    
            c = squeeze(Ccrs(I,J,K,:));
            R = V*(L\U')*c;
            P_svd(I,J,K) = max(R) / dt;
        
        end
    end
end
