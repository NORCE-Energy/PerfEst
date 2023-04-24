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
P_svd = zeros(nxComp,nyComp);
for I = 1:nxComp
    for J = 1:nyComp
               
        c = squeeze(Ccrs(I,J,:));
        R = V*(L\U')*c;
        P_svd(I,J) = max(R) / dt;
        
    end
end
