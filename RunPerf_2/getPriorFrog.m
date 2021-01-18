function getPriorFrog(prm,data,options)


% upscale level
L = getOption(options,'L',[158,128]);

% find poro in arteries
dimFine = prm.dim(1:2);
dimCrs = floor(dimFine/L); %L; 
poroArt = options.pA*ones(dimFine);
art = data.arterial.tree.bw;
poroArt(art == 1) = options.pVeins;
poroArt = options.fineMask.*poroArt;

% find poro in veins
poroVen = options.pV*ones(dimFine);
ven = data.venous.tree.bw;
poroVen(ven == 1) = options.pVeins;
poroVen = options.fineMask.*poroVen;

% find perm in arteries (based on Poiseuille’s equation)
permArt = options.kA*ones(dimFine);
permArt = options.fineMask.*permArt;
D = 2*data.arterial.tree.radius.im / 1000;
scale = 1;
[I,J] = find(D>0);
for k = 1:length(I)
    xRange = find(art(I(k),:) == 1);
    xDiff = diff(xRange) == 1;
    xDiff = [xDiff,false]; %#ok<*AGROW>
    index = find(xRange == J(k));
    upRange = find(xDiff(index:end)==false,1,'first')-1;
    if isempty(upRange), upRange = 0; end
    upRange = J(k)+upRange;
    downRange = find(xDiff(1:index-1)==false,1,'last')+1;
    if isempty(downRange), downRange = 1; end
    downRange = J(k)-(index-downRange);
    range = downRange:upRange;
    permArt(I(k),range) = max(permArt(I(k),range),scale * D(I(k),J(k))^2 / 32);
end

% find perm in veins (based on Poiseuille’s equation)
permVen = options.kV*ones(dimFine);
permVen = options.fineMask.*permVen;
D = 2*data.venous.tree.radius.im / 1000;
[I,J] = find(D>0);
for k = 1:length(I)
    xRange = find(ven(I(k),:) == 1);
    xDiff = diff(xRange) == 1;
    xDiff = [xDiff,false]; %#ok<*AGROW>
    index = find(xRange == J(k));
    upRange = find(xDiff(index:end)==false,1,'first')-1;
    if isempty(upRange), upRange = 0; end
    upRange = J(k)+upRange;
    downRange = find(xDiff(1:index-1)==false,1,'last')+1;
    if isempty(downRange), downRange = 1; end
    downRange = J(k)-(index-downRange);
    range = downRange:upRange;
    permVen(I(k),range) = max(permVen(I(k),range),scale * D(I(k),J(k))^2 / 32);
end

% inter-compartment perm and poro 
permQ = options.kQ*ones(dimFine);
permQ(art == 1 | ven == 1) = 1e-16;
permQ = options.fineMask.*permQ;
poroQ = options.pQ*ones(dimFine);
poroQ(art == 1 | ven == 1) = 1e-16;
poroQ = options.fineMask.*poroQ;

% find boundary cells
boundaryCells = [];
N = floor(dimFine(2) / L(2));
for j=1:N:dimFine(2)-N+1
    if sum(permArt(end,j:j+N-1)>0) > 0 || sum(permVen(end,j:j+N-1)>0) > 0 
            boundaryCells = [ boundaryCells, floor(1+j/N) ];
    end
end
        
% upscale
poroArtCrs = upscale(poroArt,L);
poroVenCrs = upscale(poroVen,L);
permArtCrs = upscale(permArt,L);
permVenCrs =upscale(permVen,L);
permQCrs = upscale(permQ,L);
poroQCrs = upscale(poroQ,L);
poroArtCrs(isnan(poroArtCrs)) = options.pA;
poroVenCrs(isnan(poroVenCrs)) = options.pV;
permArtCrs(isnan(permArtCrs)) = options.kA;
permVenCrs(isnan(permVenCrs)) = options.kV;
permQCrs(isnan(permQCrs)) = options.kQ;
poroQCrs(isnan(poroQCrs)) = options.pQ;

% total porosity and fraction (phi_a / phi_t)
poroTotCrs = poroVenCrs + poroArtCrs;
poroFCrs = poroArtCrs ./ poroTotCrs;  %#ok<*NASGU>
% phi_v = (1-F)*phi_t
% 0 < F < 1

% save 
save('priorFrog.mat','poroArtCrs','poroVenCrs','permArtCrs','permVenCrs','permQCrs','poroQCrs','boundaryCells');
disp('priorFrog.mat is created.')