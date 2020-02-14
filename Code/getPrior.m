function getPrior(prm,data,options)


% upscale level
L = getOption(options,'L',prm.dim);

% find poro in arteries
dimFine = prm.dim;
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
D(D<=0) = NaN;
D = mean(D(:),'omitnan');
K = D^2 / 32;
permArt(art == 1) = K;

% find perm in veins (based on Poiseuille’s equation)
permVen = options.kV*ones(dimFine);
permVen = options.fineMask.*permVen;
D = 2*data.venous.tree.radius.im / 1000;
D(D<=0) = NaN;
D = mean(D(:),'omitnan');
K = D^2 / 32;
permVen(ven == 1) = K;

% inter-compartment perm and poro 
permQ = options.kQ*ones(dimFine);
permQ(art == 1 | ven == 1) = 1e-16;
permQ = options.fineMask.*permQ;
poroQ = options.pQ*ones(dimFine);
poroQ(art == 1 | ven == 1) = 1e-16;
poroQ = options.fineMask.*poroQ;
        
% upscale
poroArtCrs = upscale(poroArt,L);
poroVenCrs = upscale(poroVen,L);
permArtCrs = upscale(permArt,L);
permVenCrs =upscale(permVen,L);
permQCrs = upscale(permQ,L);
poroQCrs = upscale(poroQ,L);
poroArtCrs(isnan(poroArtCrs)) = options.pA; %#ok<*NASGU>
poroVenCrs(isnan(poroVenCrs)) = options.pV;
permArtCrs(isnan(permArtCrs)) = options.kA;
permVenCrs(isnan(permVenCrs)) = options.kV;
permQCrs(isnan(permQCrs)) = options.kQ;
poroQCrs(isnan(poroQCrs)) = options.pQ;

% save 
save('priorTissue.mat','poroArtCrs','poroVenCrs','permArtCrs','permVenCrs','permQCrs','poroQCrs');
disp('priorTissue.mat is created.')