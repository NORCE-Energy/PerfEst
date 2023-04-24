function h = plotField(field,varargin)


cl = setProperty(varargin,'cl',[]);
load trueSolutionSmoother.mat;
load inputData.mat;
if exist('initialState.mat','file')
    load('initialState.mat','options','initialState','ref','prm');
    initialOptions = options;
end
if exist('finalState.mat','file')
    load finalState.mat;
end
visc = options.visc0;

actnumMask = reshape(options.actnum,options.dim(1),options.dim(2),options.dim(3));
actnumMask(actnumMask==0) = NaN;

% Estimated fields
switch field
    case 'PERMXART'
        scale = 1e-6;
        T = prm.comp.arterial.gm.perm(1);
        T = scale*T;
        V{1} = initialOptions.permX(:,:,1,1);
        V{2} = options.permX(:,:,1,1);
   case 'PERMXVEN'
        scale = 1e-6;
        T = prm.comp.venous.gm.perm(1);
        T = scale*T;
        V{1} = initialOptions.permX(:,:,1,2);
        V{2} = options.permX(:,:,1,2);
    case 'PERMYART'
        scale = 1e-6;
        T = prm.comp.arterial.gm.perm(2);
        T = scale*T;
        V{1} = initialOptions.permY(:,:,1,1);
        V{2} = options.permY(:,:,1,1);
   case 'PERMYVEN'
        scale = 1e-6;
        T = prm.comp.venous.gm.perm(1);
        T = scale*T;
        V{1} = initialOptions.permY(:,:,1,2);
        V{2} = options.permY(:,:,1,2);        
    case 'PERMQ'
        scale = 1e-6;
        T = 1.5e-9;
        T = scale*T;
        V{1} = initialOptions.permQ;
        V{2} = options.permQ;
    case 'TRANXART'
        scale = 1e-6;
        T = prm.comp.arterial.gm.perm(1);
        T = scale*T/visc; 
        V{1} = initialOptions.permX(2:end,:,1,1);
        V{2} = options.permX(2:end,:,1,1);
    case 'TRANXVEN'
        scale = 1e-6;
        T = prm.comp.venous.gm.perm(1);
        T = scale*T/visc; 
        V{1} = initialOptions.permX(2:end,:,1,2);
        V{2} = options.permX(2:end,:,1,2);
    case 'TRANYART'
        scale = 1e-6;
        T = prm.comp.arterial.gm.perm(2);
        T = scale*T/visc; 
        V{1} = initialOptions.permY(:,1:end-1,1,1);
        V{2} = options.permY(:,1:end-1,1,1);
    case 'TRANYVEN'
        scale = 1e-6;
        T = prm.comp.venous.gm.perm(2);
        T = scale*T/visc; 
        V{1} = initialOptions.permY(:,1:end-1,1,2);
        V{2} = options.permY(:,1:end-1,1,2);
    case 'TRANQ'
        scale = 1e-6;
        T = prm.f*visc/1000;
        T = scale*T/visc; 
        V{1} = initialOptions.permQ;
        V{2} = options.permQ;
    case 'POROART'
        T = 0.05;
        V{1} = initialOptions.porosityArt(:,:,1);
        V{2} = options.porosityArt(:,:,1);
    case 'POROVEN'
        T = 0.1;
        V{1} = initialOptions.porosityVen(:,:,1);
        V{2} = options.porosityVen(:,:,1);
    case 'POROQ'
        T = 0;
        V{1} = initialOptions.porosityQ;
        V{2} = options.porosityQ;
end
   
   
h = figure('Name',['Frog ',field],'NumberTitle','off',...
       'position',[1, 776, 1165, 335]);
hs(1) = subplot(1,3,1);
v =  T.*actnumMask;
b = imagesc(v);
set(b,'AlphaData',~isnan(v))
title([field,' True'])
colorbar; 
set(gca,'fontsize',14,'fontweight','bold')
hs(2) = subplot(1,3,2);
v = V{1}.*actnumMask;
b = imagesc(v);
set(b,'AlphaData',~isnan(v))
title([field,' Initial'])
colorbar; 
set(gca,'fontsize',14,'fontweight','bold')
hs(3) = subplot(1,3,3);
v = V{2}.*actnumMask;
minmax = [min(v(:)),max(v(:))];
b = imagesc(v);
set(b,'AlphaData',~isnan(v))
title([field,' Final'])
colorbar;
set(gca,'fontsize',14,'fontweight','bold')
if isempty(cl), cl = [0,3]*T;
elseif strcmp(cl,'minmax'), cl = minmax; end
for I = 1:3
    try
        set(hs(I),'clim',cl);
    catch
    end
end