function plotResults(layer)


if nargin < 1
    layer = 1;
end

if ~exist('Figures','dir')
    mkdir('Figures');
end

load inputData.mat; %#ok<*LOAD>
trans = false;
prefix = 'perm';
if strcmp(options.staticVar(1,1:4),'TRAN')
    trans = true;
    prefix = 'trans';
end

h = plotTof(layer);
prtFile = 'Figures/tissue_tof';
print(h,'-r0',prtFile,'-depsc2');
print(h,'-r0',prtFile,'-dpng');

h = plotPerf('layer',layer);
prtFile = 'Figures/tissue_perf';
print(h,'-r0',prtFile,'-depsc2');
print(h,'-r0',prtFile,'-dpng');

% Estimated x- permeability or transmissibility fields
if trans
    plotField('TRANXART','layer',layer); 
else
    plotField('PERMXART','layer',layer); 
end % art
prtFile = ['Figures/tissue_',prefix,'xart'];
print('-r0',prtFile,'-depsc2');
print('-r0',prtFile,'-dpng');
if trans
    plotField('TRANXVEN','layer',layer); 
else
    plotField('PERMXVEN','layer',layer); 
end % ven
prtFile = ['Figures/tissue_',prefix,'xven'];
print('-r0',prtFile,'-depsc2');
print('-r0',prtFile,'-dpng');

% Estimated q-permeability fields
if trans, plotField('TRANQ','cl','minmax','layer',layer); else
          plotField('PERMQ','cl','minmax','layer',layer); end% q
prtFile = ['Figures/tissue_',prefix,'q'];
print('-r0',prtFile,'-depsc2');
print('-r0',prtFile,'-dpng');

% Estimated y- permeability or transmissibility fields
if trans 
    plotField('TRANYART','layer',layer); 
else
    plotField('PERMYART','layer',layer); 
end % art
prtFile = ['Figures/tissue_',prefix,'yart'];
print('-r0',prtFile,'-depsc2');
print('-r0',prtFile,'-dpng');
if trans 
    plotField('TRANYVEN','layer',layer); 
else
    plotField('PERMYVEN','layer',layer); 
end % ven
prtFile = ['Figures/tissue_',prefix,'yven'];
print('-r0',prtFile,'-depsc2');
print('-r0',prtFile,'-dpng');

% Estimated porosity fields
plotField('POROART','layer',layer);
prtFile = 'Figures/tissue_poroart';
print('-r0',prtFile,'-depsc2');
print('-r0',prtFile,'-dpng');
plotField('POROVEN','layer',layer);
prtFile = 'Figures/tissue_poroven';
print('-r0',prtFile,'-depsc2');
print('-r0',prtFile,'-dpng');
plotField('POROQ','cl','minmax','layer',layer);
prtFile = 'Figures/tissue_poroq';
print('-r0',prtFile,'-depsc2');
print('-r0',prtFile,'-dpng');

% Contrast - temporal
%coord = [148, 25; 88, 68; 22, 103]; % left middle right
D = options.dim(1:2); 
if prod(D) >= 9
    C = round(D/2);
    H = round(D(2)/2) - round(D(2)/4);
    V = round(D(1)/2) - round(D(1)/4);
    coord = [ C-[V H] ; C-[V 0]; C-[V -H]; ...
        C-[0 H] ; C      ; C+[0 H] ; ...
        C+[V -H]; C+[V 0]; C+[V H] ; ]; % nine spot
else
    K = 1;
    for I = 1:D(1)
        for J = 1:D(2)
            coord(K,1) = I;  %#ok<*AGROW>
            coord(K,2) = J;
            K = K + 1;
        end
    end
end
%coord = [148, 103; 87, 103; 22, 103; ...
%         148,  25; 87,  25; 22, 25]; % six spot
%coord = [153, 52; 140, 57; 133, 67]; % inlet test
plotTemporalCon(coord,layer);
prtFile = 'Figures/tissue_con_temporal';
print('-r0',prtFile,'-depsc2');
print('-r0',prtFile,'-dpng');

% Plot spatial concentration
tstep = [8 10 12 14];
plotSpatialCon(tstep,layer);
prtFile = 'Figures/tissue_con_spatial';
print('-r0',prtFile,'-depsc2');
print('-r0',prtFile,'-dpng');

% Objective function
if exist('objRealIter0.mat','file')
    plotObjective('objRealIter','y_log_scale',1,'plotType','boxplot');
end

% Plot perfusion estimate
postProcess(options.dim,'mae');
