function plotResults


if ~exist('Figures','dir')
    mkdir('Figures');
end

load inputData.mat;
trans = false;
prefix = 'perm';
if strcmp(options.staticVar(1,1:4),'TRAN')
    trans = true;
    prefix = 'trans';
end

h = plotTof;
prtFile = 'Figures/frog_tof';
print(h,'-r0',prtFile,'-depsc2');
print(h,'-r0',prtFile,'-dpng');

h = plotPerf;
prtFile = 'Figures/frog_perf';
print(h,'-r0',prtFile,'-depsc2');
print(h,'-r0',prtFile,'-dpng');

% Estimated x- permeability or transmissibility fields
if trans, plotField('TRANXART'); else, plotField('PERMXART'); end % art
prtFile = ['Figures/frog_',prefix,'xart'];
print('-r0',prtFile,'-depsc2');
print('-r0',prtFile,'-dpng');
if trans, plotField('TRANXVEN'); else, plotField('PERMXVEN'); end % ven
prtFile = ['Figures/frog_',prefix,'xven'];
print('-r0',prtFile,'-depsc2');
print('-r0',prtFile,'-dpng');

% Estimated q-permeability fields
if trans, plotField('TRANQ','cl','minmax'); else
          plotField('PERMQ','cl','minmax'); end% q
prtFile = ['Figures/frog_',prefix,'q'];
print('-r0',prtFile,'-depsc2');
print('-r0',prtFile,'-dpng');

% Estimated y- permeability or transmissibility fields
if trans, plotField('TRANYART'); else, plotField('PERMYART'); end % art
prtFile = ['Figures/frog_',prefix,'yart'];
print('-r0',prtFile,'-depsc2');
print('-r0',prtFile,'-dpng');
if trans, plotField('TRANYVEN'); else, plotField('PERMYVEN'); end % ven
prtFile = ['Figures/frog_',prefix,'yven'];
print('-r0',prtFile,'-depsc2');
print('-r0',prtFile,'-dpng');

% Estimated porosity fields
plotField('POROART');
prtFile = 'Figures/frog_poroart';
print('-r0',prtFile,'-depsc2');
print('-r0',prtFile,'-dpng');
plotField('POROVEN');
prtFile = 'Figures/frog_poroven';
print('-r0',prtFile,'-depsc2');
print('-r0',prtFile,'-dpng');
plotField('POROQ','cl','minmax');
prtFile = 'Figures/frog_poroq';
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
plotTemporalCon(coord);
prtFile = 'Figures/frog_con_temporal';
print('-r0',prtFile,'-depsc2');
print('-r0',prtFile,'-dpng');

% Plot spatial concentration
tstep = [8 10 12 14];
plotSpatialCon(tstep);
prtFile = 'Figures/frog_con_spatial';
print('-r0',prtFile,'-depsc2');
print('-r0',prtFile,'-dpng');

% Objective function
if exist('objRealIter0.mat','file')
    plotObjective('objRealIter','y_log_scale',1,'plotType','boxplot');
end

% Plot perfusion estimate
postProcess(options.dim(1:2),'mae');
