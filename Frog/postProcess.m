function [P,N] = postProcess(comp,display)


% number of compartments (x,y) to consider
if nargin < 1
    nxComp = 1;
    nyComp = 1;
else
    nxComp = comp(1);
    nyComp = comp(2);
end
if nargin < 2
    display = '';
end

% unit conversion (use mL/min/100mL).
uc = 6000;

% load "true" solution
load('initialState.mat','prm','Ccrs','nt','Vcrs','P_true');
Ccrs(repmat(Vcrs,1,1,nt) >= 0.5) = NaN;

% arterial inlet function (aif)
ca = prm.reportaifval;
T = prm.reporttimeline;

% compute the perfusion using Maximum Slope Model
P_ms = zeros(nxComp,nyComp);
for I = 1:nxComp
    for J = 1:nyComp
   
        dC = squeeze(diff(Ccrs(I,J,:))) ./ diff(T);
        P_ms(I,J) = max(dC) / max(ca);
    
    end
end

% compute the perfusion using SVD
dt = mean(diff(T));
P_svd = perfusion_svd(ca,Ccrs,dt);

% compute initial perfusion from model
load('initialState.mat','initialState','options');
Pinitial = initialState.rateQ;
Pinitial(options.actnum==0) = NaN;
Pinitial(Vcrs >= 0.5) = NaN;
P_model = upscale(Pinitial,[nxComp,nyComp]);

% compute final perfusion from model
load('finalState.mat','state');
Pfinal = state.rateQ;
Pfinal(options.actnum==0) = NaN;
Pfinal(Vcrs >= 0.5) = NaN;
P_da = upscale(Pfinal,[nxComp,nyComp]);

% convert everything
P_true = P_true * uc;
P_ms = P_ms * uc;
P_svd = P_svd * uc;
P_model = P_model * uc;
P_da = P_da * uc;

% compute and return statistics 
N = sum(~isnan(P_da(:)));
P(1) = sum(abs(P_true(:) - P_ms(:)),'omitnan') / N;
P(2) = sum(abs(P_true(:) - P_svd(:)),'omitnan') / N;
P(3) = sum(abs(P_true(:) - P_model(:)),'omitnan') / N;
P(4) = sum(abs(P_true(:) - P_da(:)),'omitnan') / N;

% display results
if strcmp(display,'regions')
    
    Y = []; ytl = {};
    fprintf('I,J       True      MS        DC        PR        PO        \n')
    fprintf('------------------------------------------------------------\n')
    formatSpec = '%1.1i,%1.1i       %3.1e   %3.1e   %3.1e   %3.1e   %3.1e\n';
    for I = 1:nxComp
        for J = 1:nyComp
            fprintf(formatSpec,I,J,P_true(I,J),P_ms(I,J),P_svd(I,J),P_model(I,J),P_da(I,J))
            Y = [Y;P_true(I,J),P_ms(I,J),P_svd(I,J),P_model(I,J),P_da(I,J)]; %#ok<*AGROW>
            ytl = [ytl,[num2str(I),',',num2str(J)]];
        end
    end
    
    % plot bar
    %figure;
    if isrow(Y)
        Y = vertcat(Y,nan(size(Y)));
        bar(Y);
        xlim([0.5 1.5])
        set(gca,'xTickLabel','');
    else
        bar(Y);
        set(gca,'xlim',[0 size(Y,1)+1])
        set(gca,'xTickLabel',ytl);
        xlabel('Region','fontSize',10);
    end
    legend({'True','MS','DC','PR','PO'},'location','northWest','fontSize',9);
    title('Perfusion estimate','fontsize',10,'fontweight','normal');
    ylabel('{\it \fontname{Latin Modern Math} P}  [mL/min/100mL]','fontsize',10);
    
    % save plot
    prtFile = 'Figures/frog_perf_est';
    print('-r0',prtFile,'-depsc2');
    print('-r0',prtFile,'-dpng');
    
    % save plot
    prtFile = 'Figures/frog_perf_rms';
    print('-r0',prtFile,'-depsc2');
    print('-r0',prtFile,'-dpng');
    
elseif strcmp(display,'mae')
    
    %figure;
    P_ = vertcat(P,nan(size(P)));
    bar(P_);
    xlim([0.5 1.5])
    ylim([0,max(P)+ max(P)*0.1])
    title(['Partition: [',num2str(nxComp),',',num2str(nyComp),']'],'fontsize',10,'fontweight','normal');
    legend({'MS','DC','PR','PO'},'location','northeast','fontSize',10);
    set(gca,'xTickLabel','');
    ylabel('Mean absolute error','fontSize',10);
    
    % save plot
    prtFile = 'Figures/frog_perf_mae';
    print('-r0',prtFile,'-depsc2');
    print('-r0',prtFile,'-dpng');
    
end % display
