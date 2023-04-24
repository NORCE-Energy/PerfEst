function H = plotSpatialConVideo(tstep,layer)

if nargin < 2
    layer = 1;
end

% load data
load trueSolutionSmoother.mat;
load inputData.mat;
actnumMask = reshape(options.actnum,options.dim(1),options.dim(2),options.dim(3));
actnumMask(actnumMask==0) = NaN;
na = sum(options.actnum);
time = kalmanOptions.reportTime;
if getOption(kalmanOptions,'thinobs',0) > 0
    time = time(kalmanOptions.measInd);
end
nt = length(time);
time = [0;time];
nf = options.fieldSize;
nd = nf;
if ~strcmp(kalmanOptions.ES_script,'localMDA')
    nd = na;
end

% Contrast - spatial
np = length(tstep);
col = 3;
row = np;
m = reshape(measurement,nd,nt);
m = [zeros(size(m,1),1),m]; % add initial data
A=dir('simulatedDataIter*.mat');
H = {};
cVal = 0.5;
if ~isempty(A) % plot prior, posterior and measurement
    
    load(A(1).name,'simData');
    d_ini_mean = reshape(mean(simData,2),nd,nt); %#ok<*NODEF>
    d_ini_mean = [zeros(size(d_ini_mean,1),1),d_ini_mean]; % add initial data
    load(A(end).name,'simData'); %#ok<*COLND>
    d_end_mean = reshape(mean(simData,2),nd,nt);
    d_end_mean = [zeros(size(d_end_mean,1),1),d_end_mean]; % add initial data
    
    for I = 1:row
        
        H{I} = figure('Name','Tissue spatial contrast','NumberTitle','off',...
             'Units','pixels','Position',[10, 10, 606, 206]); %#ok<*AGROW>
        K=1;
        h = [];
        
        h(1) = subplot(1,col,K);
        data = actnumMask;
        data(data==1) = m(:,tstep(I));
        data = data(:,:,layer);
        %cVal = mean(data(:),'omitnan');
        b = imagesc(data);
        %colormap('bone');
        set(b,'AlphaData',~isnan(data))
        title(['Data @ T\approx ',num2str(round(time(tstep(I))))],'fontsize',10,'fontweight','normal');
        set(gca,'xtick',[])
        set(gca,'ytick',[])
        K = K + 1;
        
        h(2) = subplot(1,col,K);
        data = actnumMask;
        data(data==1) = squeeze(d_ini_mean(:,tstep(I)));
        data = data(:,:,layer);
        b = imagesc(data);
        %colormap('bone');
        set(b,'AlphaData',~isnan(data))
        title('Prior','fontsize',10,'fontweight','normal');
        set(gca,'xtick',[])
        set(gca,'ytick',[])
        K = K + 1;
        
        h(3) = subplot(1,col,K);
        data = actnumMask;
        data(data==1) = squeeze(d_end_mean(:,tstep(I)));
        data = data(:,:,layer);
        b = imagesc(data);
        %colormap('bone');
        set(b,'AlphaData',~isnan(data))
        title('Posterior','fontsize',10,'fontweight','normal');
        set(gca,'xtick',[])
        set(gca,'ytick',[])
        
        if ~isempty(cVal)
            for J = 1:col
                set(h(J),'clim',[0,cVal]);
            end
        end
        
    end
    
end
