function H = plotSpatialCon(tstep,layer)

if nargin < 1
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
nf = options.fieldSize;
nd = nf;
if ~strcmp(kalmanOptions.ES_script,'localMDA')
    nd = na;
end

% Contrast - spatial
np = length(tstep);
col = 4;
row = np;
%cVal = 0.5e-4;
m = reshape(measurement,nd,nt);
A=dir('simulatedDataIter*.mat');
if ~isempty(A) % plot prior, posterior and measurement
      
    H = figure('Name','Tissue spatial contrast','NumberTitle','off',...
       'Units','centimeters','Position',[10, 10, 18.2, 5.5*row]);
    
    load(A(1).name,'simData');
    d_ini_mean = reshape(mean(simData,2),nd,nt); %#ok<*NODEF>
    load(A(end).name,'simData'); %#ok<*COLND>
    d_end_mean = reshape(mean(simData,2),nd,nt);
    d_end_std = reshape(std(simData,0,2),nd,nt);
    
    K=1;
    for I = 1:row
        
        h = [];
        h(1) = subplot(row,col,K);
        data = actnumMask;
        data(data==1) = m(:,tstep(I));
        data = data(:,:,layer);
        %cVal = mean(data(:),'omitnan');
        b = imagesc(data);
        %colormap('bone');
        set(b,'AlphaData',~isnan(data))
        if I == 1, title('Data','fontsize',10,'fontweight','normal');  end
        ylabel(['T\approx',num2str(round(time(tstep(I))))],'fontsize',10,'fontweight','normal')
        set(gca,'xtick',[])
        set(gca,'ytick',[])
        K = K + 1;
        
        h(2) = subplot(row,col,K);
        data = actnumMask;
        data(data==1) = squeeze(d_ini_mean(:,tstep(I)));
        data = data(:,:,layer);
        cVal = mean(data(:),'omitnan');
        b = imagesc(data);
        %colormap('bone');
        set(b,'AlphaData',~isnan(data))
        if I == 1, title('Prior','fontsize',10,'fontweight','normal');  end
        set(gca,'xtick',[])
        set(gca,'ytick',[])
        K = K + 1;
        
        h(3) = subplot(row,col,K);
        data = actnumMask;
        data(data==1) = squeeze(d_end_mean(:,tstep(I)));
        data = data(:,:,layer);
        b = imagesc(data);
        %colormap('bone');
        set(b,'AlphaData',~isnan(data))
        if I == 1, title('Posterior','fontsize',10,'fontweight','normal');  end   
        set(gca,'xtick',[])
        set(gca,'ytick',[])
        K = K + 1;
        
        h(4) = subplot(row,col,K);
        data = actnumMask;
        data(data==1) = squeeze(d_end_std(:,tstep(I)));
        data = data(:,:,layer);
        b = imagesc(data);
        %colormap('bone');
        set(b,'AlphaData',~isnan(data))
        if I == 1, title('Std','fontsize',10,'fontweight','normal');  end   
        set(gca,'xtick',[])
        set(gca,'ytick',[])
        K = K + 1;
   
        if length(cVal) > 1
            for J = 1:col
                set(h(J),'clim',[0,cVal]);
            end
        end
        
    end
end
