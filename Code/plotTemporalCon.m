function f = plotTemporalCon(coord,layer)

if nargin < 1
    layer = 1;
end

load trueSolutionSmoother.mat;
load inputData.mat;

% Contrast - temporal
na = sum(options.actnum);
time = kalmanOptions.reportTime;
if getOption(kalmanOptions,'thinobs',0) > 0
    time = time(kalmanOptions.measInd);
end
nt = length(time);
ne = kalmanOptions.ensembleSize;
nf = options.fieldSize;
nd = nf;
dim = options.dim;
np = size(coord,1);
col = floor(sqrt(np));
row = ceil(np/col);
clear c;
for I = 1:np
    c(I) = sub2ind(dim,coord(I,1),coord(I,2),layer); %#ok<*AGROW>
end
if ~strcmp(kalmanOptions.ES_script,'localMDA')
    nd = na;
    for I = 1:np
      c(I) = sum(options.actnum(1:c(I)));
    end
end
m = reshape(measurement,nd,nt);
A=dir('simulatedDataIter*.mat');
f = figure('Name','Tissue temporal contrast','NumberTitle','off',...
       'Units','centimeters','Position',[10, 10, 18.2, 5.5*col]);
if ~isempty(A) % plot prior, posterior and measurement
      
    load(A(1).name,'simData');
    d_ini = reshape(simData(:,1:ne),nd,nt,ne); %#ok<*NODEF>
    load(A(end).name,'simData');
    d_end = reshape(simData(:,1:ne),nd,nt,ne);
    
    maxVal = 0;
    h = [];
    for I = 1:np
        h(I) = subplot(col,row,I);
        data = squeeze(d_ini(c(I),:,:));
        plot(time,data,'b');
        maxVal = max(maxVal,max(data(:)));
        hold on;
        data = squeeze(d_end(c(I),:,:));
        plot(time,data,'g');
        maxVal = max(maxVal,max(data(:)));
        data = m(c(I),:);
        plot(time,data,'r','linewidth',2);
        set(gca,'xlim',[0 max(time)]);
        maxVal = max(maxVal,max(data(:)));
        title(['   (',num2str(coord(I,1)),',',num2str(coord(I,2)),')'],'fontSize',10,'fontweight','normal');
    end
    for I = 1:np
        set(h(I),'ylim',[0,maxVal]);
    end
   
else % assume we have a 'ref' solution from the initialization
    
    for I = 1:np
        subplot(col,row,I);
        plot(ref(c(I),:),'g','linewidth',2);
        hold on;
        plot(time,m(c(I),:),'r','linewidth',2);
        set(gca,'xlim',[0 max(time)]);
        title(['(',num2str(coord(I,1)),',',num2str(coord(I,2)),')'],'fontSize',10,'fontweight','normal');
    end
    
end
