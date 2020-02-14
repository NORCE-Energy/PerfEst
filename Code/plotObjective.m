function H = plotObjective(objFileName,varargin)

% function H = plotObjective(objFileName,varargin)
%
% box plots of data mismatch
%
% Input
% -----
% objFileName     : Base name for files with objective measure
% varargin:
% thresholdValue  : Plot a horizonal line, e.g. the reference
% plotTrue        : Add true value to the plots (from trueStaticVar)
% xTick           : xTick values
% y_log_scale     : Use log scale on y-axis if this is 1
% visible         : if set to 'off', the figures are not shown, but saved
% plotType        : Either 'boxplot' or 'errorbar'
%
% Output
% ------
% H               : Vector with figure handles
%
% Copyright (c) 2010-2016 IRIS, All Rights Reserved.
% $Id: //depot/rfmatlab/main/postProcess/plotObjective.m#15 $
% $DateTime: 2019/09/20 13:34:32 $


thresholdValue = setProperty(varargin,'thresholdValue',[]);
plotTrue = setProperty(varargin,'plotTrue',0);
xTick = setProperty(varargin,'xTick', []);
y_log_scale = setProperty(varargin,'y_log_scale',0);
visible = setProperty(varargin,'visible','on');
plotType = setProperty(varargin,'plotType','errorbar');

I = 0;
objData = [];
objSeisData = [];
objLabel = {};
fileName = [objFileName,num2str(I),'.mat'];

while exist(fileName,'file')

    load(fileName);
    
    if exist('individualObjInfo','var')
        objSeisData = [objSeisData,individualObjInfo{2}'];
        objData = [objData,individualObjInfo{1}'];
    else
        objData = [objData,objReal];
    end
    objLabel = [objLabel,num2str(I)];
    
    I = I + 1;
    fileName = [objFileName,num2str(I),'.mat'];
    
end

% remove the scaling
load inputData.mat;
scale = getOption(kalmanOptions,'scaling',1);
objData = objData / scale(1)^2;

%if y_log_scale
%    objData = log10(objData);
%    objSeisData = log10(objSeisData);
%end

% plot production mismatch

H(1) = figure('visible',visible);
if strcmp(plotType,'boxplot')
    boxplot(objData,'Labels',objLabel,'Whisker',10);
elseif strcmp(plotType,'errorbar')
    M = mean(objData);
    S = std(objData);
    X = 0:length(M)-1;
    errorbar(X,M,3*S,'linewidth',2);
    set(gca,'xlim',[-1,length(M)]);
    yl = ylim;
    set(gca,'ylim',[yl(1)-yl(1)*0.01,yl(2)+yl(2)*0.01])
end

if y_log_scale
    yl = get(gca,'ylim');
    set(gca,'ylim',[0,yl(2)]);
    set(gca,'yScale','log');
end

if isempty(thresholdValue) && plotTrue
    load forwSim.mat;
    load trueSolutionSmoother.mat;
    y = y(:); %#ok<*NODEF>
    trueData = y(H(:,2));
    objTrue = getDataMismatch(trueData,sqrtW(W),measurement);
    thresholdValue = objTrue;  
end
if ~isempty(thresholdValue)
    %if y_log_scale
    %    thresholdValue = log10(thresholdValue);
    %end
    xl = get(gca,'xlim');
    yl = get(gca,'ylim');
    hold on;
    plot(xl,[thresholdValue,thresholdValue],'--k');
    if thresholdValue >= yl(2)
        set(gca,'ylim',[yl(1),thresholdValue + 0.1*(diff(yl))]);
    elseif thresholdValue <= yl(1)
        set(gca,'ylim',[thresholdValue - 0.1*(diff(yl)),yl(2)]);
    end
end

xlabel('Iteration number');
title('Data mismatch');
fs = 14;
%if y_log_scale
%    ylabel('log_{10}', 'interpreter', 'tex');
%end
set(gca,'fontsize',fs,'fontweight','bold');

if ~isempty(xTick) 
    xTickLabel = '';
    for i = 1:length(xTick)
        xTickLabel = [xTickLabel {num2str(xTick(i))}];
    end
    set(gca,'xTick',xTick,'xTickLabel',xTickLabel);
end

if ~exist('Figures','dir')
   mkdir('Figures'); 
end
prtFile = ['Figures/',objFileName];
print('-r0',prtFile,'-depsc2');
print('-r0',prtFile,'-dpng');

% plot seismic mismatch

if ~isempty(objSeisData)
    
    H(2) = figure('visible',visible);
    if strcmp(plotType,'boxplot')
        boxplot(objSeisData,'Labels',objLabel,'Whisker',15);
    elseif strcmp(plotType,'errorbar')
        M = mean(objSeisData);
        S = std(objSeisData);
        X = 0:length(M)-1;
        errorbar(X,M,3*S,'linewidth',2);
        set(gca,'xlim',[-1,length(M)]);
        yl = ylim;
        set(gca,'ylim',[yl(1)-yl(1)*0.01,yl(2)+yl(2)*0.01])
    end
    
    if y_log_scale
        yl = get(gca,'ylim');
        set(gca,'ylim',[0,yl(2)]);
        set(gca,'yScale','log');
    end
    
    xlabel('Iteration number');
    title('Data mismatch for seismic data');
    fs = 14;
    %if y_log_scale
    %    ylabel('log_{10}', 'interpreter', 'tex');
    %end
    set(gca,'fontsize',fs,'fontweight','bold');
    
    if ~isempty(xTick)
        xTickLabel = '';
        for i = 1:length(xTick)
            xTickLabel = [xTickLabel {num2str(xTick(i))}];
        end
        set(gca,'xTick',xTick,'xTickLabel',xTickLabel);
    end
    
    if ~exist('Figures','dir')
        mkdir('Figures');
    end
    prtFile = ['Figures/',objFileName,'_seis'];
    print('-r0',prtFile,'-depsc2');
    print('-r0',prtFile,'-dpng');
    
end
