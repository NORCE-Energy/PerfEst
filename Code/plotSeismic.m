function plotSeismic(data,data_true,data_rec,Time,Ricker,Tricker)

% function plotSeismic(data,data_rec,data_true,Time,Ricker,Tricker)
%
% Plot seismic information
%
%%%%%%% INPUT Parameters
%
% data     : seismic data including noise
% data_true: seismic data without noise
% data_rec : reconstructed data (if using wavelets)
% Time     : 3D matrix with travel time information in sec
% Ricker   : values of computed Ricker wavelet
% Tricker  : Time information corresponds to Ricker wavelet
%
% Written by T. Bhakta, 2016; Modified by Xiaodong Luo, 2016; Modified by
% Rolf Lorentzen, 2017
%
% Copyright (c) 2010-2016 IRIS, All Rights Reserved.
% $Id: //depot/rfmatlab/main/seismicHM/Rockphysics Model/get3DAmplitudeData.m#3 $
% $DateTime: 2017/02/06 19:52:27 $


% load input
load('inputData.mat');

plot_save_dir = seismicOptions.plot_save_dir;
plot_truth = seismicOptions.plot_truth;
angle = seismicOptions.Theta;
subfix= seismicOptions.subfix;
seismic_time = seismicOptions.seismic_time;
data_in_use = seismicOptions.data_in_use;

% slices to plot
Xline_number = seismicOptions.Xline_number;
Yline_number  = seismicOptions.Yline_number;
Zline_number  = seismicOptions.Zline_number;

% This information is needed to generate Ricker wavelet ,
% fr is the deominanat frequency
fr = seismicOptions.fr;

% post processing
visible = seismicOptions.visible;
clim = seismicOptions.clim;

if ~strcmp(plot_save_dir(end),'/')
    plot_save_dir = [plot_save_dir '/'];
end
if ~exist(plot_save_dir,'dir')
    mkdir(plot_save_dir);
end

if ~plot_truth
    data_true = [];
end

fs = 14; % font size

% Plot the Ricker wavelet
if nargin > 3
    T= squeeze(Time(1,1,:));
    figure('visible',visible);
    h1=plot(Ricker,Tricker,'m');
    set(gca,'fontsize',fs);
    set(gca,'YDir','reverse'); %
    grid on;
    title(['Ricker wavelet with' num2str(fr) ' Hz central freq'])
    prtFile = [plot_save_dir 'Ricker_wavelet' subfix];
    savePlotTight(prtFile,h1);
    vaflag = 1;
    fact  = 1;
    clplvl  = 1;
    flipy  = 1;
    kolor = [ 0 0 1];
end

for i = 1 : 3 % data, data_true, data_rec
    
    str1 = 'data';
    if i == 1 && ~isempty(data)
        current_data = data;
    elseif i == 2 && ~isempty(data_true)
        current_data = data_true;
        str1 = [str1,'_true'];
    elseif i == 3 && ~isempty(data_rec)
        current_data = data_rec;
        str1 = [str1,'_rec'];
    else
        continue;
    end
    
    for k = 1:length(current_data) % angle or impedance type
        
        if nargin > 3
            str2 = [str1,'_',num2str(angle(data_in_use(k)))]; %#ok<*AGROW>
        elseif strcmp(seismicOptions.dataType,'Timeshift')
            str2 = [str1,'_timshift'];
        else
            if data_in_use(k) == 1
                str2 = [str1,'_pImp'];
            else
                str2 = [str1,'_sImp'];
            end
        end
        plot_data = current_data{k};
        
        for m = 1:size(plot_data,4) % time
            
            str3 = [str2,'_T',num2str(seismic_time(m))];
            plot_data_3D = plot_data(:,:,:,m);
            % remove inactive
            
            if nargin > 3 % assume AVO data
                
                % Xline slice (x-direction fixed)
                for j = 1 : length(Xline_number)
                    slice_data = squeeze(plot_data_3D(Xline_number(j),:,:))'; %#ok<*NODEF>
                    x= 1:size(slice_data,2);
                    figure('visible',visible);
                    plotseis(slice_data,T,x,vaflag,fact,clplvl,flipy,kolor); % plot seismogram without noise
                    set(gca, 'fontsize', fs, 'linewidth', 1)
                    xlabel('Trace no.');
                    set(gca,'XAxisLocation','bottom')
                    ylabel('Time(s)');
                    title(['Inline No. ' num2str(Xline_number(j))]);
                    
                    prtFile = [plot_save_dir str3 '_X' num2str(Xline_number(j)) subfix];
                    savePlotTight(prtFile,gcf);
                end
                close all
                
                % Yline slice (y-direction fixed)
                for j = 1 : length(Yline_number)
                    slice_data = squeeze(plot_data_3D(:,Yline_number(j),:))'; %#ok<*NODEF>
                    x= 1:size(slice_data,2);
                    figure('visible',visible);
                    plotseis(slice_data,T,x,vaflag,fact,clplvl,flipy,kolor); % plot seismogram without noise
                    set(gca, 'fontsize', fs, 'linewidth', 1)
                    xlabel('Trace no.');
                    set(gca,'XAxisLocation','bottom')
                    ylabel('Time(s)');
                    title(['Xline No. ' num2str(Yline_number(j))]);
                    prtFile = [plot_save_dir str3 '_Y' num2str(Yline_number(j)) subfix];
                    savePlotTight(prtFile,gcf);
                end
                close all
                
                % Zline slice (time-direction fixed)
                for j = 1 : length(Zline_number)
                    slice_data = squeeze(plot_data_3D(:,:,Zline_number(j)))'; %#ok<*NODEF>
                    x = 1:size(slice_data,2);
                    y = 1:size(slice_data,1);
                    figure('visible',visible);
                    plotseis(slice_data,y,x,vaflag,fact,clplvl,flipy,kolor); % plot seismogram without noise
                    set(gca, 'fontsize', fs, 'linewidth', 1)
                    xlabel('Trace no.');
                    ylabel('Trace no.')
                    set(gca,'XAxisLocation','bottom')
                    title(['Time = ' num2str(T(Zline_number(j)))]);
                    title(['Time No. ' num2str(Zline_number(j))]);
                    prtFile = [plot_save_dir str3 '_Z' num2str(Zline_number(j)) subfix];
                    savePlotTight(prtFile,gcf);
                end
                close all
                
            else % assume spatial data
                
                % first remove inactive cells
                if ~isnan(getOption(seismicOptions,'avg_data'))
                    actnum_tmp = reshape(options.actnum,options.dim);
                    avgMask = seismicOptions.avg_data;
                    for K = 1:size(avgMask,2)
                        actnum(:,:,K) = max(actnum_tmp(:,:,avgMask(1,K):avgMask(2,K)),[],3);
                    end
                    plot_data_3D(actnum(:) == 0) = NaN;
                else
                    plot_data_3D(options.actnum == 0) = NaN;
                end
                
                % Xline slice (x-direction fixed)
                for j = 1 : length(Xline_number)
                    slice_data = squeeze(plot_data_3D(Xline_number(j),:,:))'; %#ok<*NODEF>
                    figure('visible',visible);
                    eplot(slice_data,0); % plot 
                    set(gca, 'fontsize', fs, 'linewidth', 1)
                    colorbar;
                    xlabel('Y');
                    ylabel('Z');
                    title(['Xline No. ' num2str(Xline_number(j))]);
                    if ~isempty(clim), set(gca,'clim',clim); end
                    prtFile = [plot_save_dir str3 '_X' num2str(Xline_number(j)) subfix];
                    savePlotTight(prtFile,gcf);
                end
                
                % Yline slice (y-direction fixed)
                for j = 1 : length(Yline_number)
                    slice_data = squeeze(plot_data_3D(:,Yline_number(j),:))'; %#ok<*NODEF>
                    figure('visible',visible);
                    eplot(slice_data,0); % plot 
                    set(gca, 'fontsize', fs, 'linewidth', 1)
                    colorbar;
                    xlabel('X');
                    ylabel('Z');
                    title(['Yline No. ' num2str(Yline_number(j))]);
                    if ~isempty(clim), set(gca,'clim',clim); end
                    prtFile = [plot_save_dir str3 '_Y' num2str(Yline_number(j)) subfix];
                    savePlotTight(prtFile,gcf);
                end              
                
                % Zline slice (time-direction fixed)
                for j = 1 : length(Zline_number)
                    slice_data = squeeze(plot_data_3D(:,:,Zline_number(j)))'; %#ok<*NODEF>
                    figure('visible',visible);
                    eplot(slice_data,0); % plot 
                    avg = mean(slice_data(:),'omitnan');
                    set(gca, 'fontsize', fs, 'linewidth', 1)
                    colorbar;
                    xlabel('X');
                    ylabel('Y')
                    title(['Zline No. ' num2str(Zline_number(j)),', Avg: ',num2str(avg,1)]);
                    if ~isempty(clim), set(gca,'clim',clim); end
                    prtFile = [plot_save_dir str3 '_Z' num2str(Zline_number(j)) subfix];
                    savePlotTight(prtFile,gcf);
                end
                
            end
            
        end % for m
    end % for k
end % for i

if strcmp(visible,'off')
    close all;
end


