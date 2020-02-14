function kalmanOptions = Localc_2D()

% Localc_2D()
%
% Output: 
% -------
% kalmanOptions:        updated kalmanOptions
%
% --
% locfunc:              string specifying localization type. Available options {'FD','GD'}. FD:Furrer-Bengtsson, GD:Gaspari-Cohn
% cl:                   a vector specifying critical length in x and y direction
% theta:                specifying rotation angle in degree based on anisotropy ellipse
% obsType:              observation type, e.g., WOPR/WWPR/WBHP etc.
% eclipse_filename:     elcipse file name to extract reservoir description information
%
% Copyright (c) 2010-2014 IRIS, All Rights Reserved.
% $Id: //depot/rfmatlab/main/preProcess/Localc_2D.m#6 $
% $DateTime: 2017/12/14 17:08:12 $

visible = 'off';

load inputData

locfunc=kalmanOptions.locfunc;
cl=kalmanOptions.locrange;
Theta=kalmanOptions.locangle;
obsType=kalmanOptions.loctype;

index=variableIndex(options);
kalmanOptions.variable = index.name;
kalmanOptions.index = [index.first; index.last]';

eclipse_filename = kalmanOptions.trueSimName;

load(kalmanOptions.initialEnsemble); %#ok<*NODEF>
Ne=size(ensemble,2); 

xdim=options.dim(1);
ydim=options.dim(2);

% get last .X file
R = dir([eclipse_filename,'.X0*']);
Res=readXfile(R(end).name);
Init=readXfile([eclipse_filename '.INIT']);

%pre-allocation
% INTEHEAD(17) =>  NWELL = number of wells
% INTEHEAD(18) =>  NCWMA = maximum number of completions per well
% INTEHEAD(28) =>  NZWEL = no of 8-character words per well in ZWEL array (= 3)
% INTEHEAD(33) => NICON = no of data elements per completion in ICON array (default 19)

%rotation angle defined based on anisotropy ellipse
Theta=(pi/180)*Theta*ones(Res.INTEHEAD(17),1);

location=[];
type=[];
F=[];
Rho_F=zeros(xdim,ydim,Res.INTEHEAD(17)); 
Rho_G=zeros(xdim,ydim,Res.INTEHEAD(17)); 
Dist=zeros(xdim,ydim,Res.INTEHEAD(17)); 
Distx=zeros(xdim,ydim,Res.INTEHEAD(17));
Disty=zeros(xdim,ydim,Res.INTEHEAD(17));

wn=1;
for w=1:double(Res.INTEHEAD(17))
    actnum=options.actnum;
    wname=Res.ZWEL(wn,:);
    %wname=Res.ZWEL(1+3*(w-1),:);
    
    %--
    % in ICON, 
    % Item 2 - I-coordinate 
    % Item 3 - J-coordinate 
    % Item 4 - K-coordinate 
    %--
    WI=double(Res.ICON((w-1)*Res.INTEHEAD(18)*Res.INTEHEAD(33)+2)); 
    WJ=double(Res.ICON((w-1)*Res.INTEHEAD(18)*Res.INTEHEAD(33)+3));
    
    % not sure how to handle inactive layers in a general way
    % maybe we could use varargin to pass a (case-dependent) function 
    % handle to deal with this situation
    % so far I just assume that all layers are active (by xilu).
    WK=double(Res.ICON((w-1)*Res.INTEHEAD(18)*Res.INTEHEAD(33)+4)); 
    %calculate average gridblock size in x and y direction
    DX=double(mean(Init.DX));
    DY=double(mean(Init.DY));
    
    %build a transformation matrix based on user knowledge about the
    %filed, wells, drainage area and prior information
    Rotate(1,1)=cos(Theta(w));
    Rotate(1,2)=sin(Theta(w));
    Rotate(2,1)=-1.*sin(Theta(w));
    Rotate(2,2)=cos(Theta(w));
    
    
    for i=1:xdim
        for j=1:ydim
            
            l(1)=(WI-i)*DX;
            l(2)=(WJ-j)*DY;
            lt=Rotate*l';
            d=sqrt(sum(lt.^2));
            
            %Furrer-Bengtsson
            h1=sqrt((lt(1)/cl(1))^2+((lt(2)/cl(2))^2));
            if h1>1
                rtemp=0;
            else
                rtemp=1-(1.5*h1)+(0.5*h1^3);
            end
            
            rho_F=1/(1+1/Ne*(1+1/(rtemp^2)));
            
            %Gaspari-Chon
            ratio=sqrt((lt(1)/cl(1))^2+((lt(2)/cl(2))^2));
            h1=sqrt((lt(1)/cl(1))^2+((lt(2)/cl(2))^2));
            h2=sqrt((lt(1)/(2*cl(1)))^2+((lt(2)/(2*cl(2)))^2));
            
            if h1<=1 && h2<=1
                rho_G=(-1/4)*ratio^5+(1/2)*ratio^4+(5/8)*ratio^3-(5/3)*ratio^2+1;
            elseif h1>1 && h2 <=1
                rho_G=(1/12)*ratio^5-(1/2)*ratio^4+(5/8)*ratio^3+(5/3)*ratio^2-5*ratio+4-(2/3)*ratio^(-1);
            elseif h1>1 && h2 >1
                rho_G=0;
            end
            
            Distx(i,j,w)=l(1);
            Disty(i,j,w)=l(2);
            Dist(i,j,w)=d;
            Rho_F(i,j,w)=rho_F;
            Rho_G(i,j,w)=rho_G;
        end
    end
    
    %plot localization functions
    if strmatch(locfunc,'FD')
        figure('visible', visible);
        set(gcf,'Name',wname);
        act=reshape(actnum((WK-1)*xdim*ydim+1:WK*xdim*ydim,1),xdim,ydim);
        %act=reshape(actnum((WK-1)*xdim*ydim+1:WK*xdim*ydim),xdim,ydim);
        Rho_F_temp=Rho_F(:,:,w);
        Rho_F_temp(act==0)=NaN;
        Rho_F_temp(WI,WJ)=NaN;
        eplot(Rho_F_temp(:,:));
        axis equal
        axis tight
        colorbar
        saveas(gca,[deblank(wname),'_F.png'])
        close
    else
        figure('visible', visible);
        set(gcf,'Name',wname);
        act=reshape(actnum((WK-1)*xdim*ydim+1:WK*xdim*ydim,1),xdim,ydim);
        %act=reshape(actnum((WK-1)*xdim*ydim+1:WK*xdim*ydim),xdim,ydim);
        Rho_G_temp=Rho_G(:,:,w);
        Rho_G_temp(act==0)=NaN;
        Rho_G_temp(WI,WJ)=NaN;
        eplot(Rho_G_temp(:,:));
        axis equal
        axis tight
        colorbar
        saveas(gca,[deblank(wname),'_G.png'])
        close 
    end
    
    wn=wn+Res.INTEHEAD(28); % 
    
    %alphaDB is used in EnKF.m to perform a distance-based localization on Kalman gain.
    location=strvcat(location,wname);
    %alphaDB.type=strvcat(alphaDB.type,'WGOR WWCT');
    type=strvcat(type,obsType);
    if strmatch(locfunc,'FD')
        % Furrer-bengtsson
        func=reshape(Rho_F(:,:,w),xdim*ydim,1);
    else
        %Gaspari-chon
        func=reshape(Rho_G(:,:,w),xdim*ydim,1);
    end
    
    % vertical localization
    if isfield(kalmanOptions,'wells') && isfield(kalmanOptions,'comp') && ...
            isfield(kalmanOptions,'layers')
        
        wells = kalmanOptions.wells;
        comp = kalmanOptions.comp;
        layers = kalmanOptions.layers;
        
        idx = strfind(deblank(wells),deblank(wname));
        idx = find(~cellfun(@isempty,idx));
        
        vec = [];
        for I = 1:length(layers)
            vec = [vec;comp(idx,I)*ones(layers(I)*xdim*ydim,1)];
        end
        func=repmat(func,options.dim(3),1);
        func = vec.*func;
    end
    
    F=[F,func(actnum==1)];
    
end

%In addtion to alphaDB following information should be saved into
%Kalmanoptions in the inputData (due to current implementation of
%localization scheme(could be improved))
%kalmanOptions.distanceLoc = 'alphaDB';
%updated_kalmanOptionis = kalmanOptions;

%save inputData kalmanOptions -append

kalmanOptions.alphaDB.location = location;
kalmanOptions.alphaDB.type = type;
kalmanOptions.alphaDB.func = F;
save('inputData.mat', 'kalmanOptions','-append');
%save alphaDB alphaDB
save localization Distx Disty Dist Rho_F Rho_G


