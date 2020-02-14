if ~existfile('setupForEnKF.mat')
    %% results from Erlend
    load simfullbrainindicator-346x448x61 prm results
    %load simfullbraindisc-346x448x61 disc
    %% calculate contrast
    contrast=results.Cmat.arterial.im+results.Cmat.venous.im;
    
    %% upscale contrast
    upscaleContrast
    
    %% measurements and total porosity
    meas=reshape(upscaledContrast,62*66*20,90);
    totPor=sum(meas,2);
    totPor=totPor/max(totPor);
    totPor=reshape(totPor,62,66,20);
    time=prm.reporttimeline;
    
    %% prepare for simulation
    setupOptions
    clear prm results disc
    %% if time to peak is low, we assume it is an artery otherwise it is a vein in layer 1
    % might consider to follow veins/arteries upward
    [a,b]=max(meas,[],2);
    
    mask=activeBlocks;
    arteriesBottom=reshape(b(1:62*66),62,66);
    arteriesBottom(arteriesBottom<11)=1;
    arteriesBottom(arteriesBottom>1)=0;
    arteriesBottom(mask(:,:,1)==0)=0;
    veinsBottom=reshape(b(1:62*66),62,66);
    veinsBottom(veinsBottom<11)=0;
    veinsBottom(veinsBottom>0)=1;
    
    
    options.mask = zeros(options.nx,options.ny,options.nz,2);
    options.maskQ = zeros(options.nx,options.ny,options.nz);
    
    maskArt=mask;
    maskArt(:,:,1)=false;
    maskArt(:,:,1)=arteriesBottom;
    maskVen=mask;
    maskVen(:,:,1)=false;
    maskVen(:,:,1)=veinsBottom;
    options.mask(:,:,:,1) = maskArt;
    options.mask(:,:,:,2) = maskVen;
    options.maskQ = mask;
    options.maskQ(:,:,1)=false;
    
    %% more preparations
    options.porosityArt=totPor/3;
    options.porosityVen=totPor/3;
    options.porosityQ=totPor/3;
    options.porosityArt(:,:,1)=arteriesBottom.*totPor(:,:,1);%totPor(arteries==1);
    options.porosityVen(:,:,1)=veinsBottom.*totPor(:,:,1);
    options.porosityQ(:,:,1)=0;
    
    % options.porosityArt()=0;
    % options.porosityVen(veins==1)=1;%totPor(veins==1);;
    % options.porosityQ(veins==1)=0;
    
    options.permX(:,:,:,1)=1e-8*totPor+eps;
    options.permX(:,:,:,2)=1e-8*totPor+eps;
    options.permY(:,:,:,1)=1e-8*totPor+eps;
    options.permY(:,:,:,2)=1e-8*totPor+eps;
    options.permZ(:,:,:,1)=1e-8*totPor+eps;
    options.permZ(:,:,:,2)=1e-8*totPor+eps;
    options.permQ=1e-8*totPor+eps;%5e-9*ones(62,66,20); % starte with 5e-9
    
    options.porosityArt(options.mask(:,:,:,1)==0)=0;
    options.porosityVen(options.mask(:,:,:,2)==0)=0;
    options.porosityQ(options.maskQ==0)=0;
    
    options.empiricalAIF=[[0; time(1:end-1); 2*time(end)] [0; meas(2537,2:end)'/totPor(2537); 0]]; % adjust with one second
    options.time=0:0.1:90;
    options.variableList=char('permX','permY','permZ','permQ');
    gridSize=62*66*20
    options.variableLength=[gridSize*2;gridSize*2;gridSize*2;gridSize];
    %% make an ensemble
    save setupForEnKF
end
load setupForEnKF
%% NEED to prepare below this
for I=1:100
    v=zeros(7*62*66*20,1);
    for II=1:5
        corrLen=66*rand(1);
        x= fastGaussian3d([62 66 20], 1, corrLen );
        v(1+(II-1)*gridSize:II*gridSize)=x;
    end
    ensemble(:,I)=v(:);
end   

options.time=0:0.025:90;
reportTimes=round(time*10)/10;
reportTimes=reportTimes(reportTimes>0);
options.reportTime=reportTimes;

%% test funWrapper
for I=1:100
    param=exp(ensemble(:,I));
    [Y,volTracerArt,volTracerCap,volTracerVen,velX,velY,velZ,rateQ,pres,modifiedOptions]=funWrapper(options,param);
    output(:,I)=sum(Y)';%Y(:,end);
    figure(1),clf
    %plot(Y(:,end)),hold on,plot(meas(:,1),'r'),title(num2str(I)),drawnow
    plot(reportTimes,output(:,I)),hold on,plot(time,sum(meas),'r'),title(num2str(I)),drawnow
    I
    if max(output(:,I))>0.1
        keyboard
    end
end

%% EnKF
measunc=0.01*max(sum(meas)')
[meanstate,updEnsemble]=EnKF([ensemble;output],(measunc).^2*eye(size(output,1)),[zeros(size(output,1),size(ensemble,1)) eye(size(output,1))],sum(meas)');

%% run on updated ensemble
figure(1),clf
for I=1:100
    param=exp(updEnsemble(1:size(ensemble,1),I));
    [Y,volTracerArt,volTracerCap,volTracerVen,velX,velY,velZ,rateQ,pres]=funWrapper(options,param);
    newoutput(:,I)=sum(Y)';%Y(:,end);
    figure(1)
    plot(sum(Y),'b'),hold on,plot(sum(meas),'r'),title(num2str(I)),drawnow
    I
end

%% run on meanstate
param=exp(meanstate(1:size(ensemble,1)));
[Y,volTracerArt,volTracerCap,volTracerVen,velX,velY,velZ,rateQ,pres]=funWrapper(options,param);
%figure(1),clf
plot(sum(Y),'k','linewidth',2),hold on,plot(sum(meas),'r'),title('With mean estiamte'),drawnow
%% time development, last ensemble member
% close all
% for I=1:120
%     figure(1),imagesc(reshape(Y(:,I),29,29)),title(num2str(options.reportTime(I))),colorbar,drawnow
%     figure(2),imagesc(Cmattot(:,:,1,I)),title(num2str(time(I))),colorbar,drawnow
%     pause(1)
% end
%%    
save fromWorkflow2 -v7.3
