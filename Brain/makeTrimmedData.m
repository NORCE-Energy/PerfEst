%function initializeBrain
%% set up for contrast data from Erlend's simulation
if ~exist('trimmedData.mat','file')
    disp('start loading data')
    tic
    warning off;
    load('/home/gen/Data/simfullbrainindicator-384x384x256.mat','results','prm');
    load('/home/gen/Data/simfullbrainsolve-384x384x256','disc')
    warning on;
    toc
    disp('data is loaded')
    
    Cart = results.Cmat.arterial.im;
    [nxF,nyF,nzF,nt] = size(Cart); %#ok<*ASGLU>
    Cven = results.Cmat.venous.im;
    % ioVeins = logical((Cart(:,:,7) > 7e-11) + (Cven(:,:,100) > 1e-10));
    % art = data.arterial.tree.bw;
    % art(ioVeins) = 0;
    % art = repmat(art,1,1,nt);
    % ven = data.venous.tree.bw;
    % ven(ioVeins) = 0;
    % ven = repmat(ven,1,1,nt);
    %Cart(art == 1) = Cart(art==1)/0.05;
    %Cven(ven == 1) = Cven(ven == 1)/0.1;
    %Cven(art == 1) = 0;
    Cfine = Cart + Cven;
    perfObs = results.perf;
    [maxCfine,argMaxFine] = max(Cfine,[],4);
    options.fineMask = NaN(nxF,nyF,nzF);
    options.fineMask(argMaxFine>1) = 1;
    for I=1:size(Cfine,4)
        Cfine(:,:,:,I) = options.fineMask.*Cfine(:,:,:,I);
    end
    perfObs = options.fineMask.*perfObs;
%     totConc=sum(abs(Cfine),4);
%     [numIremoved,numIendRemoved,numJremoved,numJendRemoved,numKremoved,numKendRemoved,comprTotConc]=removeUnactiveLayers(totConc);
%     Cfine=Cfine(numIremoved+1:end-numIendRemoved,numJremoved+1:end-numJendRemoved,numKremoved+1:end-numKendRemoved,:);
%     perfObs=perfObs(numIremoved+1:end-numIendRemoved,numJremoved+1:end-numJendRemoved,numKremoved+1:end-numKendRemoved);
    save('trimmedData','-v7.3','Cfine','perfObs','Cart','Cven','prm','disc') 
end
