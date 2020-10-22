function ensemble = getFrogEnsemble(kalmanOptions,options)


% generate initial ensemble for poro and perm
ensemble = generateInitialEnsemble(kalmanOptions,options);

% generate transmissibilities if given in options
if strcmp(options.staticVar(1,1:4),'TRAN')
    
    ne = size(ensemble,2);
    na = sum(options.actnum);
    nf = prod(options.dim);
    nx = options.nx;
    ny = options.ny;
    threshold = getOption(kalmanOptions,'threshold',1e-10);
    
    % compute transmissibilities in X direction
    permXArt = zeros(nf,ne);
    permXArt(options.actnum==1,:) = exp(ensemble(1:na,:));
    permXArt = reshape(permXArt,nx,ny,ne);
    permXArtMean = mean(permXArt,3);
    permXVen = zeros(nf,ne);
    permXVen(options.actnum==1,:) = exp(ensemble(na+1:2*na,:));
    permXVen = reshape(permXVen,nx,ny,ne);
    permXVenMean = mean(permXVen,3);
    transXArt = zeros(nx+1,ny,ne);
    transXVen = zeros(nx+1,ny,ne);
    transXArt(1,:,:) = permXArt(1,:,:);
    transXArt(nx+1,:,:) = permXArt(nx,:,:); % inlet boundary 
    transXVen(1,:,:) = permXVen(1,:,:);
    transXVen(nx+1,:,:) = permXVen(nx,:,:);
    for I = 2:nx
        for J = 1:ny
            if threshold > 0 && ( (permXArtMean(I,J) > threshold && permXArtMean(I-1,J) < 1e-11) ||...
                    (permXArtMean(I-1,J) > threshold && permXArtMean(I,J) < 1e-11) )
                transXArt(I,J,:) = 1e-16;
            else
                transXArt(I,J,:) = harmmean([permXArt(I-1,J,:),permXArt(I,J,:)],2);
            end
            if threshold > 0 && ( (permXVenMean(I,J) > threshold && permXVenMean(I-1,J) < 1e-11) ||...
                    (permXVenMean(I-1,J) > threshold && permXVenMean(I,J) < 1e-11) )
                transXVen(I,J,:) = 1e-16;
            else
                transXVen(I,J,:) = harmmean([permXVen(I-1,J,:),permXVen(I,J,:)],2);
            end
        end
    end
    transXArt = transXArt / options.visc0;
    transXVen = transXVen / options.visc0;
    index = 1:ny; index(options.boundaryCells) = [];
    transXArt(end,index) = 0; % null out boundary outside arteries
    transXVen(end,index) = 0; % null out boundary outside veins
    actnumX = [zeros(1,ny);reshape(options.actnum,nx,ny)];
    actnumX = actnumX(:);
    transXArt = reshape(transXArt,(nx+1)*ny,ne);
    transXArt = transXArt(actnumX==1,:);
    transXVen = reshape(transXVen,(nx+1)*ny,ne);
    transXVen = transXVen(actnumX==1,:);
    
    % compute transmissibilities in Y direction
    permYArt = zeros(nf,ne);
    permYArt(options.actnum==1,:) = exp(ensemble(2*na+1:3*na,:));
    permYArt = reshape(permYArt,nx,ny,ne);
    permYArtMean = mean(permYArt,3);
    permYVen = zeros(nf,ne);
    permYVen(options.actnum==1,:) = exp(ensemble(3*na+1:4*na,:));
    permYVen = reshape(permYVen,nx,ny,ne);
    permYVenMean = mean(permYVen,3);
    transYArt = zeros(nx,ny+1,ne);
    transYVen = zeros(nx,ny+1,ne);
    transYArt(:,1,:) = permYArt(:,1,:); 
    transYArt(:,ny+1,:) = permYArt(:,ny,:); 
    transYVen(:,1,:) = permYVen(:,1,:);
    transYVen(:,ny+1,:) = permYVen(:,ny,:);
    for I = 1:nx
        for J = 2:ny
            if threshold > 0 && ( (permYArtMean(I,J) > threshold && permYArtMean(I,J-1) < 1e-11) ||...
                    (permYArtMean(I,J-1) > threshold && permYArtMean(I,J) < 1e-11) )
                transYArt(I,J,:) = 1e-16;
            else
                transYArt(I,J,:) = harmmean([permYArt(I,J-1,:),permYArt(I,J,:)],2);
            end
            if threshold > 0 && ( (permYVenMean(I,J) > threshold && permYVenMean(I,J-1) < 1e-11) ||...
                    (permYVenMean(I,J-1) > threshold && permYVenMean(I,J) < 1e-11) )
                transYVen(I,J,:) = 1e-16;
            else
                transYVen(I,J,:) = harmmean([permYVen(I,J-1,:),permYVen(I,J,:)],2);
            end
        end
    end
    transYArt = transYArt / options.visc0;
    transYVen = transYVen / options.visc0;
    actnumY = [reshape(options.actnum,nx,ny),zeros(nx,1)];
    actnumY = actnumY(:);
    transYArt = reshape(transYArt,nx*(ny+1),ne);
    transYArt = transYArt(actnumY==1,:);
    transYVen = reshape(transYVen,nx*(ny+1),ne);
    transYVen = transYVen(actnumY==1,:);
    
    % inter-compartment transmissibilities
    transQ = exp(ensemble(4*na+1:5*na,:));
    transQ = transQ / options.visc0;
    
    % chop out the useful porosity part from ensemble
    poroEns = ensemble(5*na+1:end,:);
    
    % build new ensemble
    ensemble = [log(transXArt);log(transXVen);log(transYArt);log(transYVen);log(transQ);poroEns];
    
end

