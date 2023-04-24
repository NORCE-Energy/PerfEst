function options = setOptions(options,state)
    

% set the options based on values in a state vector
nx = options.nx;
ny = options.ny;
nz = options.nz;
actnum = options.actnum;
actnum = actnum(:);
na = sum(options.actnum);
if strcmp(options.staticVar(1,1:4),'TRAN')
    actnumX = [zeros(1,ny);reshape(options.actnum,nx,ny)];
    actnumX = actnumX(:);
    actnumY = [reshape(options.actnum,nx,ny),zeros(nx,1)];
    actnumY = actnumY(:);
    naX = sum(actnumX);
    naY = sum(actnumY);
    nX = nx+1;
    nY = ny+1;
    options.permX = zeros(nx+1,ny,nz,2);
    options.permY = zeros(nx,ny+1,nz,2);
    options.permZ = 2e-10*ones(nx,ny,nz+1,2)/options.visc0; % not used
else
    actnumX = actnum;
    actnumY = actnum;
    naX = na;
    naY = na;
    nX = nx;
    nY = ny;
    options.permX = zeros(nx,ny,nz,2);
    options.permY = zeros(nx,ny,nz,2);
    options.permZ = 2e-10*ones(nx,ny,nz,2); % not used
end

V = actnumX;
V(V==1) = state(1:naX);
options.permX(:,:,:,1)=reshape(exp(V),nX,ny,nz);
V = actnumX;
V(V==1) = state(naX+1:2*naX);
options.permX(:,:,:,2)=reshape(exp(V),nX,ny,nz);
V = actnumY;
V(V==1) = state(2*naX+1:2*naX+naY);
options.permY(:,:,:,1)=reshape(exp(V),nx,nY,nz);
V = actnumY;
V(V==1) = state(2*naX+naY+1:2*naX+2*naY);
options.permY(:,:,:,2)=reshape(exp(V),nx,nY,nz);
V = options.actnum;
V(V==1) = state(2*naX+2*naY+1:2*naX+2*naY+na);
options.permQ=reshape(exp(V),nx,ny,nz);
V = options.actnum;
V(V==1) = state(2*naX+2*naY+na+1:2*naX+2*naY+2*na);
options.porosityArt = reshape(V,nx,ny,nz);
V = options.actnum;
V(V==1) = state(2*naX+2*naY+2*na+1:2*naX+2*naY+3*na);
options.porosityVen = reshape(V,nx,ny,nz);
V = options.actnum;
V(V==1) = state(2*naX+2*naY+3*na+1:2*naX+2*naY+4*na);
options.porosityQ = reshape(V,nx,ny,nz);
