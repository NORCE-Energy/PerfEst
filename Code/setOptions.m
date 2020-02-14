function options = setOptions(options,state)


% set the options based on values in a state vector
nx = options.nx;
ny = options.ny;
nz = options.nz;
actnum = options.actnum;
na = sum(options.actnum);
if strcmp(options.staticVar(1,1:4),'TRAN')
    actnum = reshape(options.actnum,nx,ny,nz);
    actnumX = cat(1,zeros(1,ny,nz),actnum);
    actnumX = actnumX(:);
    actnumY = cat(2,actnum,zeros(nx,1,nz));
    actnumY = actnumY(:);
    actnumZ = cat(3,actnum,zeros(nx,ny,1));
    actnumZ = actnumZ(:);
    naX = sum(actnumX);
    naY = sum(actnumY);
    naZ = sum(actnumZ);
    nX = nx+1;
    nY = ny+1;
    nZ = nz+1;
else
    actnumX = actnum;
    actnumY = actnum;
    actnumZ = actnum;
    naX = na;
    naY = na;
    naZ = na;
    nX = nx;
    nY = ny;
    nZ = nz;
end

stop = 0;
for I = 1:size(options.staticVar,1)
    
    field = deblank(options.staticVar(I,:));
    switch field
        
        case 'PERMXART'
            V = actnumX;
            start = stop + 1;
            stop = start + naX - 1;
            V(V==1) = state(start:stop);
            options.permX(:,:,:,1)=reshape(exp(V),nX,ny,nz);
        case 'PERMXVEN'
            V = actnumX;
            start = stop + 1;
            stop = start + naX - 1;
            V(V==1) = state(start:stop);
            options.permX(:,:,:,2)=reshape(exp(V),nX,ny,nz);
        case 'PERMYART'
            V = actnumY;
            start = stop + 1;
            stop = start + naY - 1;
            V(V==1) = state(start:stop);
            options.permY(:,:,:,1)=reshape(exp(V),nx,nY,nz);
        case 'PERMYVEN'
            V = actnumY;
            start = stop + 1;
            stop = start + naY - 1;
            V(V==1) = state(start:stop);
            options.permY(:,:,:,2)=reshape(exp(V),nx,nY,nz);
        case 'PERMZART'
            V = actnumZ;
            start = stop + 1;
            stop = start + naZ - 1;
            V(V==1) = state(start:stop);
            options.permZ(:,:,:,1)=reshape(exp(V),nx,ny,nZ);
        case 'PERMZVEN'
            V = actnumZ;
            start = stop + 1;
            stop = start + naZ - 1;
            V(V==1) = state(start:stop);
            options.permZ(:,:,:,2)=reshape(exp(V),nx,ny,nZ);
        case 'PERMQ'
            V = actnum;
            start = stop + 1;
            stop = start + na - 1;
            V(V==1) = state(start:stop);
            options.permQ=reshape(exp(V),nx,ny,nz);
        case 'TRANXART'
            V = actnumX;
            start = stop + 1;
            stop = start + naX - 1;
            V(V==1) = state(start:stop);
            options.permX(:,:,:,1)=reshape(exp(V),nX,ny,nz);
        case 'TRANXVEN'
            V = actnumX;
            start = stop + 1;
            stop = start + naX - 1;
            V(V==1) = state(start:stop);
            options.permX(:,:,:,2)=reshape(exp(V),nX,ny,nz);
        case 'TRANYART'
            V = actnumY;
            start = stop + 1;
            stop = start + naY - 1;
            V(V==1) = state(start:stop);
            options.permY(:,:,:,1)=reshape(exp(V),nx,nY,nz);
        case 'TRANYVEN'
            V = actnumY;
            start = stop + 1;
            stop = start + naY - 1;
            V(V==1) = state(start:stop);
            options.permY(:,:,:,2)=reshape(exp(V),nx,nY,nz);
        case 'TRANZART'
            V = actnumZ;
            start = stop + 1;
            stop = start + naZ - 1;
            V(V==1) = state(start:stop);
            options.permZ(:,:,:,1)=reshape(exp(V),nx,ny,nZ);
        case 'TRANZVEN'
            V = actnumZ;
            start = stop + 1;
            stop = start + naZ - 1;
            V(V==1) = state(start:stop);
            options.permZ(:,:,:,2)=reshape(exp(V),nx,ny,nZ);
        case 'TRANQ'
            V = actnum;
            start = stop + 1;
            stop = start + na - 1;
            V(V==1) = state(start:stop);
            options.permQ=reshape(exp(V),nx,ny,nz);
        case 'POROART'
            V = actnum;
            start = stop + 1;
            stop = start + na - 1;
            V(V==1) = state(start:stop);
            options.porosityArt=reshape(V,nx,ny,nz);
        case 'POROVEN'
            V = actnum;
            start = stop + 1;
            stop = start + na - 1;
            V(V==1) = state(start:stop);
            options.porosityVen=reshape(V,nx,ny,nz);
        case 'POROQ'
            V = actnum;
            start = stop + 1;
            stop = start + na - 1;
            V(V==1) = state(start:stop);
            options.porosityQ=reshape(V,nx,ny,nz);
    end
    
end
