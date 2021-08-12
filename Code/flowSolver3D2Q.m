%Compute velocities for a 2-comparment model over a 3-dimensional 
% rectangular region given porosites and permeabilities.  
%Boundary conditions are unit pressure at "down" for comparment one, 
% zero pressure at "down" for compartment two, and no-flow for remaining
% boundaries.
%
%%%%Input:
%
% xL: Extension in x-direction [m]
% yL: Extension in y-direction [m]
% zL: Extension in z-direction [m]
% nx: Number of grid cells in x-direction
% ny: Number of grid cslls in y-direction
% nz: Number of grid cslls in z-direction
% visc: matrix of 2*nx*ny*nz viscosities [Pa s]
% permX: matrix of 2*nx*ny*nz x-permeabilities [m^2]
% permY: matrix of 2*nx*ny*nz y-permeabilities [m^2]
% permZ: matrix of 2*nx*ny*nz z-permeabilities [m^2]
% permQ: matrix of nx*ny*nz inter-compartment conductivities
%        (dimensionless, "permeability per interface area")
%
%%%%Optional input:
%
% mask: matrix of 2*nx*ny*nz of active/inactive flags (1/0)
% maskQ: matrix of nx*ny*nz of active/inactive flags (1/0)
%
%           north
%        _________
%       |         |
%  west |  up     | east
%       |    /    |
%       |    down |
%       |_________|
%           south
%
%%%%Output:
%
% pres: matrix of 2*nx*ny*nz cell-center pressures [Pa]
% velX: matrix of 2*(nx+1)*ny*nz cell-edge velocities [m/s]
% velY: matrix of 2*nx*(ny+1)*nz cell-edge velocities [m/s]
% velY: matrix of 2*nx*ny*(nz+1) cell-edge velocities [m/s]
% rateQ: matrix of nx*ny*nx interface rates [1/s]
%
%             velY(i,j+1,k)
%              ___________    velZ(i,j,k+1)
%             |           |
%             |           |
% velX(i,j,k) |  p(i,j,k) | velX(i+1,j,k)
%             |           |
%             |           |
%             |___________|
%  velZ(i,j,k)
%              velY(i,j,k)
%
%
function [pres, velX, velY, velZ, rateQ, termSrcFlux, termSinkFlux] = flowSolver3D2Q(options, extSrcTerm, sampleMask)

  nx = options.nx;
  ny = options.ny;
  nz = options.nz;
  
  dx = options.xL/nx;
  dy = options.yL/ny;
  dz = options.zL/nz;
  
  velX = zeros(nx+1,ny,nz,2);
  velY = zeros(nx,ny+1,nz,2);
  velZ = zeros(nx,ny,nz+1,2);
  
  xCon = zeros(nx+1,ny,nz,2);
  yCon = zeros(nx,ny+1,nz,2);
  zCon = zeros(nx,ny,nz+1,2);
  qCon = zeros(nx,ny,nz);
  
  if (size(options.permX,1) == nx && size(options.permY,2) == ny && size(options.permZ,3) == nz)
    %disp('### flowSolver3D2Q:  permeabilities detected ###');
    isTransNotPerm=0;
    mobX = options.permX./options.visc;
    mobY = options.permY./options.visc;
    mobZ = options.permZ./options.visc;
    mobQ = 2.0*options.permQ./(options.visc(:,:,:,1)+options.visc(:,:,:,2));
    xCon(2:nx,:,:,:) = (dy*dz/dx)*2*mobX(1:nx-1,:,:,:).*mobX(2:nx,:,:,:)./(mobX(1:nx-1,:,:,:)+mobX(2:nx,:,:,:));
    yCon(:,2:ny,:,:) = (dx*dz/dy)*2*mobY(:,1:ny-1,:,:).*mobY(:,2:ny,:,:)./(mobY(:,1:ny-1,:,:)+mobY(:,2:ny,:,:));
    zCon(:,:,2:nz,:) = (dx*dy/dz)*2*mobZ(:,:,1:nz-1,:).*mobZ(:,:,2:nz,:)./(mobZ(:,:,1:nz-1,:)+mobZ(:,:,2:nz,:));
    %qCon(:,:) = 0.5*(dy/dx+dx/dy)*mobQ;
    qCon(:,:,:) = dy*dx*dz*mobQ;
  elseif (size(options.permX,1) == nx+1 && size(options.permY,2) == ny+1 && size(options.permZ,3) == nz+1)
    %disp('### flowSolver3D2Q:  transmissibilities detected ###');
    isTransNotPerm=1;
    xCon(2:nx,:,:,:) = (dy*dz/dx)*options.permX(2:nx,:,:,:);
    yCon(:,2:ny,:,:) = (dx*dz/dy)*options.permY(:,2:ny,:,:);
    zCon(:,:,2:nz,:) = (dx*dy/dz)*options.permZ(:,:,2:nz,:);
    qCon(:,:,:) = dy*dx*dz*options.permQ;
  else
    disp('### flowSolver3D2Q:  inconsistency detected ###');
    assert(1 == 0);
  end
  
  if (isfield(options,'mask'))
    mask = options.mask;
    xCon(2:nx,:,:,:) = xCon(2:nx,:,:,:).*mask(1:nx-1,:,:,:).*mask(2:nx,:,:,:);
    yCon(:,2:ny,:,:) = yCon(:,2:ny,:,:).*mask(:,1:ny-1,:,:).*mask(:,2:ny,:,:);
    zCon(:,:,2:nz,:) = zCon(:,:,2:nz,:).*mask(:,:,1:nz-1,:).*mask(:,:,2:nz,:);
    qCon(:,:,:) = qCon(:,:,:).*mask(:,:,:,1).*mask(:,:,:,2);
  end
  if (isfield(options,'maskQ'))
    qCon(:,:,:) = qCon(:,:,:).*options.maskQ(:,:,:);
  end
  
  if size(sampleMask,1) > 0
%     sampleMask=reshape(sampleMask,nx,ny,nz);
%     xCon(2:nx,:,:,1) = xCon(2:nx,:,:,1).*sampleMask(1:nx-1,:,:).*sampleMask(2:nx,:,:);
%     yCon(:,2:ny,:,1) = yCon(:,2:ny,:,1).*sampleMask(:,1:ny-1,:).*sampleMask(:,2:ny,:);
%     zCon(:,:,2:nz,1) = zCon(:,:,2:nz,1).*sampleMask(:,:,1:nz-1).*sampleMask(:,:,2:nz);
    
%     xCon(2:nx,:,:,2) = 0;
%     yCon(:,2:ny,:,2) = 0;
%     zCon(:,:,2:nz,2) = 0;
    %qCon(:,:,:) = 0;
  end
  
  nn = nx*ny*nz;
  
  B = zeros(2*nn,9);
  d = [-nn -nx*ny -nx -1 0 1 nx nx*ny nn]';
  
  B(1:nn,1) = -reshape(qCon(:,:,:),1,[])';
  B(nn+1:2*nn,9) = -reshape(qCon(:,:,:),1,[])';
  
  B(1:nn,2) = -reshape(zCon(:,:,1:nz,1),1,[])';
  B(nn+1:2*nn,2) = -reshape(zCon(:,:,1:nz,2),1,[])';
  B(1:nn,8) = -reshape(zCon(:,:,2:nz+1,1),1,[])';
  B(nn+1:2*nn,8) = -reshape(zCon(:,:,2:nz+1,2),1,[])';
  
  B(1:nn,3) = -reshape(yCon(:,1:ny,:,1),1,[])';
  B(nn+1:2*nn,3) = -reshape(yCon(:,1:ny,:,2),1,[])';
  B(1:nn,7) = -reshape(yCon(:,2:ny+1,:,1),1,[])';
  B(nn+1:2*nn,7) = -reshape(yCon(:,2:ny+1,:,2),1,[])';
  
  B(1:nn,4) = -reshape(xCon(1:nx,:,:,1),1,[])';
  B(nn+1:2*nn,4) = -reshape(xCon(1:nx,:,:,2),1,[])';
  B(1:nn,6) = -reshape(xCon(2:nx+1,:,:,1),1,[])';
  B(nn+1:2*nn,6) = -reshape(xCon(2:nx+1,:,:,2),1,[])';
  
% %   if size(extSrcTerm,1) > 0
% %     B(:,1) = 0;
% %     B(:,9) = 0;
% %   end
  B(:,5) = -B(:,1)-B(:,2)-B(:,3)-B(:,4)-B(:,6)-B(:,7) -B(:,8)-B(:,9);
  %ddMax=max(B(:,5));
  B(:,5) = B(:,5) + (B(:,5) == 0);
  
  %The "strange" behaviour of spdiags (m=n or m>n versus m<n) ...
  %B(:,1) = circshift(B(:,1),-nn);
  B(:,2) = circshift(B(:,2),-nx*ny);
  B(:,3) = circshift(B(:,3),-nx);
  B(:,4) = circshift(B(:,4),-1);
  B(:,6) = circshift(B(:,6),1);
  B(:,7) = circshift(B(:,7),nx);
  B(:,8) = circshift(B(:,8),nx*ny);
  %B(:,9) = circshift(B(:,7),nn);
 
  %pressure boundaries
  prsBndMtx = zeros(nx,ny,nz,2);
  prsBndRhs = zeros(nx,ny,nz,2);
  

  [isPrsBnd, prsBndVal] = setFaceBoundaries(options);
  % prsBndVal ~ (dim1, dim2, x/y/z, left/right, Q1/Q2)
  % isPrsBnd ~ (x/y/z, left/right, Q1/Q2)
  

  for q=1:2
    if isPrsBnd(1,1,q) %West
        if (isTransNotPerm)
          xCon(1,:,:,q) = (dy*dz/(0.5*dx))*options.permX(1,:,:,q);
        else
          xCon(1,:,:,q) = (dy*dz/(0.5*dx))*mobX(1,:,:,q);
        end
        if (isfield(options,'mask'))
          xCon(1,:,:,q) = xCon(1,:,:,q).*mask(1,:,:,q);
        end
        prsBndMtx(1,:,:,q) = prsBndMtx(1,:,:,q) + xCon(1,:,:,q);
        prsBndRhs(1,:,:,q) = prsBndRhs(1,:,:,q) + permute(prsBndVal(1:ny,1:nz,1,1,q),[3,1,2]).*xCon(1,:,:,q);
    end
    if isPrsBnd(1,2,q) %East
        if (isTransNotPerm)
          xCon(nx+1,:,:,q) = (dy*dz/(0.5*dx))*options.permX(nx+1,:,:,q);
        else
          xCon(nx+1,:,:,q) = (dy*dz/(0.5*dx))*mobX(nx,:,:,q);
        end
        if (isfield(options,'mask'))
          xCon(nx+1,:,:,q) = xCon(nx+1,:,:,q).*mask(nx,:,:,q);
        end
        prsBndMtx(nx,:,:,q) = prsBndMtx(nx,:,:,q) + xCon(nx+1,:,:,q);
        prsBndRhs(nx,:,:,q) = prsBndRhs(nx,:,:,q) + permute(prsBndVal(1:ny,1:nz,1,2,q),[3,1,2]).*xCon(nx+1,:,:,q);
    end
    if isPrsBnd(2,1,q) %South
        if (isTransNotPerm)
          yCon(:,1,:,q) = (dx*dz/(0.5*dy))*options.permY(:,1,:,q);
        else
          yCon(:,1,:,q) = (dx*dz/(0.5*dy))*mobY(:,1,:,q);
        end
        if (isfield(options,'mask'))
          yCon(:,1,:,q) = yCon(:,1,:,q).*mask(:,1,:,q);
        end
        prsBndMtx(:,1,:,q) = prsBndMtx(:,1,:,q) + yCon(:,1,:,q);
        prsBndRhs(:,1,:,q) = prsBndRhs(:,1,:,q) + permute(prsBndVal(1:nx,1:nz,2,1,q),[1,3,2]).*yCon(:,1,:,q);
    end
    if isPrsBnd(2,2,q) %North
        if (isTransNotPerm)
          yCon(:,ny+1,:,q) = (dx*dz/(0.5*dy))*options.permY(:,ny+1,:,q);
        else
          yCon(:,ny+1,:,q) = (dx*dz/(0.5*dy))*mobY(:,ny,:,q);
        end
        if (isfield(options,'mask'))
          yCon(:,ny+1,:,q) = yCon(:,ny+1,:,q).*mask(:,ny,:,q);
        end
        if (isfield(options,'faceMaskNorth'))
          yCon(:,ny+1,:,q) = yCon(:,ny+1,:,q).*reshape(options.faceMaskNorth(:,:,q),nx,1,nz,1);
        end
        prsBndMtx(:,ny,:,q) = prsBndMtx(:,ny,:,q) + yCon(:,ny+1,:,q);
        prsBndRhs(:,ny,:,q) = prsBndRhs(:,ny,:,q) + permute(prsBndVal(1:nx,1:nz,2,2,q),[1,3,2]).*yCon(:,ny+1,:,q);
    end
    if isPrsBnd(3,1,q) %Down
        if (isTransNotPerm)
          zCon(:,:,1,q) = (dx*dy/(0.5*dz))*options.permZ(:,:,1,q);
        else
          zCon(:,:,1,q) = (dx*dy/(0.5*dz))*mobZ(:,:,1,q);
        end
        if (isfield(options,'mask'))
          zCon(:,:,1,q) = zCon(:,:,1,q).*mask(:,:,1,q);
        end
        prsBndMtx(:,:,1,q) = prsBndMtx(:,:,1,q) + zCon(:,:,1,q);
        prsBndRhs(:,:,1,q) = prsBndRhs(:,:,1,q) + prsBndVal(1:nx,1:ny,3,1,q).*zCon(:,:,1,q);
    end
    if isPrsBnd(3,2,q) %Up
        if (isTransNotPerm)
          zCon(:,:,nz+1,q) = (dx*dy/(0.5*dz))*options.permZ(:,:,nz+1,q);
        else
          zCon(:,:,nz+1,q) = (dx*dy/(0.5*dz))*mobZ(:,:,nz,q);
        end
        if (isfield(options,'mask'))
          zCon(:,:,nz+1,q) = zCon(:,:,nz+1,q).*mask(:,:,nz,q);
        end
        prsBndMtx(:,:,nz,q) = prsBndMtx(:,:,nz,q) + zCon(:,:,nz+1,q);
        prsBndRhs(:,:,nz,q) = prsBndRhs(:,:,nz,q) + prsBndVal(1:nx,1:ny,3,2,q).*zCon(:,:,nz+1,q);
    end
  end
  B(:,5) = B(:,5) + reshape(prsBndMtx,1,[])';
  
  % Terminal sources (arterial)
  termCoeff=zeros(nn,1);
  termSrcRhs=zeros(nn,1);
  if (isfield(options,'termSrc'))
    i=options.termSrc.i;
    j=options.termSrc.j;
    k=options.termSrc.k;
    nglob=nx*ny*(k-1)+nx*(j-1)+i;
    for n=1:size(i,2)
      termCoeff(nglob(n))=termCoeff(nglob(n))+options.termSrc.perm(n);
      termSrcRhs(nglob(n))=termSrcRhs(nglob(n))+options.termSrc.perm(n)*options.termSrc.p0(n);
    end
  end 
  B(1:nn,5) = B(1:nn,5) + termCoeff;
  
  % Terminal sinks (venular)
  termCoeff=zeros(nn,1);
  termSinkRhs=zeros(nn,1);
  if (isfield(options,'termSink'))
    i=options.termSink.i;
    j=options.termSink.j;
    k=options.termSink.k;
    nglob=nx*ny*(k-1)+nx*(j-1)+i;
    for n=1:size(i,2)
      termCoeff(nglob(n))=termCoeff(nglob(n))+options.termSink.perm(n);
      termSinkRhs(nglob(n))=termSinkRhs(nglob(n))+options.termSink.perm(n)*options.termSink.p0(n);
    end
  end 
  B(nn+1:2*nn,5) = B(nn+1:2*nn,5) + termCoeff;
  
  %msg=['--- B-mtx-diag:  min=', num2str(min(B(:,5)),'%10.5e'), ' max=', num2str(max(B(:,5)),'%10.5e')];
  %disp(msg);
  
  mtx = spdiags(B,d,2*nn,2*nn);
  
  rhs = zeros(2*nn,1);

  if size(extSrcTerm,1) > 0
    rhs(1:nn) = extSrcTerm(:,1);
    rhs(nn+1:2*nn) = extSrcTerm(:,2);
  else
    rhs(1:nn) = (reshape(prsBndRhs(:,:,:,1),1,[]))'+termSrcRhs;
    rhs(nn+1:2*nn) = (reshape(prsBndRhs(:,:,:,2),1,[]))' + termSinkRhs;
  end
  
  x = mtx\rhs;
  
  pres(:,:,:,1) = reshape(x(1:nn),nx,ny,nz);
  pres(:,:,:,2) = reshape(x(nn+1:2*nn),nx,ny,nz);
  
  if size(extSrcTerm,1) > 0
    rateQ=zeros(0);
    return
  end
  
  % Compute velocities:
  
  velX(1,:,:,:) = 0.0;
  velX(nx+1,:,:,:) = 0.0;
  velX(2:nx,:,:,:) = -xCon(2:nx,:,:,:).*(pres(2:nx,:,:,:)-pres(1:nx-1,:,:,:));
  for q=1:2
    if isPrsBnd(1,1,q) %West
      velX(1,:,:,q) = -xCon(1,:,:,q).*(pres(1,:,:,q) - permute(prsBndVal(1:ny,1:nz,1,1,q),[3,1,2]).*ones(1,ny,nz));
    end
    if isPrsBnd(1,2,q) %East
      velX(nx+1,:,:,q) = -xCon(nx+1,:,:,q).*(permute(prsBndVal(1:ny,1:nz,1,2,q),[3,1,2]).*ones(1,ny,nz) - pres(nx,:,:,q));
    end
  end
  velX = velX/(dy*dz);
  
  velY(:,1,:,:) = 0.0;
  velY(:,ny+1,:,:) = 0.0;
  velY(:,2:ny,:,:) = -yCon(:,2:ny,:,:).*(pres(:,2:ny,:,:)-pres(:,1:ny-1,:,:));
  for q=1:2
    if isPrsBnd(2,1,q) %South
      velY(:,1,:,q) = -yCon(:,1,:,q).*(pres(:,1,:,q) - permute(prsBndVal(1:nx,1:nz,2,1,q),[1,3,2]).*ones(nx,1,nz));
    end
    if isPrsBnd(2,2,q) %North
      velY(:,ny+1,:,q) = -yCon(:,ny+1,:,q).*(permute(prsBndVal(1:nx,1:nz,2,2,q),[1,3,2]).*ones(nx,1,nz) - pres(:,ny,:,q));
    end
  end
  velY = velY/(dx*dz);
  
  velZ(:,:,1,:) = 0.0;
  velZ(:,:,nz+1,:) = 0.0;
  velZ(:,:,2:nz,:) = -zCon(:,:,2:nz,:).*(pres(:,:,2:nz,:)-pres(:,:,1:nz-1,:));
  for q=1:2
    if isPrsBnd(3,1,q) %Down
      velZ(:,:,1,q) = -zCon(:,:,1,q).*(pres(:,:,1,q) - prsBndVal(1:nx,1:ny,3,1,q).*ones(nx,ny));
    end

    if isPrsBnd(3,2,q) %Up
      velZ(:,:,nz+1,q) = -zCon(:,:,nz+1,q).*(prsBndVal(1:nx,1:ny,3,2,q).*ones(nx,ny) - pres(:,:,nz,q));
    end
  end
  velZ = velZ/(dx*dy);
  
  rateQ = -qCon.*(pres(:,:,:,2)-pres(:,:,:,1));
  rateQ = rateQ/(dx*dy*dz);
  
  %Flux for terminal sources
  termSrcFlux = zeros(0);
  if (isfield(options,'termSrc'))
    i=options.termSrc.i;
    j=options.termSrc.j;
    k=options.termSrc.k;
    termSrcFlux=zeros(size(i,2),1);
    for n=1:size(i,2)  
      termSrcFlux(n)= - options.termSrc.perm(n) * (pres(i(n),j(n),k(n),1) - options.termSrc.p0(n));
    end
  end 
  
  %Flux for terminal sinks
  termSinkFlux = zeros(0);
  if (isfield(options,'termSink'))
    i=options.termSink.i;
    j=options.termSink.j;
    k=options.termSink.k;
    termSinkFlux=zeros(size(i,2),1);
    for n=1:size(i,2)  
      termSinkFlux(n)= options.termSink.perm(n) * (options.termSink.p0(n) - pres(i(n),j(n),k(n),2));
    end
  end 
  
  end

  
function [isPrsBnd, prsBndVal] = setFaceBoundaries(options)

  nx = options.nx;
  ny = options.ny;
  nz = options.nz;
  nmax = max([nx,ny,nz]);

  isPrsBnd = zeros(3,2,2); % x/y/z, left/right, Q1/Q2
  prsBndVal = zeros(nmax,nmax,3,2,2);
  
  %Note that no-flow condition overides any prs cnd given ...
  isNoFlow = zeros(3,2,2); % x/y/z, left/right, Q1/Q2
  isNoFlow(1,1,1) = isfield(options,'noFlowWestQ1') && options.noFlowWestQ1 == 1;
  isNoFlow(1,2,1) = isfield(options,'noFlowEastQ1') && options.noFlowEastQ1 == 1;
  isNoFlow(2,1,1) = isfield(options,'noFlowSouthQ1') && options.noFlowSouthQ1 == 1;
  isNoFlow(2,2,1) = isfield(options,'noFlowNorthQ1') && options.noFlowNorthQ1 == 1;
  isNoFlow(3,1,1) = isfield(options,'noFlowDownQ1') && options.noFlowDownQ1 == 1;
  isNoFlow(3,2,1) = isfield(options,'noFlowUpQ1') && options.noFlowUpQ1 == 1;
  isNoFlow(1,1,2) = isfield(options,'noFlowWestQ2') && options.noFlowWestQ2 == 1;
  isNoFlow(1,2,2) = isfield(options,'noFlowEastQ2') && options.noFlowEastQ2 == 1;
  isNoFlow(2,1,2) = isfield(options,'noFlowSouthQ2') && options.noFlowSouthQ2 == 1;
  isNoFlow(2,2,2) = isfield(options,'noFlowNorthQ2') && options.noFlowNorthQ2 == 1;
  isNoFlow(3,1,2) = isfield(options,'noFlowDownQ2') && options.noFlowDownQ2 == 1;
  isNoFlow(3,2,2) = isfield(options,'noFlowUpQ2') && options.noFlowUpQ2 == 1;

  %Face-wise pressure conditions:
  if isfield(options,'prsValueWestQ1') && ~isNoFlow(1,1,1)
    isPrsBnd(1,1,1) = 1;
    prsBndVal(1:ny,1:nz,1,1,1) = options.prsValueWestQ1*ones(ny,nz);
  end
  if isfield(options,'prsValueWestQ2') && ~isNoFlow(1,1,2)
    isPrsBnd(1,1,2) = 1;
    prsBndVal(1:ny,1:nz,1,1,2) = options.prsValueWestQ2*ones(ny,nz);
  end
  if isfield(options,'prsValueEastQ1') && ~isNoFlow(1,2,1)
    isPrsBnd(1,2,1) = 1;
    prsBndVal(1:ny,1:nz,1,2,1) = options.prsValueEastQ1*ones(ny,nz);
  end
  if isfield(options,'prsValueEastQ2') && ~isNoFlow(1,2,2)
    isPrsBnd(1,2,2) = 1;
    prsBndVal(1:ny,1:nz,1,2,2) = options.prsValueEastQ2*ones(ny,nz);
  end
  if isfield(options,'prsValueSouthQ1') && ~isNoFlow(2,1,1)
    isPrsBnd(2,1,1) = 1;
    prsBndVal(1:nx,1:nz,2,1,1) = options.prsValueSouthQ1*ones(nx,nz);
  end
  if isfield(options,'prsValueSouthQ2') && ~isNoFlow(2,1,2)
    isPrsBnd(2,1,2) = 1;
    prsBndVal(1:nx,1:nz,2,1,2) = options.prsValueSouthQ2*ones(nx,nz);
  end
  if isfield(options,'prsValueNorthQ1') && ~isNoFlow(2,2,1)
    isPrsBnd(2,2,1) = 1;
    prsBndVal(1:nx,1:nz,2,2,1) = options.prsValueNorthQ1*ones(nx,nz);
  end
  if isfield(options,'prsValueNorthQ2') && ~isNoFlow(2,2,2)
    isPrsBnd(2,2,2) = 1;
    prsBndVal(1:nx,1:nz,2,2,2) = options.prsValueNorthQ2*ones(nx,nz);
  end
  if isfield(options,'prsValueDownQ1') && ~isNoFlow(3,1,1)
    isPrsBnd(3,1,1) = 1;
    prsBndVal(1:nx,1:ny,3,1,1) = options.prsValueDownQ1*ones(nx,ny);
  end
  if isfield(options,'prsValueDownQ2') && ~isNoFlow(3,1,2)
    isPrsBnd(3,1,2) = 1;
    prsBndVal(1:nx,1:ny,3,1,2) = options.prsValueDownQ2*ones(nx,ny);
  end
  if isfield(options,'prsValueUpQ1') && ~isNoFlow(3,2,1)
    isPrsBnd(3,2,1) = 1;
    prsBndVal(1:nx,1:ny,3,2,1) = options.prsValueUpQ1*ones(nx,ny);
  end
  if isfield(options,'prsValueUpQ2') && ~isNoFlow(3,2,2)
    isPrsBnd(3,2,2) = 1;
    prsBndVal(1:nx,1:ny,3,2,2) = options.prsValueUpQ2*ones(nx,ny);
  end

  %Vertex-wise pressure conditions:
  isPrsVtx = zeros(2,2,2,2); % xL/xR, yL/yR, zL,zR, Q1/Q2
  prsVtx = zeros(2,2,2,2);
  if (isfield(options,'prsValQ1LLL'))
    isPrsVtx(1,1,1,1) = 1;
    prsVtx(1,1,1,1) = options.prsValQ1LLL;
  end
  if (isfield(options,'prsValQ1LLR'))
    isPrsVtx(1,1,2,1) = 1;
    prsVtx(1,1,2,1) = options.prsValQ1LLR;
  end
  if (isfield(options,'prsValQ1LRL'))
    isPrsVtx(1,2,1,1) = 1;
    prsVtx(1,2,1,1) = options.prsValQ1LRL;
  end
  if (isfield(options,'prsValQ1LRR'))
    isPrsVtx(1,2,2,1) = 1;
    prsVtx(1,2,2,1) = options.prsValQ1LRR;
  end
  if (isfield(options,'prsValQ1RLL'))
    isPrsVtx(2,1,1,1) = 1;
    prsVtx(2,1,1,1) = options.prsValQ1RLL;
  end
  if (isfield(options,'prsValQ1RLR'))
    isPrsVtx(2,1,2,1) = 1;
    prsVtx(2,1,2,1) = options.prsValQ1RLR;
  end
  if (isfield(options,'prsValQ1RRL'))
    isPrsVtx(2,2,1,1) = 1;
    prsVtx(2,2,1,1) = options.prsValQ1RRL;
  end
  if (isfield(options,'prsValQ1RRR'))
    isPrsVtx(2,2,2,1) = 1;
    prsVtx(2,2,2,1) = options.prsValQ1RRR;
  end
  if (isfield(options,'prsValQ2LLL'))
    isPrsVtx(1,1,1,2) = 1;
    prsVtx(1,1,1,2) = options.prsValQ2LLL;
  end
  if (isfield(options,'prsValQ2LLR'))
    isPrsVtx(1,1,2,2) = 1;
    prsVtx(1,1,2,2) = options.prsValQ2LLR;
  end
  if (isfield(options,'prsValQ2LRL'))
    isPrsVtx(1,2,1,2) = 1;
    prsVtx(1,2,1,2) = options.prsValQ2LRL;
  end
  if (isfield(options,'prsValQ2LRR'))
    isPrsVtx(1,2,2,2) = 1;
    prsVtx(1,2,2,2) = options.prsValQ2LRR;
  end
  if (isfield(options,'prsValQ2RLL'))
    isPrsVtx(2,1,1,2) = 1;
    prsVtx(2,1,1,2) = options.prsValQ2RLL;
  end
  if (isfield(options,'prsValQ2RLR'))
    isPrsVtx(2,1,2,2) = 1;
    prsVtx(2,1,2,2) = options.prsValQ2RLR;
  end
  if (isfield(options,'prsValQ2RRL'))
    isPrsVtx(2,2,1,2) = 1;
    prsVtx(2,2,1,2) = options.prsValQ2RRL;
  end
  if (isfield(options,'prsValQ2RRR'))
    isPrsVtx(2,2,2,2) = 1;
    prsVtx(2,2,2,2) = options.prsValQ2RRR;
  end
  
  for q=1:2
    for i=1:2
      if (all(reshape(isPrsVtx(i,:,:,q),[],1)) && isNoFlow(1,i,q)==0)
        isPrsBnd(1,i,q) = 1;
        prsValFace = bilinPrs(reshape(prsVtx(i,:,:,q),2,2),ny,nz);
        prsBndVal(1:ny,1:nz,1,i,q) = prsValFace .* ones(ny,nz);
      end
    end
  end
  
  for q=1:2
    for j=1:2
      if (all(reshape(isPrsVtx(:,j,:,q),[],1)) && isNoFlow(2,j,q)==0)
        isPrsBnd(2,j,q) = 1;
        prsValFace = bilinPrs(reshape(prsVtx(:,j,:,q),2,2),nx,nz);
        prsBndVal(1:nx,1:nz,2,j,q) = prsValFace .* ones(nx,nz);
      end
    end
  end 
  
  for q=1:2
    for k=1:2
      if (all(reshape(isPrsVtx(:,:,k,q),[],1)) && isNoFlow(3,k,q)==0)
        isPrsBnd(3,k,q) = 1;
        prsValFace = bilinPrs(reshape(prsVtx(:,:,k,q),2,2),nx,ny);
        prsBndVal(1:nx,1:ny,3,k,q) = prsValFace .* ones(nx,ny);
      end
    end
  end

end

function [prsVal] = bilinPrs(prsVtx,n,m)

  dn = 0.5/n;
  dm = 0.5/m;
  axn = linspace(dn,1.0-dn,n);
  axm = linspace(dm,1.0-dm,m);
  
  p11 = prsVtx(1,1);
  p21 = prsVtx(2,1);
  p12 = prsVtx(1,2);
  p22 = prsVtx(2,2);
  
  pn1 = p11 + (p21-p11)*axn;
  pnm = p12 + (p22-p12)*axn;
  
  prsVal = pn1'*ones(1,m) + (pnm-pn1)'*axm;
end
