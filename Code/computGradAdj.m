function [gradPerm, gradPor, sampleMask, lamTau, lamSrc, lamPrs]=computGradAdj(options, state, tofObserved, alphaPermReg, vb)

  
  gradPerm=0;
  gradPor=0;
  
  nx = options.nx;
  ny = options.ny;
  nz = options.nz;

  dx = options.xL/nx;
  dy = options.yL/ny;
  dz = options.zL/nz;
  
  dV=dx*dy*dz;
  porvol=state.porosityArt*dV;
   
  nn = nx*ny*nz;
  
  tol = 1.0e-20;
  indexNatural = 1:nn;
  [I,J,K] = ind2sub([nx ny nz],indexNatural);
  idNb = zeros(nn,6);
  idNb(:,1) = (indexNatural' - 1) .* (I(:) > 1);
  idNb(:,2) = (indexNatural' + 1) .* (I(:) < nx);
  idNb(:,3) = (indexNatural' - nx) .* (J(:) > 1);
  idNb(:,4) = (indexNatural' + nx) .* (J(:) < ny);
  idNb(:,5) = (indexNatural' - nx*ny) .* (K(:) > 1);
  idNb(:,6) = (indexNatural' + nx*ny) .* (K(:) < nz);

  lamTau=zeros(nn,1);
  
  %sampleMask=(I(:) > 4).*(I(:) < 7) .* (J(:) > 4).*(J(:) < 7) .* (K(:) >= 1).*(K(:) <= 10);
  %sampleMask=(I(:) > 4).*(I(:) < 7) .* (J(:) > 4).*(J(:) < 7) .* (K(:) >= 1).*(K(:) <= 5);
  %sampleMask=(I(:) >= 4).*(I(:) <= 7) .* (J(:) >= 4).*(J(:) <= 7) .* (K(:) >= 2).*(K(:) <= 10);
  %sampleMask=(I(:) >= 4).*(I(:) <= 7) .* (J(:) >= 4).*(J(:) <= 7) .* (K(:) >= 1).*(K(:) <= 5);
  
  
  %%sampleMask=(I(:) >= 3).*(I(:) <= 8) .* (J(:) >= 3).*(J(:) <= 8) .* (K(:) >= 2).*(K(:) <= 10);
  %sampleMask=(I(:) >= 2).*(I(:) <= 9) .* (J(:) >= 2).*(J(:) <= 9) .* (K(:) >= 2).*(K(:) <= 9);
  %sampleMask=(I(:) >= 1).*(I(:) <= 10) .* (J(:) >= 1).*(J(:) <= 10) .* (K(:) >= 1).*(K(:) <= 10);
  sampleMask=reshape(options.mask(:,:,1,1),1,[]);
  
  
  %sampleMask=(I(:) > 0).*(I(:) < 11).*(J(:) > 0).*(J(:) < 11).*(K(:) < 11);
  %sampleMask=(I(:) > 0).*(I(:) < 11).*(J(:) > 0).*(J(:) < 11).*(K(:) < 11);
  %sampleMask=(I(:) > 1).*(I(:) < 10).*(J(:) > 1).*(J(:) < 10).*(K(:) < 11);
  
  %sampleMask=(I(:) > 2).*(I(:) < 9).*(J(:) > 2).*(J(:) < 9).*(K(:) < 11);
  %sampleMask=(I(:) > 3).*(I(:) < 8).*(J(:) > 3).*(J(:) < 8).*(K(:) < 9);
  %sampleMask=(I(:) > 3).*(I(:) < 8).*(J(:) > 3).*(J(:) < 8).*(K(:) < 11);
  %reshape(mask,nx,ny,[])
 
  %---------------------------------------------
  %%% d/dT
  for ii = nn:-1:1
    vxl = state.tos(ii);
    
    if sampleMask(vxl)==1
        i0=I(vxl);
        j0=J(vxl);
        k0=K(vxl);
    end
    
    accTerm = 0;
    factor = 0;
    
    dTau = state.tof(vxl) - tofObserved(vxl);
    accTerm = -dTau*sampleMask(vxl);
    if (state.source(vxl) < - tol)
      factor = -state.source(vxl);
    end
    
    for face=1:6  % flux negativ out of vxl
      if (state.flux(vxl,face) < -tol) && (idNb(vxl,face) > 0)
          accTerm = accTerm - state.flux(vxl,face)*lamTau(idNb(vxl,face));
          factor = factor - state.flux(vxl,face);
      end
    end
    
    if (factor > tol)
      lamTau(vxl) = accTerm/factor;
    end
    
  end
  
  %lamTauMask=reshape((lamTau ~= 0),nx,ny,nz)
  if (vb)
    disp('dtau:');
    (state.tofArt - tofObserved) .* reshape(sampleMask,nx,ny,nz)
    lamTauArt = reshape(lamTau,nx,ny,nz)
  end
  
  %---------------------------------------------
  %%% d/dp
  
  mobX = options.permX./options.visc;
  mobY = options.permY./options.visc;
  mobZ = options.permZ./options.visc;
  mobQ = 2.0*options.permQ./(options.visc(:,:,:,1)+options.visc(:,:,:,2));
  
  xCon = zeros(nx+1,ny,nz,2);
  yCon = zeros(nx,ny+1,nz,2);
  zCon = zeros(nx,ny,nz+1,2);
  qCon = zeros(nx,ny,nz);
  %Interior
  xCon(2:nx,:,:,:) = (dy*dz/dx)*2*mobX(1:nx-1,:,:,:).*mobX(2:nx,:,:,:)./(mobX(1:nx-1,:,:,:)+mobX(2:nx,:,:,:));
  yCon(:,2:ny,:,:) = (dx*dz/dy)*2*mobY(:,1:ny-1,:,:).*mobY(:,2:ny,:,:)./(mobY(:,1:ny-1,:,:)+mobY(:,2:ny,:,:));
  zCon(:,:,2:nz,:) = (dx*dy/dz)*2*mobZ(:,:,1:nz-1,:).*mobZ(:,:,2:nz,:)./(mobZ(:,:,1:nz-1,:)+mobZ(:,:,2:nz,:));
  %Boundary
%   xCon(1,:,:,:) = (dy*dz/dx)*2*mobX(1,:,:,:);
%   xCon(nx+1,:,:,:) = (dy*dz/dx)*2*mobX(nx,:,:,:);
%   yCon(:,1,:,:) = (dx*dz/dy)*2*mobY(:,1,:,:);
%   yCon(:,ny+1,:,:) = (dx*dz/dy)*2*mobY(:,ny,:,:);
%   zCon(:,:,1,:) = (dx*dy/dz)*2*mobZ(:,:,1,:);
%   zCon(:,:,nz+1,:) = (dx*dy/dz)*2*mobZ(:,:,nz,:);
  %qCon(:,:) = 0.5*(dy/dx+dx/dy)*mobQ;
  qCon(:,:,:) = dy*dx*dz*mobQ;
  
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
  
  nbCon = zeros(nx,ny,nz,6);
  nbCon(1:nx,:,:,1) = xCon(1:nx,:,:,1);
  nbCon(1:nx,:,:,2) = xCon(2:nx+1,:,:,1);
  nbCon(:,1:ny,:,3) = yCon(:,1:ny,:,1);
  nbCon(:,1:ny,:,4) = yCon(:,2:ny+1,:,1);
  nbCon(:,:,1:nz,5) = zCon(:,:,1:nz,1);
  nbCon(:,:,1:nz,6) = zCon(:,:,2:nz+1,1);   
  nbCon = reshape(nbCon,[nn,6]);
  
  srcTerm = zeros(nn,2);
  for vxl=1:nn
    srcTerm(vxl,1) =  lamTau(vxl)*state.tof(vxl)*qCon(vxl);
    srcTerm(vxl,2) = -lamTau(vxl)*state.tof(vxl)*qCon(vxl);
    for nFace=1:6
      if (idNb(vxl,nFace) > 0)
          if state.flux(vxl,nFace) > 0 % flux positive into vxl
            srcTerm(vxl,1) = srcTerm(vxl,1) + state.tof(idNb(vxl,nFace))*(lamTau(vxl)-lamTau(idNb(vxl,nFace)))*nbCon(vxl,nFace);
          elseif state.flux(vxl,nFace) < 0 % flux negative out of vxl
            srcTerm(vxl,1) = srcTerm(vxl,1) + state.tof(vxl)*(lamTau(vxl)-lamTau(idNb(vxl,nFace)))*nbCon(vxl,nFace);
          end
      end
    end
  end
  
  %srcTerm=srcTerm.*sampleMask;
  %reshape(srcTerm,nx,ny,nz);
  %maskNull=zeros(0);
  
%   srcTerm(:,1) = srcTerm(:,1) .* (lamTau ~= 0);
%   srcTerm(:,2) = srcTerm(:,2) .* (lamTau ~= 0);
  
  if (vb)
    lamSrc=-reshape(srcTerm,nx,ny,nz,2)
  else
    lamSrc=-reshape(srcTerm,nx,ny,nz,2);
  end
  
  [lamPrs, ~, ~, ~, ~] = flowSolver3D2Q(options, -srcTerm, sampleMask);
  
   
  %lamPrs(:,:,:,1) = lamPrs(:,:,:,1).*reshape(sampleMask,nx,ny,nz);
  if (vb)
    lamPrsArt=lamPrs(:,:,:,1)
  end
%   lamPrsVen=lamPrs(:,:,:,2)
  
  gradPerm = zeros(nx,ny,nz,2,3);
  
  dxRight = zeros(nx,ny,nz,2);
  dxRight(1:nx-1,:,:,:) = (dy*dz/dx)*2*mobX(2:nx,:,:,:).*mobX(2:nx,:,:,:)./(mobX(1:nx-1,:,:,:)+mobX(2:nx,:,:,:))./(mobX(1:nx-1,:,:,:)+mobX(2:nx,:,:,:))./options.visc(1:nx-1,:,:,:);
  dxLeft = zeros(nx,ny,nz,2);
  dxLeft(2:nx,:,:,:) = (dy*dz/dx)*2*mobX(1:nx-1,:,:,:).*mobX(1:nx-1,:,:,:)./(mobX(1:nx-1,:,:,:)+mobX(2:nx,:,:,:))./(mobX(1:nx-1,:,:,:)+mobX(2:nx,:,:,:))./options.visc(2:nx,:,:,:);
  
  dyRight = zeros(nx,ny,nz,2);
  dyRight(:,1:ny-1,:,:) = (dx*dz/dy)*2*mobY(:,2:ny,:,:).*mobY(:,2:ny,:,:)./(mobY(:,1:ny-1,:,:)+mobY(:,2:ny,:,:))./(mobY(:,1:ny-1,:,:)+mobY(:,2:ny,:,:))./options.visc(:,1:ny-1,:,:);
  dyLeft = zeros(nx,ny,nz,2);
  dyLeft(:,2:ny,:,:) = (dx*dz/dy)*2*mobY(:,1:ny-1,:,:).*mobY(:,1:ny-1,:,:)./(mobY(:,1:ny-1,:,:)+mobY(:,2:ny,:,:))./(mobY(:,1:ny-1,:,:)+mobY(:,2:ny,:,:))./options.visc(:,2:ny,:,:);
  
  dzRight = zeros(nx,ny,nz,2);
  dzRight(:,:,1:nz-1,:) = (dx*dy/dz)*2*mobZ(:,:,2:nz,:).*mobZ(:,:,2:nz,:)./(mobZ(:,:,1:nz-1,:)+mobZ(:,:,2:nz,:))./(mobZ(:,:,1:nz-1,:)+mobZ(:,:,2:nz,:))./options.visc(:,:,1:nz-1,:);
  dzLeft = zeros(nx,ny,nz,2);
  dzLeft(:,:,2:nz,:) = (dx*dy/dz)*2*mobZ(:,:,1:nz-1,:).*mobZ(:,:,1:nz-1,:)./(mobZ(:,:,1:nz-1,:)+mobZ(:,:,2:nz,:))./(mobZ(:,:,1:nz-1,:)+mobZ(:,:,2:nz,:))./options.visc(:,:,2:nz,:);
  
  if (isfield(options,'mask'))
    mask = options.mask;
    dxLeft(2:nx,:,:,:) = dxLeft(2:nx,:,:,:).*mask(1:nx-1,:,:,:).*mask(2:nx,:,:,:);
    dxRight(1:nx-1,:,:,:) = dxRight(1:nx-1,:,:,:).*mask(1:nx-1,:,:,:).*mask(2:nx,:,:,:);
    
    dyLeft(:,2:ny,:,:) = dyLeft(:,2:ny,:,:).*mask(:,1:ny-1,:,:).*mask(:,2:ny,:,:);
    dyRight(:,1:ny-1,:,:) = dyRight(:,1:ny-1,:,:).*mask(:,1:ny-1,:,:).*mask(:,2:ny,:,:);
    
    dzLeft(:,:,2:nz,:) = dzLeft(:,:,2:nz,:).*mask(:,:,1:nz-1,:).*mask(:,:,2:nz,:);
    dzRight(:,:,1:nz-1,:) = dzRight(:,:,1:nz-1,:).*mask(:,:,1:nz-1,:).*mask(:,:,2:nz,:);
  end
  
% %   sm = reshape(sampleMask,nx,ny,nz);
% %   dxLeft(2:nx,:,:,:) = dxLeft(2:nx,:,:,:).*sm(1:nx-1,:,:,:).*sm(2:nx,:,:,:);
% %   dxRight(1:nx-1,:,:,:) = dxRight(1:nx-1,:,:,:).*sm(1:nx-1,:,:,:).*sm(2:nx,:,:,:);
% % 
% %   dyLeft(:,2:ny,:,:) = dyLeft(:,2:ny,:,:).*sm(:,1:ny-1,:,:).*sm(:,2:ny,:,:);
% %   dyRight(:,1:ny-1,:,:) = dyRight(:,1:ny-1,:,:).*sm(:,1:ny-1,:,:).*sm(:,2:ny,:,:);
% % 
% %   dzLeft(:,:,2:nz,:) = dzLeft(:,:,2:nz,:).*sm(:,:,1:nz-1,:).*sm(:,:,2:nz,:);
% %   dzRight(:,:,1:nz-1,:) = dzRight(:,:,1:nz-1,:).*sm(:,:,1:nz-1,:).*sm(:,:,2:nz,:);
  
  gradPerm(1:nx-1,:,:,:,1) = (lamPrs(1:nx-1,:,:,:)-lamPrs(2:nx,:,:,:)).*dxRight(1:nx-1,:,:,:).*(state.pres(1:nx-1,:,:,:)-state.pres(2:nx,:,:,:));
  gradPerm(2:nx,:,:,:,1) = gradPerm(2:nx,:,:,:,1) + (lamPrs(1:nx-1,:,:,:)-lamPrs(2:nx,:,:,:)).*dxLeft(2:nx,:,:,:).*(state.pres(1:nx-1,:,:,:)-state.pres(2:nx,:,:,:));
  
  gradPerm(:,1:ny-1,:,:,2) = (lamPrs(:,1:ny-1,:,:)-lamPrs(:,2:ny,:,:)).*dyRight(:,1:ny-1,:,:).*(state.pres(:,1:ny-1,:,:)-state.pres(:,2:ny,:,:));
  gradPerm(:,2:ny,:,:,2) = gradPerm(:,2:ny,:,:,2) + (lamPrs(:,1:ny-1,:,:)-lamPrs(:,2:ny,:,:)).*dyLeft(:,2:ny,:,:).*(state.pres(:,1:ny-1,:,:)-state.pres(:,2:ny,:,:));
  
  gradPerm(:,:,1:nz-1,:,3) = (lamPrs(:,:,1:nz-1,:)-lamPrs(:,:,2:nz,:)).*dzRight(:,:,1:nz-1,:).*(state.pres(:,:,1:nz-1,:)-state.pres(:,:,2:nz,:));
  gradPerm(:,:,2:nz,:,3) = gradPerm(:,:,2:nz,:,3) + (lamPrs(:,:,1:nz-1,:)-lamPrs(:,:,2:nz,:)).*dzLeft(:,:,2:nz,:).*(state.pres(:,:,1:nz-1,:)-state.pres(:,:,2:nz,:));
  
  % missing prs boundaries ...
  
  tauxLeft=zeros(nx,ny,nz);
  tauxRight=zeros(nx,ny,nz);
  
  tauyLeft=zeros(nx,ny,nz);
  tauyRight=zeros(nx,ny,nz);
  
  tauzLeft=zeros(nx,ny,nz);
  tauzRight=zeros(nx,ny,nz);
  
  tauCenter=reshape(state.tof,nx,ny,nz);
  lamTau = reshape(lamTau,nx,ny,nz);
  
  tauxLeft(2:nx,:,:) = (state.pres(1:nx-1,:,:,1)>=state.pres(2:nx,:,:,1)).*tauCenter(1:nx-1,:,:) + (state.pres(1:nx-1,:,:,1)<state.pres(2:nx,:,:,1)).*tauCenter(2:nx,:,:);
  tauxRight(1:nx-1,:,:) = (state.pres(1:nx-1,:,:,1)>=state.pres(2:nx,:,:,1)).*tauCenter(1:nx-1,:,:) + (state.pres(1:nx-1,:,:,1)<state.pres(2:nx,:,:,1)).*tauCenter(2:nx,:,:);
  
  tauyLeft(:,2:ny,:) = (state.pres(:,1:ny-1,:,1)>=state.pres(:,2:ny,:,1)).*tauCenter(:,1:ny-1,:) + (state.pres(:,1:ny-1,:,1)<state.pres(:,2:ny,:,1)).*tauCenter(:,2:ny,:);
  tauyRight(:,1:ny-1,:) = (state.pres(:,1:ny-1,:,1)>=state.pres(:,2:ny,:,1)).*tauCenter(:,1:ny-1,:) + (state.pres(:,1:ny-1,:,1)<state.pres(:,2:ny,:,1)).*tauCenter(:,2:ny,:);
  
  tauzLeft(:,:,2:nz) = (state.pres(:,:,1:nz-1,1)>=state.pres(:,:,2:nz,1)).*tauCenter(:,:,1:nz-1) + (state.pres(:,:,1:nz-1,1)<state.pres(:,:,2:nz,1)).*tauCenter(:,:,2:nz);
  tauzRight(:,:,1:nz-1) = (state.pres(:,:,1:nz-1,1)>=state.pres(:,:,2:nz,1)).*tauCenter(:,:,1:nz-1,:) + (state.pres(:,:,1:nz-1,1)<state.pres(:,:,2:nz,1)).*tauCenter(:,:,2:nz);
  
  gradPerm(1:nx-1,:,:,1,1) = gradPerm(1:nx-1,:,:,1,1) + (lamTau(1:nx-1,:,:)-lamTau(2:nx,:,:)).*tauxRight(1:nx-1,:,:).*dxRight(1:nx-1,:,:,1).*(state.pres(1:nx-1,:,:,1)-state.pres(2:nx,:,:,1));
  gradPerm(2:nx,:,:,1,1) = gradPerm(2:nx,:,:,1,1) + (lamTau(1:nx-1,:,:)-lamTau(2:nx,:,:)).*tauxLeft(2:nx,:,:).*dxLeft(2:nx,:,:,1).*(state.pres(1:nx-1,:,:,1)-state.pres(2:nx,:,:,1));
  
  gradPerm(:,1:ny-1,:,1,2) = gradPerm(:,1:ny-1,:,1,2) + (lamTau(:,1:ny-1,:)-lamTau(:,2:ny,:)).*tauyRight(:,1:ny-1,:).*dyRight(:,1:ny-1,:,1).*(state.pres(:,1:ny-1,:,1)-state.pres(:,2:ny,:,1));
  gradPerm(:,2:ny,:,1,2) = gradPerm(:,2:ny,:,1,2) + (lamTau(:,1:ny-1,:)-lamTau(:,2:ny,:)).*tauyLeft(:,2:ny,:).*dyLeft(:,2:ny,:,1).*(state.pres(:,1:ny-1,:,1)-state.pres(:,2:ny,:,1));
  
  gradPerm(:,:,1:nz-1,1,3) = gradPerm(:,:,1:nz-1,1,3) + (lamTau(:,:,1:nz-1)-lamTau(:,:,2:nz)).*tauzRight(:,:,1:nz-1).*dzRight(:,:,1:nz-1,1).*(state.pres(:,:,1:nz-1,1)-state.pres(:,:,2:nz,1));
  gradPerm(:,:,2:nz,1,3) = gradPerm(:,:,2:nz,1,3) + (lamTau(:,:,1:nz-1)-lamTau(:,:,2:nz)).*tauzLeft(:,:,2:nz).*dzLeft(:,:,2:nz,1).*(state.pres(:,:,1:nz-1,1)-state.pres(:,:,2:nz,1));
  
  %regularization terms
  gradPerm(:,:,:,1,1) = gradPerm(:,:,:,1,1) + alphaPermReg*(options.permX(:,:,:,1)-2e-10);
  gradPerm(:,:,:,1,2) = gradPerm(:,:,:,1,2) + alphaPermReg*(options.permY(:,:,:,1)-2e-10);
  gradPerm(:,:,:,1,3) = gradPerm(:,:,:,1,3) + alphaPermReg*options.permZ(:,:,:,1);
  
%   gradPerm(:,:,:,1,1) = reshape(gradPerm(:,:,:,1,1),nx,ny,nz) .* (lamTau ~= 0);
%   gradPerm(:,:,:,1,2) = reshape(gradPerm(:,:,:,1,2),nx,ny,nz) .* (lamTau ~= 0);
%   gradPerm(:,:,:,1,3) = reshape(gradPerm(:,:,:,1,3),nx,ny,nz) .* (lamTau ~= 0);
  
%   dfdk1Art=gradPerm(:,:,1:3,1,1)
%   dfdk2Art=gradPerm(:,:,1:3,1,2)
%   dfdk3Art=gradPerm(:,:,1:10,1,3)
%   dfdk1Ven=gradPerm(:,:,:,2,1)
  
  
  gradPor = -dx*dy*dx*lamTau;
%   dfdPorArt=reshape(gradPor,nx,ny,nz)
  
  return
end