function [state]=runFlowSolver3D2QSort(options)


while (1) %%% enabel goto hack ...

  %Domain
  xL = options.xL;
  yL = options.yL;
  zL = options.zL;
  
  %Discretization
  nx = options.nx;
  ny = options.ny;
  nz = options.nz;
  
  %Time
  time = options.time;
  nTime = length(time);
  dTime = time(2)-time(1);
  
  nface = 6;
  nn = nx*ny*nz;

  dx = xL/nx;
  dy = yL/ny;
  dz = zL/nz;
  
  tol = 1.0e-12;

  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  %%% Pressure and velocities for all compartments
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  srcNull=zeros(0);
  maskNull=zeros(0);
  [pres, velX, velY, velZ, rateQ, termSrcFlux, termSinkFlux] = flowSolver3D2Q(options, srcNull, maskNull);
  
  rateQ = rateQ .* (rateQ > 0);
  
%   presQ1_ = pres(:,:,:,1)
%   presQ2_ = pres(:,:,:,2)
%   rateQ_ = rateQ
%   velZ_ = velZ(:,:,:,1)
  
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  %%% Contrast concentration for comparment one (arterial)
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  %%%Order cells according to upstream priority:
  
  [tosArt] = upstreamSort3D(nx, ny, nz, velX(:,:,:,1), velY(:,:,:,1), velZ(:,:,:,1));
  
  %tosMtx(tos) = 1:nn;
  %tosMtx_ = reshape(tosMtx,nx,ny,nz);
   
  %%%Compute time of flight from influx boundaries:

  %Face-wise fluxes for each voxel, positive into voxel:
  flux = zeros(nx,ny,nz,6);
  flux(1:nx,:,:,1) = dy*dz*velX(1:nx,:,:,1);
  flux(1:nx,:,:,2) = -dy*dz*velX(2:nx+1,:,:,1);
  flux(:,1:ny,:,3) = dx*dz*velY(:,1:ny,:,1);
  flux(:,1:ny,:,4) = -dx*dz*velY(:,2:ny+1,:,1);
  flux(:,:,1:nz,5) = dx*dy*velZ(:,:,1:nz,1);
  flux(:,:,1:nz,6) = -dx*dy*velZ(:,:,2:nz+1,1);   
  flux = reshape(flux,[nn,6]);
  
  %Voxel-wise source terms:
  source = -rateQ*dx*dy*dz;
  source = reshape(source,1,[]);
  
  % External source term for terminal sources
  extSource = zeros(nx,ny,nz);
  if (isfield(options,'termSrc'))
    i=options.termSrc.i;
    j=options.termSrc.j;
    k=options.termSrc.k;
    for n=1:size(i,2)  
      extSource(i(n),j(n),k(n)) = termSrcFlux(n);
    end
  end
  extSource = reshape(extSource,1,[]);

  %Pore-volume 
  if isfield(options,'porosityArt')
    porosityArt = options.porosityArt;
  else
    %Default: decay linearly with distance to BOTTOM (z=0) boundary.
    exy = ones(nx*ny,1);
    zM = linspace(zL-0.5*dz,0.5*dz,nz);
    porosityArt = 0.25*(exy*zM);
  end
  porosityArt = porosityArt .* options.mask(:,:,:,1);
  porosityArt = reshape(porosityArt,nx,ny,nz);
  porVolArt = reshape(dx*dx*dz*porosityArt,1,[]);
  
  % Time of flight
  [tof, dtof] = upstreamTof3D(nx, ny, nz, tosArt,flux, source, porVolArt,extSource);
  tofArt = reshape(tof,nx,ny,nz);
  %dtofArt_ = reshape(dtof,nx,ny,nz);
  
  %tofMaxArt_ = permute(max(tofArt(:,:,:),[],3),[1,2,3])
  
  %if isfield(options,'doCAC')
      % disp("### rundFlowSolver3D2QSort: do tracer calculations ###");
  %else
  %   disp("### rundFlowSolver3D2QSort: skip tracer calculations ###");
  %   break; %%% skip tracer calculations
  %end
  
  %%%Dynamic tracer distribution:
  
  %Tracer boundary conditions:
  if isfield(options,'cValBnd')
    cValBnd = options.cValBnd;
  else
    cValBnd = 1.0;
  end
  if isfield(options,'timeStart')
    timeStart = options.timeStart;
  else
    timeStart = 0.0;
  end
  if isfield(options,'timeEnd')
    timeEnd = min(options.timeEnd,time(nTime));
  else
    timeEnd = time(floor(nTime/7));
  end
  nStart = int64(1+floor((timeStart+tol)/dTime));
  nStop = int64(1+floor((timeEnd+tol)/dTime));
  maxFace = max([nx*ny,nx*nz,ny*nz]);
  if isfield(options,'empiricalAIF') % empiricalAIF is entered as two columns, giving time and AIF(time)
      aifcurve=0.0:dTime:(time(nTime)-timeStart+dTime);
      aifcurve=interp1(options.empiricalAIF(:,1),options.empiricalAIF(:,2),aifcurve);
      tracerBndCnd(:,nStart:nTime,:) = repmat(repmat(aifcurve(1:nTime-nStart+1),[maxFace,1]),[1,1,nface]);
  elseif (isfield(options,'gammaShape') && isfield(options,'gammaScale'))
      a=options.gammaShape;
      b=options.gammaScale;
      const=gamma(a)*b.^a;
      gamtab=0.0:dTime:(time(nTime)-timeStart+dTime);
      gamtab=cValBnd*const*exp(-b*gamtab).*gamtab.^(a-1);
      tracerBndCnd = zeros(maxFace,nTime,nface);
      %     size(tracerBndCnd(:,nStart:nTime,:))
      %     size(gamtab(1:nTime-nStart+1))
      %     size(repmat(repmat(gamtab(1:nTime-nStart+1),[maxFace,1]),[1,1,nface]))
      tracerBndCnd(:,nStart:nTime,:) = repmat(repmat(gamtab(1:nTime-nStart+1),[maxFace,1]),[1,1,nface]);
  else
      tracerBndCnd = zeros(max([nx*ny,nx*nz,ny*nz]),nTime,nface);
      tracerBndCnd(:,nStart:nStop,:) = cValBnd;
      %tracerBndBottom = zeros(nx,ny,nTime);
      %tracerBndBottom(3,2,11:50) = 1.0;
      %tracerBndBottom(1:nx,1:ny,floor(nTime/10):floor(nTime/7)) = 1.0;
      %tracerBndCnd(1:nx*ny,:,5) = reshape(tracerBndBottom,nx*ny,[]);
  end
  %tracerBndCnd_ = permute(tracerBndCnd(1,:,:),[2,3,1])
  
  %Tracer source
  tracerSourceArt = zeros(nn,nTime);
  termTracerProfile = [];
  if (isfield(options,'termSrc') && size(termSrcFlux,1) == size(options.termSrc.i,2))
      sz=size(termSrcFlux,1);
      if isfield(options,'empiricalAIF')
          termTracerProfile = repmat(aifcurve(1:nTime-nStart+1),sz,1);
      else
          assert(size(gamtab,2) >= nTime-nStart+1, 'gamtab not available')
          termTracerProfile = repmat(gamtab(1:nTime-nStart+1),sz,1);
      end
  end
  tracerSourceArt = reshape(tracerSourceArt,nn,[]);
  
%   tofMin=min(tof(:))
%   dtofMin=min(dtof(:))

  mask=reshape(options.mask(:,:,:,1),1,[]);
  [upTracerArt, dwnTracerArt, volTracerArt] = tracerDynForward3D(nx, ny, nz, tosArt,flux,source,tracerBndCnd,tracerSourceArt,time,tof,dtof,1.0,porVolArt,mask,options.termSrc,termSrcFlux,termTracerProfile);
% %   volTracerArt_ = reshape(volTracerArt,nx,ny,nz,[]);
% %   volTracerArt__ = permute(volTracerArt_(1,1,1:nz,:),[4,3,1,2]);
% %   volTracerMaxArt_ = max(max(volTracerArt_,[],4),[],3)
% %   dwnTracerArt_ = reshape(dwnTracerArt,nx,ny,nz,[]);
% %   dwnTracerArt_ = permute(dwnTracerArt_(1,1,1:nz,:),[4,3,1,2]);
% %   
% %   upTracerArt_ = reshape(upTracerArt,nx,ny,nz,[]);
% %   upracerArt_ = permute(upTracerArt_(1,1,1:nz,:),[4,3,1,2]);
  
  totCACArt = reshape(porVolArt,1,[]) * volTracerArt;
  
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  %%% Contrast concentration for capillary bed:
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  if isfield(options,'porosityQ')
    porosityQ = options.porosityQ;
  else
    porosityQ = 0.05*ones(nx,ny,nz);
  end
  poreVolumeQ = dx*dy*dz*porosityQ;
  volumeRateQ = rateQ*dx*dy*dz;
  transitTimeQ = (options.maskQ .* poreVolumeQ)./(volumeRateQ + 1.0e-18); % or simply porosity/rate
  
  tofCap = (tofArt + transitTimeQ) .* options.maskQ;
  
  %tofMaxCap_= permute(max(tofCap(:,:,:),[],3),[1,2,3])
  
  tosCap = 1:nn;
  fluxCap = zeros(nx*ny*nz,6);
  sourceCap = -source;
  tracerBndCndCap = zeros(max([nx*ny,nx*nz,ny*nz]),nTime,nface);
  tracerSourceCap = volTracerArt; % Er nå i orden??? NB: Denne må fikses når når vxl-konsentrasjoner er beregnet ...
  %tracerSourceCap = 0.5*(upTracerArt+dwnTracerArt);
  tofCap0 = reshape(tofCap,nn,1);
  dtofCap = reshape(transitTimeQ,nn,1);
  
  termSrcDummy = [];
  termTracerProfile = [];
  
  mask=reshape(options.maskQ(:,:,:),1,[]);
  [upTracerCap, dwnTracerCap, volTracerCap] = tracerDynForward3D(nx, ny, nz, tosCap,fluxCap,sourceCap,tracerBndCndCap,tracerSourceCap,time,tofCap0,dtofCap,-1.0,poreVolumeQ,mask,options.termSrc,termSrcDummy,termTracerProfile);
  
% %   volTracerCap_ = reshape(volTracerCap,nx,ny,nz,[]);
% %   volTracerCap__ = permute(volTracerCap_(1,1,1:nz,:),[4,3,1,2]);
% %   volTracerMaxCap_ = max(max(volTracerCap_,[],4),[],3);
  
  totCACCap = reshape(poreVolumeQ,1,[]) * volTracerCap;

  %Set boundary for venular concentration for debugging purposes
  %   for kkk=1:nn 
  %     dwnTracerCap(kkk,:) = tracerBndCnd(1,:,1);
  %   end

  volTracerCapOut = dwnTracerCap;

  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  %%% Compartment two (venular)
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  %Order cells according to upstream priority:
  [tosVen] = upstreamSort3D(nx, ny, nz, velX(:,:,:,2), velY(:,:,:,2), velZ(:,:,:,2));
  
  %%tosVenMtx(tosVen) = 1:nn;
  
  %Pore-volume, similarly to art:
  if isfield(options,'porosityVen')
    porosityVen = options.porosityVen;
  else
    porosityVen = 2*porosityArt;
  end
  porosityVen = porosityVen .* options.mask(:,:,:,2);
  porVolVen = porosityVen*dx*dy*dz;
  
  %Face-wise fluxes for each voxel, positive into voxel, ...
  fluxVen = zeros(nx,ny,nz,6);
  fluxVen(1:nx,:,:,1) = dy*dz*velX(1:nx,:,:,2);
  fluxVen(1:nx,:,:,2) = -dy*dz*velX(2:nx+1,:,:,2);
  fluxVen(:,1:ny,:,3) = dx*dz*velY(:,1:ny,:,2);
  fluxVen(:,1:ny,:,4) = -dx*dz*velY(:,2:ny+1,:,2);
  fluxVen(:,:,1:nz,5) = dx*dy*velZ(:,:,1:nz,2);
  fluxVen(:,:,1:nz,6) = -dx*dy*velZ(:,:,2:nz+1,2);   
  fluxVen = reshape(fluxVen,[nn,6]);
  
% %   %Outflux at bottom:
% %   fluxVen(:,:,1,5) = dx*dy*velZ(:,:,1,2);
% %   fluxVen = reshape(fluxVen,[nn,6]);
  
  % Time of flight
  sourceVen = zeros(nx,ny,nz);
  sourceVen = sourceVen + tofCap .* rateQ*dx*dy*dz;
  sourceVen = reshape(sourceVen,1,[]);
  
  %External source term
  extSource = -(sum(flux,2)'+reshape(rateQ*dx*dy*dz,1,[]));
  
  % External source term for terminal sinks
  extSource = zeros(nx,ny,nz);
  if (isfield(options,'termSink'))
    i=options.termSink.i;
    j=options.termSink.j;
    k=options.termSink.k;
    for n=1:size(i,2)  
      extSource(i(n),j(n),k(n)) = termSinkFlux(n);
    end
  end
  extSource = reshape(extSource,1,[]);
  
  [tofVen, dtofVen] = upstreamTof3D(nx, ny, nz, tosVen, fluxVen, sourceVen, porVolVen, extSource);
  
  % Dynamic trace distribution
  tracerBndCndVen = zeros(max([nx*ny,nx*nz,ny*nz]),nTime,nface);
  sourceVen = rateQ*dx*dy*dz;
  tracerSourceVen = dwnTracerCap; % Er denne rett? Tracer UT av cap ... Ser ok ut.
% %   tracerSourceVen = dwnTracerCap;
% %   tracerSourceVen(:,1:nTime-1) = 0.5*(dwnTracerCap(:,2:nTime)+dwnTracerCap(:,1:nTime-1));
  termTracerProfile = [];
  
  mask=reshape(options.mask(:,:,:,2),1,[]);
  [upTracerVen, dwnTracerVen, volTracerVen] = tracerDynForward3D(nx,ny,nz, tosVen,fluxVen,sourceVen,tracerBndCndVen,tracerSourceVen,time,tofVen,dtofVen,1.0,porVolVen,mask,options.termSink,termSinkFlux,termTracerProfile);
% %   volTracerVen_ = reshape(volTracerVen,nx,ny,nz,[]);
% %   volTracerVen__ = permute(volTracerVen_(1,1,1:nz,:),[4,3,1,2]);
% %   volTracerMaxVen_ = max(max(volTracerVen_,[],4),[],3);
  
  totCACVen = reshape(porVolVen,1,[]) * volTracerVen;  
  dwnCACVen = reshape(porVolVen,1,[]) * dwnTracerVen;  


%%% Arterial CAC
% % % % % %   figure('Name','Arterial concentration, T=10','NumberTitle','off')
% % % % % %   %colormap(jet(128));  
% % % % % %   cMtx= reshape(volTracerArt,nx,ny,nz,[]);
% % % % % %   %inc = floor(size(cMtx,4)/60);
% % % % % %   jj = 10;
% % % % % %   climsVen = [0 max(volTracerArt(:,jj))];
% % % % % %   for j=1:10
% % % % % %     subplot(4,4,j);
% % % % % %     imagesc(cMtx(:,:,j,jj),climsVen);
% % % % % %     set(gca,'YDir','normal')
% % % % % %     %colormap(hot)
% % % % % %   end
% % % % % %   subplot(4,4,16);
% % % % % %   colorbar;
% % % % % %   
% % % % % %   %%% CapBed CAC
% % % % % %   figure('Name','CapBed concentration, T=25','NumberTitle','off')
% % % % % %   %colormap(jet(128));  
% % % % % %   cMtx= reshape(volTracerCap,nx,ny,nz,[]);
% % % % % %   %inc = floor(size(cMtx,4)/60);
% % % % % %   jj = 25;
% % % % % %   climsVen = [0 max(volTracerCap(:,jj))];
% % % % % %   for j=1:10
% % % % % %     subplot(4,4,j);
% % % % % %     imagesc(cMtx(:,:,j,jj),climsVen);
% % % % % %     set(gca,'YDir','normal')
% % % % % %     %colormap(hot)
% % % % % %   end
% % % % % %   subplot(4,4,16);
% % % % % %   colorbar;
% % % % % %   
% % % % % %   %%% Venular CAC
% % % % % %   figure('Name','Venular concentration, T=40','NumberTitle','off')
% % % % % %   %colormap(jet(128));  
% % % % % %   cMtx= reshape(volTracerVen,nx,ny,nz,[]);
% % % % % %   %inc = floor(size(cMtx,4)/60);
% % % % % %   jj = 40;
% % % % % %   climsVen = [0 max(volTracerVen(:,jj))];
% % % % % %   for j=1:10
% % % % % %     subplot(4,4,j);
% % % % % %     imagesc(cMtx(:,:,j,jj),climsVen);
% % % % % %     set(gca,'YDir','normal')
% % % % % %     %colormap(hot)
% % % % % %   end
% % % % % %   subplot(4,4,16);
% % % % % %   colorbar;
% % % % % %   
% % % % % %   %%% Total CAC layer 7
% % % % % %   figure('Name','Total concentration, layer 7','NumberTitle','off')
% % % % % %   %colormap(jet(128));  
% % % % % %   cMtx= reshape(volTracerVen+volTracerArt+volTracerCap,nx,ny,nz,[]);
% % % % % %   %inc = floor(size(cMtx,4)/60);
% % % % % %   jj = 7;
% % % % % %   climsVen = [0 max(reshape(cMtx,1,[]))];
% % % % % %   for j=1:15
% % % % % %     subplot(4,4,j);
% % % % % %     imagesc(cMtx(:,:,jj,j*4),climsVen);
% % % % % %     set(gca,'YDir','normal')
% % % % % %     %colormap(hot)
% % % % % %   end
% % % % % %   subplot(4,4,16);
% % % % % %   colorbar;
  
  %%% Deconvolution:
%   volTracerTot=volTracerArt(62*66*5+31.5*66,:)+volTracerCap(62*66*5+31.5*66,:)+volTracerVen(62*66*5+31.5*66,:);
%   col=volTracerArt(2162,:);
%   row=zeros(1,901);
%   row(1,1)=col(1,1);
%   mtxA=toeplitz(col,row);
%   [U,S,V]=svd(mtxA);
%   figure
%   hold on
%   plot(V(:,1))
%   plot(V(:,5))
%   plot(V(:,10))
%   hold off
%   figure
%   hold on
%   coeff=volTracerTot*U(:,:);
%   diagS=diag(S);
%   coeffSinv=coeff'./diagS(:);
%   plot(V(:,1:200)*coeffSinv(1:200))
%   legend
%   plot(V(:,1:400)*coeffSinv(1:400))
%   plot(V(:,1:600)*coeffSinv(1:600))
%   plot(V(:,1:800)*coeffSinv(1:800))
%   hold off
%   figure
  
 totCACcnv = totCACArt + totCACCap + totCACVen;
 
 
 %  figure('Name','aif')%,'NumberTitle','off')
 % plot(time,gamtab(1:end-1),'k','LineWidth',2);
 % hold on
 % ax = gca;
 % ax.XMinorTick = 'on';
 % hold off
  
 % figure('Name','MRIconv: Total contrast agent vs time for compartments')%,'NumberTitle','off')
 % plot(time,totCACArt,'b','LineWidth',2);
 % hold on
 % plot(time,totCACCap,'r','LineWidth',2);
 % plot(time,totCACVen,'g','LineWidth',2);
 % plot(time,totCACcnv,'m--','LineWidth',2);
 % ax = gca;
 % ax.XMinorTick = 'on';
 % hold off
%   
%   figure('Name','MRIconv: Contrast agent vs time for selected vxl')%,'NumberTitle','off')
%   vxl = int64(nx/2*ny*nz/2);
%   plot(time,volTracerArt(vxl,:),'b','LineWidth',2);
%   hold on
%   plot(time,volTracerCap(vxl,:),'r','LineWidth',2);
%   plot(time,volTracerVen(vxl,:),'g','LineWidth',2);
%   volTracerTotal = volTracerArt(vxl,:)+volTracerCap(vxl,:)+volTracerVen(vxl,:);
%   plot(time,volTracerTotal,'m--','LineWidth',2);
%   ax = gca;
%   ax.XMinorTick = 'on';
%   hold off
%   
%   figure('Name','Arterial concentration','NumberTitle','off')
%   %climsVen = [0 1.0];
%   %colormap(jet(128));  
%   cMtx= reshape(volTracerArt,nx,ny,nz,[]);
%   inc = floor(size(cMtx,4)/60);
%   jj = 1;
%   for j=1:60
%     subplot(6,10,j);
%     imagesc(cMtx(:,:,1,jj)'); %,climsVen)
%     jj = jj + inc;
%     set(gca,'YDir','normal')
%     %colormap(hot)
%   end

%   figure('Name','CapBedUp concentration','NumberTitle','off')
%   %climsVen = [0 1.0];
%   %colormap(jet(128));  
%   cMtx= reshape(upTracerCap,nx,ny,nz,[]);
%   inc = floor(size(cMtx,4)/60);
%   jj = 1;
%   for j=1:60
%     subplot(6,10,j);
%     imagesc(cMtx(:,:,1,jj)'); %,climsVen)
%     jj = jj + inc;
%     set(gca,'YDir','normal')
%     %colormap(hot)
%   end
% 
%   figure('Name','CapBedDwn concentration','NumberTitle','off')
%   %climsVen = [0 1.0];
%   %colormap(jet(128));  
%   cMtx= reshape(dwnTracerCap,nx,ny,nz,[]);
%   inc = floor(size(cMtx,4)/60);
%   jj = 1;
%   for j=1:60
%     subplot(6,10,j);
%     imagesc(cMtx(:,:,1,jj)'); %,climsVen)
%     jj = jj + inc;
%     set(gca,'YDir','normal')
%     %colormap(hot)
%   end
% 
%   figure('Name','Venular concentration','NumberTitle','off')
%   %climsVen = [0 1.0];
%   %colormap(jet(128));  
%   cMtx= reshape(volTracerVen,nx,ny,nz,[]);
%   inc = floor(size(cMtx,4)/60);
%   jj = 1;
%   for j=1:60
%     subplot(6,10,j);
%     imagesc(cMtx(:,:,1,jj)'); %,climsVen)
%     jj = jj + inc;
%     set(gca,'YDir','normal')
%     %colormap(hot)
%   end

  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  %%%%%%%%%%%%  END  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  %save('conv.mat', 'volTracerArt', 'volTracerCap', 'volTracerCapOut', 'volTracerVen', 'totCACcnv');
  
  state.volTracerArt = volTracerArt;
  state.volTracerCap = volTracerCap;
  state.volTracerVen = volTracerVen;
  
  break;
  
end %%% while loop to enable goto-hack (skip tracer calculations)
  
  state.velX = velX;
  state.velY = velY;
  state.velZ = velZ;
  state.rateQ = rateQ;
  state.pres = pres;
  
  state.termSinkFlux = termSinkFlux;
  state.termSrcFlux = termSrcFlux;
  
  if isfield(options,'doCAC')
    state.tosArt = tosArt;
    state.tosCap = tosCap;
    state.tosVen = tosVen;
    state.fluxVen = fluxVen;
    state.tofVen = tofVen;
    state.dtofVen = dtofVen;
  end
  state.tof = tof;
  state.dtof = dtof;
  state.flux = flux;
  state.source = source;
  state.porosityArt = porosityArt;
  %state.porosityQ = porosityQ;
  %state.tofCap = tofCap;
  state.tofArt = tofArt;
    state.tofVen = tofVen;
    state.dtofVen = dtofVen;
    state.tofCap = tofCap;
    state.dtofCap = dtofCap;
    vol=xL*yL*zL;
    state.totCACcnv = totCACcnv/vol;
    state.totCACArt = totCACArt/vol;
    state.totCACCap = totCACCap/vol;
    state.totCACVen = totCACVen/vol;
    state.dwnCACVen = dwnCACVen/vol;
     
  return
  
end
