%Computes dynamic tracer distribution from abitrary set of source(s).
%
%%%%Input:
%
% nx: Number of grid cells in x-direction
% ny: Number of grid cslls in y-direction
% nz: Number of grid cslls in z-direction
% tos: Voxel indices sorted according to upstream-before-downstream.
% flux: Array of voxel-wise, face-wise fluxes, positive into voxel.
% source: Array of voxel-wise source terms, positive into voxel.
% tracerBnd: Tracer from boundaries.
% tracerSource: Tracer from from sources.
% time:
% tof: time of flight (for full voxel).
% dtof: filling time for voxel.
% srcSign: Should be 1.0 for standard config, and -1.0 for cap-bed calc.
% porVolume:
%
%%%%Output:
%
% dwnTracer: Dynamic tracer distribution - cell outflow.
% volTracer: Dynamic tracer distribution - cell average.
%
%                north:4
%               ___________
%              |           |
%              |  top:6    |
%       west:1 |           | east:2
%              | bottom:5  |
%              |           |
%              |___________|
%                south=3
%
%
function [upTracer, dwnTracer, volTracer] = tracerDynForward3D(nx, ny, nz, tos,flux,source,tracerBnd,tracerSource,time,tof,dtof,srcSign,porVolume,mask)
  mm = length(time);
  nn = nx*ny*nz;
  upTracer = zeros(nn,mm);
  dwnTracer = zeros(nn,mm);
  volTracer = zeros(nn,mm);
  dt = time(2:mm) - time(1:mm-1);
  assert(max(dt) < min(dt) + 1e-6);
  Dt = time(2)-time(1);
  
  tol = 1.0e-15;
  indexNatural = 1:nn;
  [I,J,K] = ind2sub([nx ny nz],indexNatural);
  idUp = zeros(nn,6);
  idUp(:,1) = (indexNatural' - 1) .* (I(:) > 1) - (I(:) == 1).*(ny*(K(:)-1)+J(:));
  idUp(:,2) = (indexNatural' + 1) .* (I(:) < nx) - (I(:) == nx).*(ny*(K(:)-1)+J(:));
  idUp(:,3) = (indexNatural' - nx) .* (J(:) > 1) - (J(:) == 1).*(nx*(K(:)-1)+I(:));
  idUp(:,4) = (indexNatural' + nx) .* (J(:) < ny) - (J(:) == ny).*(nx*(K(:)-1)+I(:));
  idUp(:,5) = (indexNatural' - nx*ny) .* (K(:) > 1) - (K(:) == 1).*(nx*(J(:)-1)+I(:));
  idUp(:,6) = (indexNatural' + nx*ny) .* (K(:) < nz) - (K(:) == nz).*(nx*(J(:)-1)+I(:));
  
  for ii = 1:nn
    vxl = tos(ii);
%     if (vxl==54077)
%         vxl;
%     end
    if (mask(vxl))
        
    if (srcSign < 0)
      fluxUp = tracerSource(vxl,:)*max(0,source(vxl));
    else
      fluxUp = zeros(1,mm);
    end
    fluxDown = max(0,-srcSign*source(vxl));
    fluxUp0 = 0;
    
    for face=1:6
      if (flux(vxl,face) > tol)
        tracerUp = zeros(1,mm);
        fluxUp0 = fluxUp0 + flux(vxl,face);
        if (idUp(vxl,face)>0)
          tracerUp = dwnTracer(idUp(vxl,face),:);
        elseif (idUp(vxl,face)<0)
          tracerUp = tracerBnd(-idUp(vxl,face),:,face);
        end
        fluxUp = fluxUp + tracerUp*flux(vxl,face);
      elseif (flux(vxl,face) < -tol)
        fluxDown = fluxDown - flux(vxl,face);
      end
    end
    
    if (srcSign*source(vxl) > 0)
%       if (fluxUp0 < tol)
%         ms = mm;
%       else
%         ms = min(mm,floor(dtof(vxl)/Dt)+1);
%       end
      ms = min(mm,floor(dtof(vxl)/Dt)+1);
      if (ms<0)
          ms=-ms;
      end

      qhat = srcSign*source(vxl) / porVolume(vxl);
      expqhat = exp(-qhat*Dt);
      vDecay = [1 cumprod(expqhat*ones(1,ms-1))];
      DtRes = min(Dt,max(0, dtof(vxl)-(ms-1)*Dt));
      expqhatRes = exp(-qhat*DtRes);
      vDecay = [vDecay expqhatRes*vDecay(ms)];
      wght = DtRes/Dt;
      sampleTau = qhat*Dt*vDecay;
      if (ms > 1)
        sampleTau(1) = sampleTau(1)*0.5;
        sampleTau(ms) = sampleTau(ms)*0.5*(1+wght) + sampleTau(ms+1)*wght*(1-wght)*0.5;
        sampleTau(ms+1) = sampleTau(ms+1)*wght*wght*0.5; %%% fractional dt ...
      else
        sampleTau(ms) = sampleTau(ms)*wght*0.5 + sampleTau(ms+1)*wght*(1-wght)*0.5;
        sampleTau(ms+1) = sampleTau(ms+1)*wght*wght*0.5; %%% fractional dt ...
      end
      shiftTau = zeros(1,ms+1);
      shiftTau(ms+1) = wght*vDecay(ms+1);
      shiftTau(ms) = (1-wght)*vDecay(ms+1);
    else
%       ms = floor(dtof(vxl)/Dt)+1;
      ms = min(mm,floor(dtof(vxl)/Dt)+1);
      if (ms>mm)
        ms=mm;
        shiftTau =  zeros(1,1);
      else
        DtRes = min(Dt,max(0, dtof(vxl)-(ms-1)*Dt));
        wght = DtRes/Dt;
      
        diffusion = [.1 .9 .9 .1];
        %diffusion = [1 1];
        diffStep = length(diffusion)/2;
        shiftTau = zeros(1,ms+diffStep);        
        % forward
        for I = 1:diffStep
            shiftTau(ms+I) = diffusion(diffStep+I) * wght;
        end
        % backward
        tauInd = ms;
        for I = 0:diffStep-1
            shiftTau(tauInd) = shiftTau(tauInd) + diffusion(diffStep-I) * (1-wght);
            if tauInd > 1
                tauInd = tauInd - 1;
            end
        end
      end
      
      sampleTau = (1.0/(ms))*ones(1,ms);
    end   
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    
%%%% Ser ut til at conv er såpass effektiv for shift at det ikke er grunn
%%%% til håndkoding.  En annen ting er at en sampling av volTracer
%%%% ikke er nødvendig for beregningene, men kan begrenses til tidspunkt
%%%% hvor en ønsker output...
%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    if (srcSign*source(vxl) > 0)
      upTracer(vxl,:) = fluxUp/(fluxUp0+tol);
      cProfile = conv(upTracer(vxl,:),shiftTau);
      dwnTracer(vxl,:) = cProfile(1:mm);
      cSample = conv(tracerSource(vxl,:),sampleTau);
      dwnTracer(vxl,:) = dwnTracer(vxl,:) + cSample(1:mm);
      volTracer(vxl,:) = dwnTracer(vxl,:);

    else
      fluxDown = max(fluxDown,fluxUp0);
      upTracer(vxl,:) = fluxUp/(fluxDown+tol);
      cProfile = conv(upTracer(vxl,:),shiftTau);
      dwnTracer(vxl,:) = cProfile(1:mm);
      
      cSample = conv(upTracer(vxl,:),sampleTau);
      volTracer(vxl,:) = cSample(1:mm);
%       volTracer(vxl,:) = dwnTracer(vxl,:);
    end
    
    end
  end

end