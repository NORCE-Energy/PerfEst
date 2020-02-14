%Computes steady state tracer distribution from abitrary set of source(s).
%
%%%%Input:
%
% nx: Number of grid cells in x-direction
% ny: Number of grid cslls in y-direction
% nz: Number of grid cslls in z-direction
% tos: Voxel indices sorted according to upstream-before-downstream.
% flux: Array of voxel-wise, face-wise fluxes, positive into voxel.
% source: Array of voxel-wise source terms, positive into voxel.
% tracerBnd: Tracer from boundaries
% tracerSource: Tracer from from sources
%
%%%%Output:
%
% tofTracer: Steady state tracer distribution.
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
function [tofTracer] = tracerSSForward3D(nx, ny, nz, tos,flux,source,tracerBnd,tracerSource)
  nn = nx*ny*nz;
  tofTracer = zeros(nn,1);
  
  tol = 1.0e-15;
  indexNatural = 1:nn;
  [I,J,K] = ind2sub([nx ny nz],indexNatural);
  idUp(:,1) = (indexNatural' - 1) .* (I(:) > 1) - (I(:) == 1).*(ny*(K(:)-1)+J(:));
  idUp(:,2) = (indexNatural' + 1) .* (I(:) < nx) - (I(:) == nx).*(ny*(K(:)-1)+J(:));
  idUp(:,3) = (indexNatural' - nx) .* (J(:) > 1) - (J(:) == 1).*(nx*(K(:)-1)+I(:));
  idUp(:,4) = (indexNatural' + nx) .* (J(:) < ny) - (J(:) == ny).*(nx*(K(:)-1)+I(:));
  idUp(:,5) = (indexNatural' - nx*ny) .* (K(:) > 1) - (K(:) == 1).*(nx*(J(:)-1)+I(:));
  idUp(:,6) = (indexNatural' + nx*ny) .* (K(:) < nz) - (K(:) == nz).*(nx*(J(:)-1)+I(:));
  
  for ii = 1:nn
    vxl = tos(ii);
    fluxUp = tracerSource(vxl)*max(0,source(vxl));
    fluxDown = max(0,-source(vxl));
    for face=1:6
      if (flux(vxl,face) > tol)
        tracerUp = 0.0;
        if (idUp(vxl,face)>0)
          tracerUp = tofTracer(idUp(vxl,face));
        elseif (idUp(vxl,face)<0)
          tracerUp = tracerBnd(-idUp(vxl,face),face);
        end
        fluxUp = fluxUp + tracerUp*flux(vxl,face);
      elseif (flux(vxl,face) < -tol)
        fluxDown = fluxDown - flux(vxl,face);
      end
    end
    tofTracer(vxl) = fluxUp/(fluxDown+tol);
  end

end
  
