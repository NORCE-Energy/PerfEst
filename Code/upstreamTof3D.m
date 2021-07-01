%Computes time-of-flight from upstream boundary based on a voxel ordering
%of upstream-before-downstream defined by input tos.
%
%%%%Input:
%
% nx: Number of grid cells in x-direction
% ny: Number of grid cslls in y-direction
% nz: Number of grid cslls in z-direction
% tos: Voxel indices sorted according to upstream-before-downstream.
% flux: Array of voxel-wise, face-wise fluxes, positive into voxel.
% source: Array of voxel-wise source terms, positive into voxel.
% porvol: Array of voxel-wise pore volumes.
% extSource: Correction for external source terms (ext vessel terminals)
%
%%%%Output:
%
% tof: Array of time-of-flight from inflow boundary for each voxel.
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
function [tof, dtof] = upstreamTof3D(nx, ny, nz, tos,flux,source,porvol,extSource)
  nn = nx*ny*nz;
  tof = zeros(nn,1);
  dtof = zeros(nn,1);
  
  tol = 1.0e-20;
  indexNatural = 1:nn;
  [I J K] = ind2sub([nx ny nz],indexNatural);
  idUp = zeros(nn,6);
  idUp(:,1) = (indexNatural' - 1) .* (I(:) > 1);
  idUp(:,2) = (indexNatural' + 1) .* (I(:) < nx);
  idUp(:,3) = (indexNatural' - nx) .* (J(:) > 1);
  idUp(:,4) = (indexNatural' + nx) .* (J(:) < ny);
  idUp(:,5) = (indexNatural' - nx*ny) .* (K(:) > 1);
  idUp(:,6) = (indexNatural' + nx*ny) .* (K(:) < nz);

  for ii = 1:nn
    vxl = tos(ii);
    fluxUp = 0.0;
    fluxUpTof = max(0.0,source(vxl));
    fluxDown = 0.0;
    for face=1:6
      if (flux(vxl,face) > tol)
        if (idUp(vxl,face) > 0)
          fluxUp = fluxUp + flux(vxl,face);
          fluxUpTof = fluxUpTof + tof(idUp(vxl,face))*flux(vxl,face);
        elseif (idUp(vxl,face) == 0)
          fluxUp = fluxUp + flux(vxl,face);
        end
      elseif (flux(vxl,face) < -tol)
        fluxDown = fluxDown - flux(vxl,face);
      end
    end
    if (extSource(vxl) > tol)
      fluxUp = fluxUp + extSource(vxl);
      %Assuming that incoming external source carries zero tof
      %fluxUpTof = fluxUpTof + 0.0*extSource;
    elseif (extSource(vxl) < tol)
      fluxDown = fluxDown - extSource(vxl);
    end
    tofUp = fluxUpTof/(fluxUp+max(0.0,fluxDown-fluxUp)+tol);
    %dtof(vxl) = porvol(vxl)/(0.5*(fluxDown+fluxUp)+tol);
    dtof(vxl) = porvol(vxl)/(fluxUp+max(0.0,fluxDown-fluxUp)+tol);
% % %     if (fluxDown-fluxUp > tol)
% % %       dtof(vxl) = -porvol(vxl)/(fluxDown-fluxUp)*log((fluxUp+tol)/(fluxDown+tol));
% % %     end
    tof(vxl) = tofUp + dtof(vxl);
  end

end