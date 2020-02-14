%Computes an ordered sequence of cells according to flow dependence
%starting from inflow boundary, such that an upstream cell is ordered 
%before its downstream dependants.
%
%%%%Input:
%
% nx: Number of grid cells in x-direction
% ny: Number of grid cslls in y-direction
% nz: Number of grid cslls in z-direction
% velX: matrix of (nx+1)*ny*nz cell-edge velocities
% velY: matrix of nx*(ny+1)*nz cell-edge velocities
% velZ: matrix of nx*ny*(nz+1) cell-edge velocities
%
%%%%Output:
%
% tos: Permutation of 1:nx*ny*nz representing a map from natural voxel
%      indexing to an upstream-before-downstream sequence.
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
function [tos] = upstreamSort3D(nx, ny, nz, velX, velY, velZ)
  nn = nx*ny*nz;
  tos = zeros(nn,1);
 
  tol = 1.0e-20;
  indexNatural = 1:nn;
  isUp = zeros(nx,ny,nz,6);
  isUp(2:nx,:,:,1) = (velX(2:nx,:,:) > tol);
  isUp(1:nx-1,:,:,2) = (velX(2:nx,:,:) < -tol);
  isUp(:,2:ny,:,3) = (velY(:,2:ny,:) > tol);
  isUp(:,1:ny-1,:,4) = (velY(:,2:ny,:) < -tol);
  isUp(:,:,2:nz,5) = (velZ(:,:,2:nz) > tol);
  isUp(:,:,1:nz-1,6) = (velZ(:,:,2:nz) < -tol);
  
  isUp = reshape(isUp,[nn,6]);
  isUp = repmat(indexNatural',1,6) .* isUp;
  isUp(:,1) = isUp(:,1) - (isUp(:,1) > 0);
  isUp(:,2) = isUp(:,2) + (isUp(:,2) > 0);
  isUp(:,3) = isUp(:,3) - nx*(isUp(:,3) > 0);
  isUp(:,4) = isUp(:,4) + nx*(isUp(:,4) > 0);
  isUp(:,5) = isUp(:,5) - nx*ny*(isUp(:,5) > 0);
  isUp(:,6) = isUp(:,6) + nx*ny*(isUp(:,6) > 0);

% Wikipedia:
% 1  procedure DFS-iterative(G,v):
% 2      let S be a stack
% 3      S.push(v)
% 4      while S is not empty
% 5          v = S.pop()
% 6          if v is not labeled as discovered:
% 7              label v as discovered
% 8              for all edges from v to w in G.adjacentEdges(v) do
% 9                  S.push(w)

  ntos = 0;
  vxlDone = zeros(nn,1);
  vxlStack = zeros(nn,1);
  for vxlStart = 1:nn
    if (vxlDone(vxlStart) == 0)
      vxlPt = 1;
      vxlStack(vxlPt) = vxlStart;
      while (vxlPt > 0)
        vxl = vxlStack(vxlPt);
        if (vxlDone(vxl) == 0)
          vxlPt0 = vxlPt;
          for face = 1:6
            vxlNb = isUp(vxl,face);
            if (vxlNb>0 && vxlDone(vxlNb) == 0)
              vxlPt = vxlPt + 1;
              vxlStack(vxlPt) = vxlNb;
            end
          end
          if (vxlPt0 == vxlPt)
            ntos = ntos + 1;
            tos(ntos) = vxl;
            vxlDone(vxl) = 1;
            vxlStack(vxlPt) = 0;
            vxlPt = vxlPt - 1;
          end
        else
          vxlStack(vxlPt) = 0;
          vxlPt = vxlPt - 1;
        end
      end
    end    
  end
end
  
