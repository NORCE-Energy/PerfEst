function [x] = fastGaussian3d( dimension, Sdev, Corr )
% function [x] = fastGaussian3d( dimension, Sdev, Corr )
% 
% Generates random vector from distribution satisfying Gaussian
% variogram in dimension up to 3-d.
%
% Input:
% - dimension  Dimension of grid 
% - Sdev       Standard deviation
% - Corr       Correlation length, in units of block length.
%              Corr may be replaced with a vector of length 3 with correlation
%              length in x-, y- and z-direction.
%
% Output:
% - x          Random vector.
%
%
% The parameterization of the grid is assumed to have size
% dimension, if dimension is a vector, or [dimension,dimension] if
% dimension is scalar. 
% The coefficients of the grid is assumed to be reordered
% columnwise into the parameter vector.
% The grid is assumed to have a local basis.
%
% Example of use:
%
% Want to generate a field on a 2-d grid with dimension m x n, with
% correlation length a along first coordinate axis, b along second
% coordinate axis and standard deviation sigma:
%
% x=fastGaussian([m n p],sigma,[a b c]);
%
% If the dimension is nxn one can write
%
%  x=fastGaussian(n,sigma,[a b c]);
%
% If the correlation length is the same in both directions:
%
%   x=fastGaussian([m n p],sigma,a);  or
%   x=fastGaussian(n,sigma,a); 
%
%
% The theory behind this algorithm is not published yet.
% The properties on the Kronecker product behind this algorithm can
% be found in Horn & Johnson: Topics in Matrix Analysis, Cambridge
% UP, 1991. 
%
% Note that we add a small number on the diagonal of the covariance
% matrix to avoid numerical problems with Cholesky decomposition (a
% nugget effect).
%
% Copyright(c) International Research Institute of Stavanger (IRIS)
% $Id: //depot/rfmatlab/main/stats/fastGaussian3d.m#1 $
% $DateTime: 2011/10/13 14:36:15 $


% Initialize dimension:
dimension=dimension(:); % to trap matrices in input
m = dimension(1) ;
if length(dimension) < 3
  % use the 2-d version if the dimension has length 1 or 2:
  warning('rfmatlab:stats','fastGaussian3d is called with dimension less than 3');
  [x] = fastGaussian( dimension, Sdev, Corr );
  return
elseif length(dimension)== 3
  n = dimension(2) ;
  p = dimension(3) ;
else
  error(['fastGaussian: Wrong input, dimension should have length', ...
	 ' at most 3'])
end
mxnxp = m*n*p ;

% Compute variance.
if max(size(Sdev))>1 % check input
  variance=1;
else
  variance = Sdev ; % the variance will come out through the
                  % kronecker product.
end


% Initialize correlation length:
Corr=Corr(:); % to trap matrices in input
if length(Corr)==2 || length(Corr)>3
  error(['fastGaussian3d: Wrong input, Corr should have length', ...
	 ' 1 or 3'])
end
if length(Corr)==1
  Corr(2)=Corr(1);
  Corr(3)=Corr(1);
end

if nargin<3
  error('fastGaussian3d need three input variables')
end

% 
% first generate the covariance matrix for one layer:
dist1=0:m-1;
dist1=dist1/Corr(1);
T1=toeplitz(dist1);
% to avoid problem with Cholesky factorization when the matrix is
% close to singular we add a small number on the diagonal entries. 
T1=variance*exp(-T1.^2)+1e-10*eye(m);
% Cholesky decomposition for one layer:
cholT1=chol(T1);

% generate the covariance matrix for the second layer:
% to save time - use a copy if possible:
if (Corr(1)==Corr(2)) && n==m
  cholT2=cholT1;
else 
  % same as for the first dimension:
  dist2=0:n-1;
  dist2=dist2/Corr(2);
  T2=toeplitz(dist2);
  T2=variance*exp(-T2.^2)+1e-10*eye(n);
  cholT2=chol(T2);
end

% generate the covariance matrix for the third layer:
% use variance = 1 to get the correct value
dist3=0:p-1;
dist3=dist3/Corr(3);
T3=toeplitz(dist3);
T3=exp(-T3.^2)+1e-10*eye(p);
cholT3=chol(T3);


% draw a random variable:
x=randn(mxnxp,1);

% adjust to get the correct covariance matrix,
% applying Lemma 4.3.1. in Horn & Johnson:
% we need to adjust to get the correct covariance matrix - either
% dimension 1 and 2 or 2 and 3 need to be grouped together.
if n<=p
  x=kron(cholT2',cholT1')*reshape(x,m*n,p)*cholT3;
else
  x=cholT1'*reshape(x,m,n*p)*kron(cholT3,cholT2);
end
  
% reshape back
x=x(:);




if max(size(Sdev))>1
  if min(size(Sdev))==1 && length(Sdev)==length(x)
    x=Sdev.*x;
  else
    error('fastGaussian3d: Inconsistent dimension of Sdev')
  end
end
