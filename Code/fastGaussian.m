function [x] = fastGaussian( dimension, Sdev, Corr )
% function [x] = fastGaussian( dimension, Sdev, Corr )
% 
% Generates random vector from distribution satisfying Gaussian
% variogram in 2-d. 
%
% Input:
% - dimension  Dimension of grid 
% - Sdev       Standard deviation
% - Corr       Correlation length, in units of block length.
%              Corr may be replaced with a vector of length 2 with correlation
%              length in x- and y- direction.
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
% x=fastGaussian([m n],sigma,[a b]);
%
% If the dimension is nxn one can write
%
%  x=fastGaussian(n,sigma,[a b]);
%
% If the correlation length is the same in both directions:
%
%   x=fastGaussian([m n],sigma,a);  or
%   x=fastGaussian(n,sigma,a); 


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
% $Id: //depot/rfmatlab/main/stats/fastGaussian.m#1 $
% $DateTime: 2011/10/13 14:36:15 $


% Initialize dimension:
dimension=dimension(:); % to trap matrices in input
m = dimension(1) ;
if length(dimension) == 1
  n = m ;
elseif length(dimension)== 2
  n = dimension(2) ;
else
  error(['fastGaussian: Wrong input, dimension should have length', ...
	 ' at most 2'])
end
mxn = m*n ;

% Compute variance.
if max(size(Sdev))>1 % check input
  % error('fastGaussian: Wrong input, Sdev should be a scalar')
  variance=1;
else
  variance = Sdev ; % the variance will come out through the
                  % kronecker product.
end


% Initialize correlation length:
Corr=Corr(:); % to trap matrices in input
if length(Corr)==1
  Corr(2)=Corr(1);
end
if length(Corr)>2
  error(['fastGaussian: Wrong input, Corr should have length', ...
	 ' at most 2'])
end

if nargin<3
  error('fastGaussian need three input variables')
end

% 
% first generate the covariance matrix for one layer:
dist=0:m-1;
dist=dist/Corr(1);
T=toeplitz(dist);
% to avoid problem with Cholesky factorization when the matrix is
% close to singular we add a small number on the diagonal entries. 
T=variance*exp(-T.^2)+1e-10*eye(m);
% Cholesky decomposition for one layer:
cholT=chol(T);

% generate the covariance matrix for the second layer:
% to save time - use a copy if possible:
if (Corr(1)==Corr(2)) && n==m
  cholT2=cholT;
else 
  % same as for the first dimension:
  dist2=0:n-1;
  dist2=dist2/Corr(2);
  T2=toeplitz(dist2);
  T2=variance*exp(-T2.^2)+1e-10*eye(n);
  cholT2=chol(T2);
end
  
% draw a random variable:
x=randn(mxn,1);

% adjust to get the correct covariance matrix,
% applying Lemma 4.3.1. in Horn & Johnson:
x=cholT'*reshape(x,m,n)*cholT2;

% reshape back
x=x(:);



if max(size(Sdev))>1
  if min(size(Sdev))==1 && length(Sdev)==length(x)
    x=Sdev.*x;
  else
    error('fastGaussian: Inconsistent dimension of Sdev')
  end
end
