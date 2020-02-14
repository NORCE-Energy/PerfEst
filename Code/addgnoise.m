function [Y,RTSIGMA] = addgnoise( Ytrue,SIGMA,SQ )
%addgnoise - Add Gaussian noise with given covariance structure
%
% function  [Y,RTSIGMA] = addgnoise(Ytrue,SIGMA,SQ);
%
% Add noise, normally distributed, with covariance given by SIGMA.
%
% Input:
% - Ytrue    Original signal.
% - SIGMA    Specification of covariance matrix of noise.  May be
%            entered as scalar, vector or full matrix. If it SIGMA
%            is a vector, then it is interpreted as the covariance
%            matrix is diag(SIGMA).
% - SQ       If present, determine whether SIGMA or SIGMA*SIGMA' is used
%            as the covariance matrix.  Thus, if the square root of
%            the covariance matrix has allready been calculated
%            previously, work may be saved by setting SQ to 1.
%
% Output:
% - Y        Signal with noise added.
% - RTSIGMA  The square root of SIGMA; RTSIGMA*RTSIGMA' = SIGMA.
%            (Helpful if it is cumbersome to compute).
%
% Copyright(c) International Research Institute of Stavanger (IRIS)
% $Id: //depot/rfmatlab/main/stats/addgnoise.m#1 $
% $DateTime: 2011/10/13 14:36:15 $


% Compute the normally distributed noise, with covariance matrix
% SIGMA or SIGMA*SIGMA'.
%
try
  if nargin > 2 && SQ == 1
    % Use SIGMA*SIGMA' as covariance matrix
    RTSIGMA = SIGMA ;
    if min(size(SIGMA)) == 1
      % SIGMA is a scalar or vector
      error = RTSIGMA.*randn(size(Ytrue)) ;
    else
      error = RTSIGMA*randn(size(RTSIGMA,2),1) ;
    end
  else
    % Use SIGMA as covariance matrix
    if min(size(SIGMA))==1
      % SIGMA is entered as a scalar or a vector
      RTSIGMA = realsqrt(SIGMA);
      error = RTSIGMA.*randn(size(Ytrue));
    else
      [RTSIGMA,p] = chol(SIGMA); % The matrix must be transposed.
      if p>0
	disp('Problem with Cholesky factorization')
	disp(['p = ',num2str(p)]);
	RTSIGMA = real(sqrtm(SIGMA));
	disp('Finnaly - we got a square root!')
      end
      RTSIGMA = RTSIGMA';
      error = RTSIGMA*randn(size(Ytrue));
    end %if  
  end
  % Add the noise:
  Y = Ytrue+error;
catch err
  disp('Error in addgnoise');
  disp('Size Ytrue:')
  size(Ytrue)
  disp('Size SIGMA:')
  size(SIGMA)
  rethrow(err);
end
