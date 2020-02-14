function [variable,n]=adjustVariableWithInBounds(variable, ...
						lowerbound, ...
						upperbound)

% function [variable,n]=adjustVariableWithInBounds(variable, ...
%						lowerbound, ...
%						upperbound)
%
% The function returns variable such that 
%    lowerbound <= variable(i) <= upperbound,
% i.e. if variable(i)<lowerbound before calling the function, then
% variable(i)=lowerbound as output. Similarly for upperbound, both
% with inequality the opposite way. If there is no
% lowerbound/upperbound, set the corresponding matrix empty.
% -------
% INPUTS:
% -------
% variable:     Vector or matrix. Variables (or an ensemble of variable 
%               samples) to be checked and truncated (if necessary) 
% lowerbound:   Scalar or vector. Lower bound(s) for the variables to be
%               checked
% upperbound:   Scalar or vector. Upper bound(s) for the variables to be
%               checked
% --------
% OUTPUTS:
% --------
% variable:     Variables after check/truncation
% n:            Average number of truncations
% Copyright(c) International Research Institute of Stavanger (IRIS)
% $Id: //depot/rfmatlab/main/eclipseKalman/adjustVariableWithInBounds.m#2 $
% $DateTime: 2015/07/30 11:20:56 $


if nargin < 3
  error(['adjustVariableWithInBounds is only implemented for three' ...
	 ' arguments'])
end

ne=size(variable,2);
n=0;
if ~isempty(lowerbound)
    if max(size(lowerbound))==1
        n=sum(variable<lowerbound);
        variable(variable<lowerbound)=lowerbound;
    else
        lowerbound=kron(lowerbound,ones(1,size(variable,2)));
        n=sum(variable<lowerbound);
        variable(variable<lowerbound)=lowerbound(variable<lowerbound);
    end
end

if ~isempty(upperbound)
    if max(size(upperbound))==1
        n=n+sum(variable>upperbound);
        variable(variable>upperbound)=upperbound;
    else
        upperbound=kron(upperbound,ones(1,ne));
        n=n+sum(variable>upperbound);
        variable(variable>upperbound)=upperbound(variable>upperbound);
    end
end

n=mean(n);