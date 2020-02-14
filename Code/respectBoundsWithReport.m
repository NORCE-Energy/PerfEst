function [ensemble,n]=respectBoundsWithReport(ensemble,kalmanOptions,options)

% function ensemble=respectBounds(ensemble,kalmanOptions,options)
%
% If any variables (static or dynamic) of the ensemble is outside
% its bounds, put them insid.
%
% Copyright(c) International Research Institute of Stavanger (IRIS)
% $Id: //depot/rfmatlab/main/eclipseKalman/respectBounds.m#1 $
% $DateTime: 2011/10/13 14:36:15 $

% cy added this. see after exit the content of option should be gone.
if nargin<2
   load inputData;
end

numStaticVar=1:size(options.staticVar,1)*options.numGridBlocks;
if isfield(options,'freeparam')
    numStaticVar=numStaticVar+sum(options.freeparam);
end
nummeas=options.nummeas;

if isfield(kalmanOptions,'freeparamLB')
    [ensemble(1:sum(options.freeparam),:),n]= ...
      truncate(ensemble(1:sum(options.freeparam),:), ...
				 kalmanOptions.freeparamLB, ...
				 kalmanOptions.freeparamUB);
end
% adjust staticVar within bounds if required:
if isfield(kalmanOptions,'staticVarUB')
  [ensemble(numStaticVar,:),n]= ...
      truncate(ensemble(numStaticVar,:), ...
				 kalmanOptions.staticVarLB, ...
				 kalmanOptions.staticVarUB);
end
% avoid negative saturations and saturations greater than 1:
if ~isempty(numStaticVar)
  firstDyn=numStaticVar(end);
else
  if isfield(options,'freeparam')
    firstDyn=sum(options.freeparam);
  else
    firstDyn=0;
  end
end

if ( ~isfield(kalmanOptions,'iterES') || (isfield(kalmanOptions,'iterES') && kalmanOptions.iterES~=1) ) && ~isempty(options.statevar)
  [ensemble(firstDyn+1:end-nummeas,:),n]= ...
      truncate(ensemble(firstDyn+1:end-nummeas,:),options);  
end


function [variable,n]=truncate(variable,lowerbound,upperbound)

% function variable=adjustVariableWithInBounds(variable, ...
%						lowerbound, ...
%						upperbound)
%
% The function returns variable such that 
%    lowerbound <= variable(i) <= upperbound,
% i.e. if variable(i)<lowerbound before calling the function, then
% variable(i)=lowerbound as output. Similarly for upperbound, both
% with inequality the opposite way. If there is no
% lowerbound/upperbound, set the corresponding matrix empty.
%
% Copyright(c) International Research Institute of Stavanger (IRIS)
% $Id: //depot/rfmatlab/main/eclipseKalman/adjustVariableWithInBounds.m#1 $
% $DateTime: 2011/10/13 14:36:15 $


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
