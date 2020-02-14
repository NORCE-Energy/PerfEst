function [ensemble,n]=respectBounds(ensemble,kalmanOptions,options)

% function ensemble=respectBounds(ensemble,kalmanOptions,options)
%
% If any variables (static or dynamic) of the ensemble is outside
% its bounds, put them insid.
%
% Copyright(c) International Research Institute of Stavanger (IRIS)
% $Id: //depot/rfmatlab/main/eclipseKalman/respectBounds.m#5 $
% $DateTime: 2015/07/30 10:23:07 $

% cy added this. see after exit the content of option should be gone.
if nargin<2
   load inputData;
end

numStaticVar=1:size(options.staticVar,1)*options.numGridBlocks;
if isfield(options,'freeparam')
    numStaticVar=numStaticVar+sum(options.freeparam);
end
nummeas=options.nummeas;
n = 0;
if isfield(kalmanOptions,'freeparamLB')
    [ensemble(1:sum(options.freeparam),:),tmp_n]= ...
      adjustVariableWithInBounds(ensemble(1:sum(options.freeparam),:), ...
				 kalmanOptions.freeparamLB, ...
				 kalmanOptions.freeparamUB);
    n = n+tmp_n;
end
% adjust staticVar within bounds if required:
if isfield(kalmanOptions,'staticVarUB')
  [ensemble(numStaticVar,:),tmp_n]= ...
      adjustVariableWithInBounds(ensemble(numStaticVar,:), ...
				 kalmanOptions.staticVarLB, ...
				 kalmanOptions.staticVarUB);
  n = n+tmp_n;           
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
  ensemble(firstDyn+1:end-nummeas,:)= ...
      adjustStateVariables(ensemble(firstDyn+1:end- ...
				    nummeas,:),options);  
end
