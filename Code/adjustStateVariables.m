function adjustedVar=adjustStateVariables(var,options) 

% function adjustedVar=adjustStateVariables(var,options); 
% 
% Help functions while running Kalman filter. Adjust the state
% variables to be within bounds. 
%
% var     :  A list of states, each belonging to a certain state
%            variable - may be on matrix form with each column
%            being a list of states.
% options :  Options for the simulator - need to use the following
%            fields: 
%  - statevar     : list of statevariables, same order as in var.
%  - numGridBlocks: number of states belonging to each variable.
%
% Copyright(c) International Research Institute of Stavanger (IRIS)
% $Id: //depot/rfmatlab/main/eclipseKalman/adjustStateVariables.m#1 $
% $DateTime: 2011/10/13 14:36:15 $


adjustedVar=var;

% number of states of each type.
num=options.numGridBlocks;

indexStart=1;
indexStop=num;

for i=1:size(options.statevar,1)
  type=deblank(options.statevar(i,:));
  % some are bounded in interval [0,1]:
  if (strcmp(type,'SWAT') || strcmp(type,'SGAS') ...
	|| strcmp(type,'SOMAX'))
    adjustedVar(indexStart:indexStop,:)= ...
	adjustVariableWithInBounds(var(indexStart:indexStop,:),0,1);
  % some are bounded below by zero:
  elseif (strcmp(type,'RS') || strcmp(type,'RV') )
      adjustedVar(indexStart:indexStop,:)= ...
	adjustVariableWithInBounds(var(indexStart:indexStop,:),0,[]); ...
  % check if there are any tracer concentrations coming from Eclipse, 
  % if so they are bounded by zero
  elseif isfield(options,'eclTracerVar')
      for j=1:size(options.eclTracerVar,1)
        if   (strcmp(type,num2str(options.eclTracerVar(j,:))))
            adjustedVar(indexStart:indexStop,:)= ...
	adjustVariableWithInBounds(var(indexStart:indexStop,:),0,[]); ...
        end
      end
  elseif strcmp(type,'PRESSURE') 
    adjustedVar(indexStart:indexStop,:)= ...
	adjustVariableWithInBounds(var(indexStart:indexStop,:),1,[]); ...
  else
    % nothing to adjust
  end
  % move forward the counter of the state variables:
  indexStart=indexStart+num;
  indexStop=indexStop+num;
end
% if we have tracer consentrations coming from Artsim (they are always
% attached below the Eclipse state variables in var)
if isfield(options,'artsimVar') 
     for i=1:size(options.artsimVar,1)
        adjustedVar(indexStart:indexStop,:)= ...
	adjustVariableWithInBounds(var(indexStart:indexStop,:),0,[]); ...
      % move forward the counter of the state variables:
      indexStart=indexStart+num;
      indexStop=indexStop+num;
    end
end

