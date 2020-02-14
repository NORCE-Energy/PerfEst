function index=variableIndex(options)

% index=variableIndex(options)
%
% varIndex=measIndexList - reads the variable options from inputData.
% varIndex=measIndexList(options) - use submitted variable options.
%
% Returns varIndex a structued variable containing information about the
% variables in the state vector.
%
% index.name:   Name of variable in state vector.
% index.first:  First index of corresponding variable.
% index.last:   Last index of corresponding variable.
%
% Copyright(c) International Research Institute of Stavanger (IRIS)
% $Id: //depot/rfmatlab/main/eclipseKalman/variableIndex.m#2 $
% $DateTime: 2012/12/19 10:17:32 $


if nargin <1
    load('inputData','options')
end

% Variables are stored in the following manner:
%
% 1) Free parameters
% 2) Static fields
% 3) Dynamical fields
%
% In some case not all of the above are present.
%
% For each we need to find the dimension.

name=[]; % name of variable
first=[];% first index of variable
last=0; % last index of variable, we start with a zero to be removed at
% the end!

if isfield(options,'freeparam')
    if isempty(name)
        name=char('freeparam');
    else
        name=char(name,'freeparam');
    end
    first(end+1)=last(end)+1;
    last(end+1)=last(end)+sum(options.freeparam);
end
if isfield(options,'staticVar')
    for I=1:size(options.staticVar,1)
        if isempty(name)
            name=char(options.staticVar(I,:));
        else
            name=char(name,options.staticVar(I,:));
        end
        first(end+1)=last(end)+1;
        last(end+1)=last(end)+options.numGridBlocks;
    end
end
if isfield(options,'statevar')
    for I=1:size(options.statevar,1)
        if isempty(name)
            name=char(options.options.statevar(I,:));
        else
            name=char(name,options.statevar(I,:));
        end
        first(end+1)=last(end)+1;
        last(end+1)=last(end)+options.numGridBlocks;
    end
end

last(1)=[];

% property name
index.('name')=name;
index.('first')=first;
index.('last')=last;