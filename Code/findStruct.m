function [outStruct,index] = findStruct(s,field,value)


% function [outStruct,index] = findStruct(s,field,value)
%
% Find structure with given field and value.
%
% Input
% -----
%
% s                : Array of input structures.
% field            : Search in this field.
% value            : Search for this value.
%
% Output
% ------
%
% outStruct        : Output structure.
% index            : Index in array corresponding to ouput
%                    strucutre. 
%
% Copyright(c) International Research Institute of Stavanger (IRIS)
% $Id: //depot/rfmatlab/main/eclipseKalman/findStruct.m#1 $
% $DateTime: 2011/10/13 14:36:15 $


outStruct = struct([]);
index = 0;
for i = 1:length(s)
  if strcmp(field,'date')
    if isequal(datenum(getfield(s(i),field)),...
	       datenum(value)) %#ok<*GFLD>
      index = i;
      outStruct = s(i);
      return
    end
  else
    if isequal(getfield(s(i),field),value)
      index = i;
      outStruct = s(i);
      return
    end
  end
end