function Vname=convertFieldName(Vname)
  
% There are certain tags in restart files of eclipse that casuse problems
% while used as fields in a matlab structure. This function is used to
% convert these names.
%
% % Copyright(c) International Research Institute of Stavanger (IRIS)
% $Id: //depot/rfmatlab/main/eclipse/convertFieldName.m#1 $
% $DateTime: 2011/10/13 14:36:15 $

  
% The field name can not contain spaces, and not end with a plus or
% minus. The plus and minus are replace with a '_P' and '_M' respectively.
if strcmp(Vname(end),'+')
  Vname(end:end+1)='_p';
end
if strcmp(Vname(end),'-')
  Vname(end:end+1)='_m';
end
% field names are not allowed to start with a number:
if regexp(Vname(1),'\d')
  Vname=['num_',Vname];
end
% field names can not contain a "/":
if ~isempty(strfind(Vname,'/'))
  ind=strfind(Vname,'/');
  Vname=[Vname(1:ind-1),'_s_',Vname(ind+1:end)];
end
