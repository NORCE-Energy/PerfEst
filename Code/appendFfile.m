function appendFfile(filecontents,filename,options)

% Write an eclipse  file on the F***** format.
% Assuming that the the variables up to the state variables are
% already placed on the file 'COMMON.F'.
%
% function appendFfile(filecontents,filename,options)
%
% filename:         The name of the restart file.
% filecontents:     The contents of the restart file.
% options.APPEND:   If this field exists, is is assumed that the
%                   begining of the restartfile is already written
%                   to a file named COMMON.F, otherwise a file
%                   COMMON.F is written for later use.
% options.statevar: The filecontents from the first field listed in
%                   statevar is appendend to COMMON.F.
%
% See also readFfile, writeFfile, writetoFfile and forweclipse.
% 
%  Copyright(c) International Research Institute of Stavanger (IRIS)
%  $Id: //depot/rfmatlab/main/eclipse/appendFfile.m#1 $
%  $DateTime: 2011/10/13 14:36:15 $
  
  
% find first state variable:
j=0;
ok=0;
if ~isempty(filecontents)
  while ok==0
    j=j+1;
    for i=1:size(options.statevar,1)
      v(i)=strcmp(filecontents.varnames(j,:),options.statevar(i,:)); %#ok<AGROW>
    end
    ok=max(v);
  end
  firststatevar=j;
else
  firststatevar=1;
end


% If the file COMMON.F exists it is assumed that this file contains
% the proper contents, otherwise it is created.
if ~isfield(options,'APPEND')
  firstvar=1;
  fid=fopen('COMMON.F','w');
  if fid<=0
    error('Could not open file COMMON.F')
  end
  writetoFfile(fid,firstvar,firststatevar-1,filecontents);
  fclose(fid);
end
sysCom=['cp COMMON.F ',filename];
unix(sysCom);
fid=fopen(filename,'a');
if fid<=0
  error(['Could not open file ',filename])
end
writetoFfile(fid,firststatevar,size(filecontents.varnames,1),filecontents);
fclose(fid);
