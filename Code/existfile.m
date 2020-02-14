function bool=existfile(filename)

% function bool=existfile(filename)
%
% If a file with named filename exist, bool is set 1, otherwise 0.
% The aim of this function is to get a function that is compilable,
% and returns similar information as the built-in function
% exist(filename,'file') (which is not supported by the C++
% library).
%
% The function search through the hole path defined in matlab,
% which may give something that looks like different behaviour when
% the function is compiled. 
% 
% To check if the file 'test.mat' exist in current directory use
% 
% bool=existfile('./test.mat')
%
% or in any directory in the matlab path
%
% bool=existfile('test.mat')
%
%
% Copyright (c) 2001-2002 RF-Rogaland Research. All rights reserved.
% Version 1.1. Date 03/06/2002.
%
% Geir Naevdal 26.10.01.
% 030602 GEN: Improved documentation.

% Copyright (c) 2006 IRIS, All Rights Reserved.
% $Id: //depot/rfmatlab/main/iofun/existfile.m#2 $
% $DateTime: 2015/07/21 15:24:48 $

% assume that the file exists if it is possible to open it.
check=fopen(['./',filename]);
if check>0
    bool=1;
    fclose(check);
else
    bool=0;
end
