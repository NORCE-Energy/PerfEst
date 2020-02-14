function runstatus = start_ews(commandline,basefilename,erroraction,ecl_path)
%
% a wrapper for calling eclipse with error checking
% before the call.
%
% runstatus = start_ews(commandline,basefilename,erroraction,ecl_path)
%
% Input:
%  CommandLine  - The script to use for the simulator call.  This
%                 is either the default '@eclipse' macro, or a
%                 user-supplied script.
%  BaseFilename - The base-name of the file on which eclipse should
%                 be run.  E.g., 'BASE2D'.
%  ErrorAction  - What to do if an error occured. The input
%                 variable is a structure with one field per
%                 possible error source. Each field is assigned a
%                 value of S/R (Stop/Retry) to specify the
%                 action taken if this error is encountered.
%                 The Stop action sends an error message to the
%                    user and stops the execution of the program.
%                 The Retry action sends a warning message to the
%                    user and tris to run the simulator again after
%                    a short pause. The duration of the pause can
%                    be specified in ErrorAction.pauselength. The
%                    default pause is 1 min. In addition, a maximum
%                    retry count may be specified in
%                    ErrorAction.maxretry.
% 1   ErrorAction.nolicserv - What to do when the contact with the
%                             license server is broken.
% 2   ErrorAction.nolicense - What to do when there is no free license
%  ecl_path       - The path to the eclipse installation directory
%
% Output:
%  RunStatus    - The function returns 0 if no error was detected.
%                 The function returns -1 if an error was detected.
%

% Copyright(c) RF-Rogaland Research 2004.
% Version 1.0. Kristian Thulin.

% Copyright (c) 2006 IRIS, All Rights Reserved.
% $Id: //depot/rfmatlab/main/eclipse/start_ews.m#1 $
% $DateTime: 2011/10/13 14:36:15 $

%====================
% Input check
%====================
if nargin < 2
  error('Need name of Eclipse data file as second argument.')
end
if nargin < 3
  erroraction=[];
end

% commandline
%-------------
% This should be an executable file. We rely on the user to specify
% this correctly. An error here will probably give a category 4 error.
% basefilename
%--------------
% This should point to a file in the current directory.
a=existfile([basefilename,'.DATA']);
if ~a
  error(['Eclipse input file ',basefilename,' not found in current directory'])
end

% erroraction
%-------------
% If nonempty, this should be a structure.
if isstruct(erroraction)
  if isfield(erroraction,'nolicserv')
    onerr1 = erroraction.nolicserv;
  else
    onerr1 = 'R';
  end
  if isfield(erroraction,'nolicense')
    onerr2 = erroraction.nolicense;
  else
    onerr2 = 'R';
  end
  if isfield(erroraction,'pauselength')
    user_delay = erroraction.pauselength;
  else
    user_delay = 1;
  end
  if isfield(erroraction,'maxretry')
    maxretry = erroraction.maxretry;
  else
    maxretry = 5;
  end
else
  onerr1 = 'R';
  onerr2 = 'R';
  user_delay = 1;
  maxretry = 5;
end

runstatus = 0;
errflag = 0;

%================================================
% Set up strings for the program calls
%================================================
% The program path is different on
% tor and loke (and on the Windows computers).
% In RF we have:
%  machine  return string from 'computer'
%   tor         SOL2
%   loke        GLNX86
%   frigg/odin  HP700
%
% frigg and odin will not run eclipse

if isunix
    c = computer;
    if strcmp(c,'SOL2')
        lm_path = '/opt/ecl/tools/flexlm72/lmutil';
        lm_string = [lm_path,' lmstat -f eclipse > ./lm.txt 2>&1'];
        %cfg_file = '/opt/ecl/macros/CONFIG.ECL';
        %ecl_path = '/opt/ecl/2001a_1/eclipse/source/ultra/eclipse.exe';
    elseif strcmp(c,'GLNXA64')
        if nargin < 4
            ecl_path = '/appl/ecl';
        end
        lm_path = [ecl_path,'/tools/linux_x86_64/flexlm113/lmutil'];
        lm_string = [lm_path,' lmstat -c ',ecl_path,'/macros/license.dat -f eclipse > ./lm.txt 2>&1'];
    elseif strcmp(c,'GLNX86')
        lm_path = '/appl/ecl/tools/linux/flexlm72/lmutil';
        lm_string = [lm_path,' lmstat -c /appl/ecl/macros/license.dat -f eclipse > ./lm.txt 2>&1'];
        %cfg_file = 'appl/ecl/macros/CONFIG.ECL';
        %ecl_path = '/appl/ecl/tools/linux/flexlm72/eclipse.exe';
    else
        %cfg_file = '';
        %ecl_path = '';
        error('Eclipse will not run on this machine')
    end
else
    % Windows machine
    lm_path = 'c:\ecl\tools\flexlm72\lmutil.exe';
    %cfg_file = 'c:\ecl\macros\CONFIG.ECL';
    %ecl_path = 'c:\ecl\2001a_1\eclipse\source\eclipse.exe';
    lm_string = [lm_path,' lmstat -f eclipse > lm.txt'];
end

if isempty(commandline)
  %ecl_string = ['cp ',cfg_file,' ./CFG.ECL;',...
  %     ecl_path,' ',basefilename,...
  %     ' Run 1728000 1728000 eclipse.exe >./ecl.log;'
  %     'rm ./ECL.CFG'];
  if isunix
    ecl_string = ['@eclipse ', basefilename,' > ./ecl.log 2>&1'];
  else
    ecl_string = ['$eclipse ', basefilename,' > ecl.log'];
  end
else
  if isunix
    ecl_string = [commandline, ' ', basefilename,' > ./ecl.log 2>&1'];
  else
    ecl_string = [commandline, ' ', basefilename, ' > ecl.log'];
  end
end

% Initialize loop
repeat_loop=1;
rep_count=0;

while repeat_loop
% $$$     %=================================
% $$$     % If on unix, we can check the
% $$$     % connection to the licecse server
% $$$     %=================================
% $$$     if isunix
% $$$       s=unix('ping -c1 -w1 octane >/dev/null');
% $$$       if s ~= 0
% $$$ 	% ping reported an error. License server is unreachable.
% $$$ 	errflag = 1;
% $$$       end
% $$$     end
 


    %=================================
    % Check status of license server
    %=================================
    if isunix
        cmd=['unset LM_LICENSE_FILE;',... % To prevent the matlab
            lm_string];                  % server from being asked...
        % [s,w]=unix(cmd);
        unix(cmd);
        f = fopen('lm.txt');
        if f<0
            errflag = 1;
        else
            % Read results file from license check
            w = fscanf(f,'%c');
            fclose(f);
            % check results
            b=strfind(w,'Users of eclipse');
            if isempty(b)
                % The string was not found in the output.
                % This means that we could not read the status of the eclipse
                % licenses from the server.
                errflag = 1;
            else
                % The string was found.
                % This means that the licence server is reachable, and that we
                % know, by counting the number in use, whether any of the licences
                % are free for us to use.
                % Still, since experience has shown that we in some cases may try
                % to start several simulations per second, we should perhaps
                % not try to start any more simulations if only one or two
                % licences are free.
                % Currently we have 9 licenses availible in total.
                % This number may be checked by evaluating
                a=find(isspace(w));c=a(a>b+27);licnum=str2num(w(c(1):c(2))); %#ok<ST2NM>
                inuse=length(strfind(w,'start'));
                if licnum-inuse <= 1
                    % Too few licenses
                    errflag = 2;
                else
                    % Start simulation
                end
            end
        end
    end

  if errflag == 0
    %=================================
    % Call simulator
    %=================================
    % Since we got this far without errors, we are ready to call
    % the simulator.
    %[s,w]=unix(ecl_string);
    unix(ecl_string);
  end

  %=================================  
  % Process error information
  %=================================
  switch errflag
   case 0
    % No errors found. Terminate loop.
    repeat_loop=0;
   case 1
    % The connection to the license server is broken.
    % It may be restored if we wait long enough.
    disp('Could not check license status.')
    % Either we wait a while and try again, or we quit,
    % depending on the values in the input variable erroraction.
    if strcmp(onerr1,'S')
      disp('Exiting')
      runstatus=-1;
      repeat_loop = 0;
    else
      repeat_loop = 1;
      rep_count = rep_count+1;
      rep_delay = user_delay;
    end
   case 2
    % The number of availible licenses is too low.
    disp('No free licenses.')
    if strcmp(onerr2,'S')
      disp('Exiting')
      runstatus=-1;
      repeat_loop = 0;
    else
      repeat_loop = 1;
      % This will probably improve faster than a broken connection.
      rep_delay = 0.25; % 1/4 minute
    end
  end

  if repeat_loop ~= 0
    if rep_count<=maxretry
      disp(['Trying again in ',num2str(rep_delay),...
		    ' minutes.'])
      pause(60*rep_delay)
      errflag = 0;
    else
      disp('Max number of retries reached. Exiting')
      runstatus = -1;
      repeat_loop=0;
    end
  end
end
