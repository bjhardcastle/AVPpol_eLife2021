function varargout = git(varargin)
% GIT Execute a git command.
%
% GIT <ARGS>, when executed in command style, executes the git command and
% displays the git outputs at the MATLAB console.
%
% STATUS = GIT(ARG1, ARG2,...), when executed in functional style, executes
% the git command and returns the output status STATUS.
%
% [STATUS, CMDOUT] = GIT(ARG1, ARG2,...), when executed in functional
% style, executes the git command and returns the output status STATUS and
% the git output CMDOUT.

% Check output arguments.
nargoutchk(0,2)

% original version (doesn't handle space in filepath):
%{
% Specify the location of the git executable.
gitexepath = "C:\Program Files\Git\bin\git.exe";

% Construct the git command.
cmdstr = strjoin([gitexepath, varargin]);

% Execute the git command.
[status, cmdout] = system(cmdstr);
%}

% messages must be enclosed in quotation marks to work correctly with
% 'system' and git. If we forget to enter them in the command line just
% insert them here:
if any(strcmp(varargin,'-m'))
    msgIdx = find( strcmp(varargin,'-m')) + 1;
    if ~strcmp(varargin{msgIdx}(1),'"')
        varargin{msgIdx} = ['"' varargin{msgIdx}];
    end
    if ~strcmp(varargin{msgIdx}(end),'"')
        varargin{msgIdx} = [varargin{msgIdx} '"'];
    end
end

% Execute the command (fixed to handle space in path):
[status, cmdout] = system(['"C:\Program Files\Git\bin\git.exe" ' strjoin(varargin)]);

switch nargout
    case 0
        disp(cmdout)
    case 1
        varargout{1} = status;
    case 2
        varargout{1} = status;
        varargout{2} = cmdout;
end