
function [colMap] = colMapGen(topCol,botCol,numCol,varargin)
    % Creates a colormap using two boundary colors and one middle
    % color. Both boundary colors blend into the middle color.
    % By default, the middle color is white. This can be changed using
    % the 'midCol' name-value pair argument. Input colors must be in RGB
    % triplet format (e.g. [0 0 0] for black). The user defines the number
    % of colors (segments) to make up the colormap (i.e from the first 
    % boundary color to the second boundary color).
    
    %%%% Example usage %%%%
    
    % colMap = colMapGen([1 1 0],[0 0 0],20);
        % creates a colormap going from yellow to black, with a white
        % center. Colormap consists of 20 segments going from yellow to
        % black.
       
    % colMap = colMapGen([0 1 1],[0 0 1],100,'midCol',[0 0 0]);
        % creates a colormap going from red to blue in 100 segments. The
        % middle color is set to black.
        
    % Once a colormap is generated, upload to current figure using:
    % colormap(gca,colMap)
        
% Copyright (c) 2019, Timothy Olsen
% All rights reserved.
% 
% Redistribution and use in source and binary forms, with or without
% modification, are permitted provided that the following conditions are met:
% 
% * Redistributions of source code must retain the above copyright notice, this
%   list of conditions and the following disclaimer.
% 
% * Redistributions in binary form must reproduce the above copyright notice,
%   this list of conditions and the following disclaimer in the documentation
%   and/or other materials provided with the distribution
% * Neither the name of UCSF nor the names of its
%   contributors may be used to endorse or promote products derived from this
%   software without specific prior written permission.
% THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
% AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
% IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
% DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS BE LIABLE
% FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
% DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
% SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
% CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY,
% OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
% OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

    % parse inputs - all inputs required except midCol, which is a
    % name-value argument
    p = inputParser; 
    addRequired(p,'topCol',@isnumeric); % required function input
    addRequired(p,'botCol',@isnumeric); % required function input
    addRequired(p,'numCol',@isnumeric); % required function input
    addParameter(p,'midCol',[1 1 1],@isnumeric); % varargin input  
    parse(p,topCol,botCol,numCol,varargin{:}); % parse inputs
    topCol  = p.Results.topCol; % define function variable from inputs
    botCol  = p.Results.botCol; % define function variable from inputs
    numCol  = p.Results.numCol; % define function variable from inputs
    midCol  = p.Results.midCol; % define function variable from inputs
    %%%% creates the upper portion of the colormap
    col1_topCol = linspace(topCol(1),midCol(1),round(numCol/2,0)); % 1st RGB triplet of topCol, from topCol(1) to midCol(1)
    col2_topCol = linspace(topCol(2),midCol(2),round(numCol/2,0)); % 2nd RGB triplet of topCol, from topCol(2) to midCol(2)
    col3_topCol = linspace(topCol(3),midCol(3),round(numCol/2,0)); % 3rd RGB triplet of topCol, from topCol(3) to midCol(3)
    RGB1 = [col1_topCol',col2_topCol',col3_topCol']; % matrix of RGB triplets (in rows), going from topCol to midCol 
    %%%% creates the lower portion of the colormap
    col1_botCol  = linspace(botCol(1),midCol(1),round(numCol/2,0)); % 1st RGB triplet of botCol, from botCol(1) to midCol(1)
    col2_botCol  = linspace(botCol(2),midCol(2),round(numCol/2,0)); % 2nd RGB triplet of botCol, from botCol(2) to midCol(2)
    col3_botCol  = linspace(botCol(3),midCol(3),round(numCol/2,0)); % 3rd RGB triplet of botCol, from botCol(3) to midCol(3)
    RGB2 = [col1_botCol',col2_botCol',col3_botCol']; % matrix of RGB triplets (in rows), going from botCol to midCol 
    %%%% Final colormap
    colMap = [RGB2;flipud(RGB1)]; % creates the final colMap by concatenating the bottom colormap with the flipped top colormap 
    
end
    
