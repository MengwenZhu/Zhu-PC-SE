% % MIT License
% % 
% % Copyright (c) [2022] [Mengwen Zhu]
% % 
% % Permission is hereby granted, free of charge, to any person obtaining a copy
% % of this software and associated documentation files (the "Software"), to deal
% % in the Software without restriction, including without limitation the rights
% % to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
% % copies of the Software, and to permit persons to whom the Software is
% % furnished to do so, subject to the following conditions:
% % 
% % The above copyright notice and this permission notice shall be included in all
% % copies or substantial portions of the Software.
% % 
% % THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
% % IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
% % FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
% % AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
% % LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
% % OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
% % SOFTWARE.

%%
machine=lower(getenv('computername'));
disp(['Machine''s name is ' machine]);
%% Change according to the directory that has those files
% For instance, if you put the Engram_Analysis folder on desktop, click
% into the folder to find a folder called "engram", and right click on anything
% within "engram " and go to 'properties', and then copy the directory
% below

engram_basepath = 'C:\Users\pearcelabuser\Desktop\Engram_Analysis\matlab\engram';
addpath C:\Users\pearcelabuser\Desktop\Engram_Analysis\matlab
addpath C:\Users\pearcelabuser\Desktop\Engram_Analysis\matlab\utilities
addpath C:\Users\pearcelabuser\Desktop\Engram_Analysis\matlab\utilities\system

%% No need to change anything from below
addpath(fullfile(engram_basepath, 'src'));
addpath(fullfile(engram_basepath, 'src\Filter Utilities'));
addpath(fullfile(engram_basepath, 'src\Get File Names and Rec Times'));
addpath(fullfile(engram_basepath, 'src\userParameters'));
addpath(fullfile(engram_basepath, 'src\Make_STAT_tables'));
addpath(fullfile(engram_basepath, 'src\For_Prism_Graphing'));
addpath(fullfile(engram_basepath, 'src\For_Prism_Graphing\C57 data'));
addpath(fullfile(engram_basepath, 'src\For_Prism_Graphing\GAD data'));
addpath(fullfile(engram_basepath, 'src\For_Prism_Graphing\CCK data'));
addpath(fullfile(engram_basepath, 'src\Preprocessing and CNMFE Analysis'));
addpath(fullfile(engram_basepath, 'src\Secondary Analysis Functions'));
addpath(fullfile(engram_basepath, 'src\Graphing for Specific Dataset'));
addpath(fullfile(engram_basepath, 'src\Graphing for Specific Dataset\GAD ETOM graphics'));
addpath(fullfile(engram_basepath, 'src\Graphing for Specific Dataset\GAD Midazolam graphics'));
addpath(fullfile(engram_basepath, 'src\Graphing for Specific Dataset\C57 CPP graphics'));
addpath(fullfile(engram_basepath, 'src\Graphing for Specific Dataset\C57 4h_vs_24h graphics'));
addpath(fullfile(engram_basepath, 'src\Graphing for Specific Dataset\C57 Midazolam graphics'));
addpath(fullfile(engram_basepath, 'src\Graphing for Specific Dataset\Plotting Utilities'));

curWorkDir=fullfile(engram_basepath, 'src');
cd(curWorkDir);
disp(['Current working directory is ' curWorkDir]);