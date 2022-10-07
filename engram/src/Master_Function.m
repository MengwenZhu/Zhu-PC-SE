% This script puts all the analysis functions into a complete workflow

% IMPORTANT: before you run anything, make sure the following:
% 1. In getdsap, set file names & file paths you want, or else things will
% run for some datasets unexpectedly and it will get really messy
% 2. Have metadata excel sheet correctly filled out

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

%% Primary Analysis
% Output: inscopix raw recording fully processed, 
% plus all the inscopix raw data csv files

% if you need to re-run primary analysis, you have to delete
% previously-generated results; to do this, go to Preprocessing and CNMFE
% Analysis folder and find Delete_results.m, and simply run that function.
runPrimaryBatchAnalysis;

% Get Recording Start Times (posix time)
% Output: a new secondary metadata excel sheet
[ds, ap]=getdsap;
getSessionPairStartTimes(ds, ap);

%% Secondary Analysis

% Place Cell Analysis
% Output: a MAT file that summarizes place cell information 
% (as defined by Mutual Information) for all recordings
Get_Place_Cell_Info_For_All_Cellsets;



% Spatial Engram Analysis (Event-Rate, Rate-Map, Population-Vector)
% Output: a MAT file that summarizes ER, RM, and PV analysis results, in
% addition to many other parameters (e.g. cell recurring probability, 
% mobiltiy, proportion of coactive cells, total number of cells active, etc.)
runEngramAnalysis;



% Cross-Day Analysis
% Output: two MAT files summarizing cross-day analysis based on PM/day2 sessions
% that use different contexts. It is used as a control and exploration of how
% memory drifts over time in calcium imaging recordings.
runXDayAnalysis;

% clear all variables (data are all saved on local drive, no need to have them in workspace)
clear