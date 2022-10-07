% This script is used to generate a test dataset for code testing and
% debugging purposes.

%% Run this section to get dataset for place cell analysis debugging
% first you have to define which session you want to look at, refer to 
% metadata excel sheets for more information
session = 1; 

% get analysis parameters, file paths, etc.
[ds, ~] = getdsap; 

% find the session metadata 
metadata = readtable(fullfile(ds.metadataPath, ds.metadataFileName_forSecondaryAnalysis));
META = table2struct(metadata(session,:));

% load the test data you want and you can start testing codes
[Inscopix, ~, Noldus, ~, Start_T, ~] = Load_RawData_For_Cellset_Pair(META, ds);

% clear up command window and remove the overall summarizing table
clc
clear metadata

%% Run this section to get dataset for engram analysis debugging

session = 1; 

% get analysis parameters, file paths, etc.
[ds, ap] = getdsap; 

% find the session metadata 
metadata = readtable(fullfile(ds.metadataPath, ds.metadataFileName_forSecondaryAnalysis));
META = table2struct(metadata(session,:));

% load the test data you want and you can start testing codes
[Inscopix1, Inscopix2, Noldus1, Noldus2, Start_T1, Start_T2] = Load_RawData_For_Cellset_Pair(META, ds); 
% load place cell information
load(fullfile(ds.metadataPath, 'Placeness_MetaData.mat'));
[PlaceCell_Summary1, PlaceCell_Summary2] = Find_Placeness_Summary_For_Cellset_Pair(META, Placeness_MetaData);

% clear up command window and remove the overall summarizing table
clc
clear metadata Placeness_MetaData
