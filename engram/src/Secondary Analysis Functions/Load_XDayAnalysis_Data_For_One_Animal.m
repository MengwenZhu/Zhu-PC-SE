function [XDay_Rate_Maps, XDay_Placeness_Summary, XDay_Event_Rates] = Load_XDayAnalysis_Data_For_One_Animal(animalName, ds, ap)

% Notes: this function loads raw data required by XDay-Analysis

%% Load useful XDay-Analysis raw data into structure
% read metadata
metadata = readtable(fullfile(ds.metadataPath, ds.metadataFileName_forSecondaryAnalysis));
metadata = table2struct(metadata);

% create empty structure
XDay_Raw_Data = struct('Date', [], 'Inscopix', [], 'Noldus', [], 'Start_T', [], 'Dose', [], 'Context', []);

% load raw data for all PM sessions, or day2 sessions
for n = 1:size(metadata,1)
    if strcmp(metadata(n).animalName, animalName) && metadata(n).session2_dose == 0 
        META = metadata(n);
        [XDay_Raw_Data(n).Date, XDay_Raw_Data(n).Inscopix, XDay_Raw_Data(n).Noldus, XDay_Raw_Data(n).Start_T, XDay_Raw_Data(n).Dose, XDay_Raw_Data(n).Context] = Load_RawData_For_XDay_Cellset_Pair(META, ds);
    end
end

% clean up
XDay_Raw_Data = XDay_Raw_Data(all(~cellfun(@isempty,struct2cell(XDay_Raw_Data))));

% exclude pm/day2 data where animal is injected with drug
Dose_info = extractfield(XDay_Raw_Data, 'Dose');
Dose_info = Dose_info';
XDay_Raw_Data(Dose_info~=0) = [];

% exclude repeated contexts (e.g. if context A is repeated three times, 
% then two of these experiments will be excluded from XDay analysis)
Context_info = extractfield(XDay_Raw_Data, 'Context');
Context_info = Context_info';
[~, uniqueIdx] = unique(Context_info);
XDay_Raw_Data = XDay_Raw_Data(uniqueIdx);

% sort the date of the extracted data
[~,index] = sortrows({XDay_Raw_Data.Date}.');
XDay_Raw_Data = XDay_Raw_Data(index); 
clear index
        
%% Calculate occupancy and smoothed rate maps for each day and organize into a structure (get ready for PV correlation)
% also calculate event rate

XDay_Occ_Maps = struct();
XDay_Rate_Maps = struct();
XDay_Event_Rates = struct();
for ii = 1:size(XDay_Raw_Data,2)
    [XDay_Raw_Data(ii).Behavior_Tracking, XDay_Raw_Data(ii).Cell_Events] = XDay_RawData_Preparation(XDay_Raw_Data(ii).Inscopix, XDay_Raw_Data(ii).Noldus, XDay_Raw_Data(ii).Start_T, ds, ap);
    [XDay_Event_Rates(ii).event_rate] = Compute_Event_Rate(XDay_Raw_Data(ii).Inscopix, XDay_Raw_Data(ii).Start_T, ds);
    [XDay_Occ_Maps(ii).occ_map] = Compute_Occupancy_Map(XDay_Raw_Data(ii).Behavior_Tracking, ds, ap);
    [XDay_Rate_Maps(ii).smoothed_rate_map] = Compute_Smoothed_Ca_Event_Rate_Map(XDay_Raw_Data(ii).Behavior_Tracking, XDay_Raw_Data(ii).Cell_Events, XDay_Occ_Maps(ii).occ_map, ap);
end

%% Find all the placeness_summary for relevant experiments

load(fullfile(ds.metadataPath, 'Placeness_MetaData.mat'), 'Placeness_MetaData');

% create empty structure
XDay_Placeness_Summary = struct('Date', [], 'Dose', [], 'Context', [], 'placeness', []);

for kk = 1:size(Placeness_MetaData,2)
    if strcmp(Placeness_MetaData(kk).exptParadigm, '4h_memory_test') && strcmp(Placeness_MetaData(kk).animalName, animalName) && strcmp(Placeness_MetaData(kk).session, 'PM')
        XDay_Placeness_Summary(kk).placeness = Placeness_MetaData(kk).Place_Cell_Summary;
        XDay_Placeness_Summary(kk).Dose = Placeness_MetaData(kk).session2_dose;
        XDay_Placeness_Summary(kk).Context = Placeness_MetaData(kk).session2_context;
        XDay_Placeness_Summary(kk).Date = Placeness_MetaData(kk).exptDate1;
    elseif strcmp(Placeness_MetaData(kk).exptParadigm, '24h_memory_test') && strcmp(Placeness_MetaData(kk).animalName, animalName) && strcmp(Placeness_MetaData(kk).session, 'day2')
        XDay_Placeness_Summary(kk).placeness = Placeness_MetaData(kk).Place_Cell_Summary;
        XDay_Placeness_Summary(kk).Dose = Placeness_MetaData(kk).session2_dose;
        XDay_Placeness_Summary(kk).Context = Placeness_MetaData(kk).session2_context;
        XDay_Placeness_Summary(kk).Date = Placeness_MetaData(kk).exptDate2;
    end
end

% clean up
XDay_Placeness_Summary = XDay_Placeness_Summary(~cellfun(@isempty,{XDay_Placeness_Summary.placeness}));

% exclude pm/day2 data where animal is injected with drug
Dose_info = extractfield(XDay_Placeness_Summary, 'Dose');
Dose_info = Dose_info';
XDay_Placeness_Summary(Dose_info~=0) = [];

% exclude repeated contexts (e.g. if context A is repeated three times, 
% then two of these experiments will be excluded from XDay analysis)
Context_info = extractfield(XDay_Placeness_Summary, 'Context');
Context_info = Context_info';
[~, uniqueIdx] = unique(Context_info);
XDay_Placeness_Summary = XDay_Placeness_Summary(uniqueIdx);

% sort the date of the extracted placeness data
[~,index] = sortrows({XDay_Placeness_Summary.Date}.');
XDay_Placeness_Summary = XDay_Placeness_Summary(index); 
clear index

% check if sizes of output structures match
if size(XDay_Rate_Maps,2) == size(XDay_Placeness_Summary,2)
    disp('output correct, ready to go')
else
    disp('warning: output structures have different sizes')
end

end