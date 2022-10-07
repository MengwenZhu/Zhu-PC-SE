% This function gives you a summary of cross-day analysis for all animals
% within a dataset

%% Get basic parameters & load metadata
[ds, ap] = getdsap;
metadata = readtable(fullfile(ds.metadataPath, ds.metadataFileName_forSecondaryAnalysis));
metadata = table2struct(metadata);

%% Extract all animal names into a cell
A = struct();
for n = 1:size(metadata,1)
    A(n).names = metadata(n).animalName;
end
A = struct2cell(A);
animalNames = unique(A(cellfun('isclass',A,'char')));

%% Run cross-day analysis through all animals
t_start = tic;

XDay_PV_Summary = struct();
XDay_PV_Matrix = struct();
XDay_RM_Summary = struct();
XDay_RM_Matrix = struct();
XDay_ER_Summary = struct();
XDay_ER_Matrix = struct();

% since structure field names can't have brackets, need to delete those in
% animal names
Simplified_Animal_Names = strrep(animalNames, '(', '');
Simplified_Animal_Names = strrep(Simplified_Animal_Names, ')', '');

for n = 1:size(animalNames,1)
    [XDay_Rate_Maps, XDay_Placeness_Summary, XDay_Event_Rates] = Load_XDayAnalysis_Data_For_One_Animal(string(animalNames(n)), ds, ap);
    [XDay_PV_Summary.(string(Simplified_Animal_Names(n))), XDay_PV_Matrix.(string(Simplified_Animal_Names(n)))] = XDay_PV_Analysis_For_One_Animal(XDay_Rate_Maps, XDay_Placeness_Summary, ap);
    [XDay_RM_Summary.(string(Simplified_Animal_Names(n))), XDay_RM_Matrix.(string(Simplified_Animal_Names(n)))] = XDay_RM_Analysis_For_One_Animal(XDay_Rate_Maps, XDay_Placeness_Summary, ap);
    [XDay_ER_Summary.(string(Simplified_Animal_Names(n))), XDay_ER_Matrix.(string(Simplified_Animal_Names(n)))] = XDay_ER_Analysis_For_One_Animal(XDay_Event_Rates, XDay_Placeness_Summary);
end

% save workspace
save(fullfile(ds.metadataPath, 'XDay Data Summary.mat'), '-v7.3', 'XDay_PV_Summary', 'XDay_PV_Matrix', 'XDay_RM_Summary', 'XDay_RM_Matrix', 'XDay_ER_Summary', 'XDay_ER_Matrix')
% report
t_stop = toc(t_start);
disp("Cross-Day analysis took " + string(t_stop) + " s")
% clear all variables
clear
