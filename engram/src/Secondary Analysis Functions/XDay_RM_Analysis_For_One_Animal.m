function [XDay_RM_Summary, XDay_Analysis_Matrix] = XDay_RM_Analysis_For_One_Animal(XDay_Rate_Maps, XDay_Placeness_Summary, ap)

% this function runs RM correlation analysis among all PM sessions or day2
% sessions of one animal and output a structure that summarizes all RM
% related data between each pair of sessions

%% First of all, compute combination ordered pairs that will be used for cross-day analysis pair assignment
% e.g. you wish to compute RM between PM1&2, PM1&3, PM1&4, etc.

Comb_Pair = nchoosek((1:size(XDay_Rate_Maps,2)),2);


%% Perform cross-day analysis for all of the chosen pairs
XDay_RM_Summary = struct();

for n = 1:size(Comb_Pair, 1)
    try
        [XDay_RM_Summary(n).Date1, XDay_RM_Summary(n).Date2, XDay_RM_Summary(n).RM_Corr_Values, XDay_RM_Summary(n).Memory_Index, XDay_RM_Summary(n).Coherent_Rot, XDay_RM_Summary(n).RM_Proportion, XDay_RM_Summary(n).Shuffled_RM_Proportion] = XDay_RM(XDay_Rate_Maps(Comb_Pair(n,1)).smoothed_rate_map, XDay_Rate_Maps(Comb_Pair(n,2)).smoothed_rate_map, XDay_Placeness_Summary(Comb_Pair(n,1)).placeness, XDay_Placeness_Summary(Comb_Pair(n,2)).placeness, XDay_Placeness_Summary(Comb_Pair(n,1)).Date, XDay_Placeness_Summary(Comb_Pair(n,2)).Date, ap);
    catch
        fprintf('error in a pair of cross-day analysis')
    end
end

% set errors as no memory
for n = 1:size(XDay_RM_Summary,2)
    if isempty(XDay_RM_Summary(n).Memory_Index) || isnan(XDay_RM_Summary(n).Memory_Index)
        XDay_RM_Summary(n).Memory_Index = 10^-4;
    end
end

%% Generate cross-day matrix

XDay_Analysis_Matrix = nan(size(XDay_Rate_Maps,2),size(XDay_Rate_Maps,2));

% diagonals are all ones
for n = 1:size(XDay_Rate_Maps, 2)
    XDay_Analysis_Matrix(n,n) = 1;
end

% assign values to the matrix
for n = 1:size(Comb_Pair, 1)
    XDay_Analysis_Matrix(Comb_Pair(n,1),Comb_Pair(n,2)) = XDay_RM_Summary(n).Memory_Index;
    XDay_Analysis_Matrix(Comb_Pair(n,2),Comb_Pair(n,1)) = XDay_RM_Summary(n).Memory_Index;
end

%% Graph cross-day matrix (optional)

% % heatmap(XDay_Analysis_Matrix,'CellLabelColor','none'); 
% % colormap jet; 
% % Ax = gca;
% % Ax.XDisplayLabels = nan(size(Ax.XDisplayData));
% % Ax.YDisplayLabels = nan(size(Ax.YDisplayData));
% % caxis([0, 1]);


end