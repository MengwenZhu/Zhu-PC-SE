function [XDay_ER_Summary, XDay_Analysis_Matrix] = XDay_ER_Analysis_For_One_Animal(XDay_Event_Rates, XDay_Placeness_Summary)

% this function runs ER correlation analysis among all PM sessions or day2
% sessions of one animal and output a structure that summarizes all PV
% related data between each pair of sessions

%% First of all, compute combination ordered pairs that will be used for cross-day analysis pair assignment
% e.g. you wish to compute PV between PM1&2, PM1&3, PM1&4, etc.

Comb_Pair = nchoosek((1:size(XDay_Event_Rates,2)),2);


%% Perform cross-day analysis for all of the chosen pairs
XDay_ER_Summary = struct();

for n = 1:size(Comb_Pair, 1)
    try
        [XDay_ER_Summary(n).Date1, XDay_ER_Summary(n).Date2, XDay_ER_Summary(n).ER_Corr_Values] = XDay_ER(XDay_Event_Rates(Comb_Pair(n,1)).event_rate, XDay_Event_Rates(Comb_Pair(n,2)).event_rate, XDay_Placeness_Summary(Comb_Pair(n,1)).Date, XDay_Placeness_Summary(Comb_Pair(n,2)).Date);
    catch
        fprintf('error in a pair of cross-day analysis')
    end
end

% set errors as no memory
for n = 1:size(XDay_ER_Summary,2)
    if isempty(XDay_ER_Summary(n).ER_Corr_Values) || isnan(XDay_ER_Summary(n).ER_Corr_Values)
        XDay_ER_Summary(n).ER_Corr_Values = 10^-4;
    end
end

%% Generate cross-day matrix

XDay_Analysis_Matrix = nan(size(XDay_Event_Rates,2),size(XDay_Event_Rates,2));

% diagonals are all ones
for n = 1:size(XDay_Event_Rates, 2)
    XDay_Analysis_Matrix(n,n) = 1;
end

% assign values to the matrix
for n = 1:size(Comb_Pair, 1)
    XDay_Analysis_Matrix(Comb_Pair(n,1),Comb_Pair(n,2)) = XDay_ER_Summary(n).ER_Corr_Values;
    XDay_Analysis_Matrix(Comb_Pair(n,2),Comb_Pair(n,1)) = XDay_ER_Summary(n).ER_Corr_Values;
end

%% Graph cross-day matrix (optional)

% % heatmap(XDay_Analysis_Matrix,'CellLabelColor','none'); 
% % colormap jet; 
% % Ax = gca;
% % Ax.XDisplayLabels = nan(size(Ax.XDisplayData));
% % Ax.YDisplayLabels = nan(size(Ax.YDisplayData));
% % caxis([0, 1]);

end