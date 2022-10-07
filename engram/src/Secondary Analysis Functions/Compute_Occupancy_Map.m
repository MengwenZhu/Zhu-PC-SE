function [Corrected_Occupancy_Map] = Compute_Occupancy_Map(Noldus, ds, ap)

% This function gives you a corrected occupancy map from any noldus data

% Calculate & extend occupancy map boundaries
A=min(Noldus(:,2),[],'omitnan');
B=max(Noldus(:,2),[],'omitnan');
C=min(Noldus(:,3),[],'omitnan');
D=max(Noldus(:,3),[],'omitnan');

% extend the arena size to make sure all behavioral tracking points fall
% within occupancy matrix
Map1_LowX_Bound=fix(A)-ap.Arena_Add_Edge; 
Map1_HighX_Bound=fix(B)+ap.Arena_Add_Edge;
Map1_LowY_Bound=fix(C)-ap.Arena_Add_Edge;
Map1_HighY_Bound=fix(D)+ap.Arena_Add_Edge;

% Compute occupancy map for session2 reference data
Occupancy_Map=hist3(Noldus(:,2:3),'ctrs',{Map1_LowX_Bound:(Map1_HighX_Bound-Map1_LowX_Bound)/ap.Map_Division_For_PV:Map1_HighX_Bound  Map1_LowY_Bound:(Map1_HighY_Bound-Map1_LowY_Bound)/ap.Map_Division_For_PV:Map1_HighY_Bound}); 
Occupancy_Map=rot90(Occupancy_Map,1); 
Corrected_Occupancy_Map=Occupancy_Map.*(1/ds.Noldus_Sampling_Frequency);

end