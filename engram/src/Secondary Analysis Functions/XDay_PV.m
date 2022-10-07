function [Date1, Date2, PV_Corr_Values, PV_Corr_Memory_Index, PV_Corr_Coherent_Rot, PV_Proportion, Shuffled_PV_Proportion, Shuffled_PV_Mean] = XDay_PV(Rate_Map1, Rate_Map2, PlaceCell_Summary1, PlaceCell_Summary2, firstDate, secondDate, ap)

% this function gives PV related outputs for a chosen pair of sessions to
% be included in cross-day analysis

Date1 = firstDate;
Date2 = secondDate;

%% Make sure that place cell matrices have the same size
p = abs(size(PlaceCell_Summary2,1) - size(PlaceCell_Summary1,1));
if  p ~= 1 && size(PlaceCell_Summary2,1) > size(PlaceCell_Summary1,1)
    PlaceCell_Summary1(end+1:size(PlaceCell_Summary2,1)) = nan;
elseif p ~= 1 && size(PlaceCell_Summary2,1) < size(PlaceCell_Summary1,1)
    PlaceCell_Summary2(end+1:size(PlaceCell_Summary1,1)) = nan;
end
if p == 1 && size(PlaceCell_Summary2,1) > size(PlaceCell_Summary1,1)
    PlaceCell_Summary1(end+1) = nan;
elseif p == 1 && size(PlaceCell_Summary1,1) > size(PlaceCell_Summary2,1)
    PlaceCell_Summary2(end+1) = nan;
end

% Compute # of Place Cells in Both Sessions
K = zeros(size(PlaceCell_Summary1,1),1);
for n = 1:size(PlaceCell_Summary1,1)
    if ~isnan(PlaceCell_Summary1(n,1)) && ~isnan(PlaceCell_Summary2(n,1))
        K(n,1) = (PlaceCell_Summary1(n,1) == 1 || PlaceCell_Summary2(n,1) == 1);
    end
end
Coactive_Num_Place_Cell = sum(K);

%% Find coative cells

A = nan(size(PlaceCell_Summary1,1),1);
for n = 1:size(PlaceCell_Summary1,1)
    if ~isnan(PlaceCell_Summary1(n,1)) && ~isnan(PlaceCell_Summary2(n,1))
        A(n,1) = 1;
    end
end

Coactive_Cell_Identity = nan(size(A,1),1);
for n = 1:size(A,1)
    if A(n,1) == 1
        Coactive_Cell_Identity(n,1) = n;
    end
end
Coactive_Cell_Identity(isnan(Coactive_Cell_Identity)) = [];

%% Make sure rate map matrices have the same size

p = abs(size(Rate_Map2,3) - size(Rate_Map1,3));
if  p ~= 1 && size(Rate_Map2,3) > size(Rate_Map1,3)
    Rate_Map1(1:size(Rate_Map2,1),1:size(Rate_Map2,2),end+1:size(Rate_Map2,3)) = zeros(size(Rate_Map2,1),size(Rate_Map2,2),size(Rate_Map2,3)-size(Rate_Map1,3));
elseif p ~= 1 && size(Rate_Map2,3) < size(Rate_Map1,3)
    Rate_Map1(1:size(Rate_Map2,1),1:size(Rate_Map2,2),end+1:size(Rate_Map1,3)) = zeros(size(Rate_Map2,1),size(Rate_Map2,2),size(Rate_Map1,3)-size(Rate_Map2,3));
end
if p == 1 && size(Rate_Map2,3) > size(Rate_Map1,3)
    Rate_Map1(:,:,end+1) = zeros(size(Rate_Map2,1), size(Rate_Map2,2));
elseif p == 1 && size(Rate_Map1,3) > size(Rate_Map2,3)
    Rate_Map2(:,:,end+1) = zeros(size(Rate_Map2,1), size(Rate_Map2,2));
end

%% Select coactive cells' smoothed rate maps

for n = 1:size(Rate_Map1,3)
    if ~(n == Coactive_Cell_Identity)
        Rate_Map1(1:size(Rate_Map1,1), 1:size(Rate_Map1,2),n) = zeros(size(Rate_Map1,1), size(Rate_Map1,2));
        Rate_Map2(1:size(Rate_Map2,1), 1:size(Rate_Map2,2),n) = zeros(size(Rate_Map2,1), size(Rate_Map2,2));
    end
end

%% PV correlation analysis

if sum(Rate_Map1(:),'all') ~= 0 || sum(Rate_Map2(:), 'all') ~= 0
    [Shuffled_PV_Corr, Shuffled_PV_Mean]=Shuffle_PV(Coactive_Num_Place_Cell, Rate_Map1, Rate_Map2, PlaceCell_Summary1, PlaceCell_Summary2, ap);
    [PV_Corr_Values,PV_Corr_Memory_Index,PV_Corr_Coherent_Rot,PV_Proportion,Shuffled_PV_Proportion]=PV_Correlation_Analysis(Coactive_Num_Place_Cell, Shuffled_PV_Corr, Rate_Map1, Rate_Map2, PlaceCell_Summary1, PlaceCell_Summary2, ap);
elseif sum(Rate_Map1(:),'all') == 0 || sum(Rate_Map2(:), 'all') == 0
    disp('no coactive cells for two sessions');
    PV_Corr_Memory_Index = nan;
    PV_Corr_Coherent_Rot = nan;
    PV_Corr_Values = [];
    PV_Proportion = []; 
    Shuffled_PV_Proportion = [];
    Shuffled_PV_Mean = [];
end

end