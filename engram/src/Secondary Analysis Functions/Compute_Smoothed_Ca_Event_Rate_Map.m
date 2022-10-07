function [Smoothed_Ca_Event_Rate_Map] = Compute_Smoothed_Ca_Event_Rate_Map(Noldus, Cell_Events, Corrected_Occupancy_Map, ap)

% This function gives you a smoothed calcium event rate map

%% Calculate Ca_Event_Map

% Calculate the boundaries that will be used
A=min(Noldus(:,2),[],'omitnan');
B=max(Noldus(:,2),[],'omitnan');
C=min(Noldus(:,3),[],'omitnan');
D=max(Noldus(:,3),[],'omitnan');

Map_LowX_Bound=fix(A)-ap.Arena_Add_Edge; 
Map_HighX_Bound=fix(B)+ap.Arena_Add_Edge;
Map_LowY_Bound=fix(C)-ap.Arena_Add_Edge;
Map_HighY_Bound=fix(D)+ap.Arena_Add_Edge;

n=(unique(Cell_Events(:,2)))';
Max_Firing=max(sum(Cell_Events(:,2)==n)); 

% Reorganize inscopix data
Reorganize=nan(Max_Firing, size(Cell_Events,2), size(unique(Cell_Events(:,2)),1));
for n=(unique(Cell_Events(:,2)))'
    y=sum(Cell_Events(:,2)==n);
    Reorganize(:,:,n)=reshape([Cell_Events(cat(1,Cell_Events(:,2))==n,:);NaN(Max_Firing-y,size(Cell_Events,2))],Max_Firing,size(Cell_Events,2));
end

% Calculate event map
Ca_Event_Map=nan(size(Corrected_Occupancy_Map,1),size(Corrected_Occupancy_Map,2),size(Reorganize,3));
for k=1:size(Reorganize,3)
    Ca_Event_Map(:,:,k)=hist3(Reorganize(:,4:5,k),'ctrs',{Map_LowX_Bound:(Map_HighX_Bound-Map_LowX_Bound)/ap.Map_Division_For_PV:Map_HighX_Bound  Map_LowY_Bound:(Map_HighY_Bound-Map_LowY_Bound)/ap.Map_Division_For_PV:Map_HighY_Bound});
end
Ca_Event_Map=rot90(Ca_Event_Map,1);

%% Then calculate Ca_Event_Rate_Map
Ca_Event_Rate_Map=Ca_Event_Map./Corrected_Occupancy_Map;
Ca_Event_Rate_Map(isnan(Ca_Event_Rate_Map))=0;
Ca_Event_Rate_Map=reshape(Ca_Event_Rate_Map,size(Corrected_Occupancy_Map,1),size(Corrected_Occupancy_Map,2),[]);

%% Lastly, calculate Smoothed_Ca_Event_Rate_Map

Smoothed_Ca_Event_Rate_Map=zeros(size(Ca_Event_Rate_Map,1), size(Ca_Event_Rate_Map,2), size(Ca_Event_Rate_Map,3));
for n=1:size(Ca_Event_Rate_Map,3)
    Smoothed_Ca_Event_Rate_Map(:,:,n)=filter2(ap.Gaussian_Filter,Ca_Event_Rate_Map(:,:,n)); 
end

end