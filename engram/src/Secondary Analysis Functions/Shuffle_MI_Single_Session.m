function [Session_Shuffled_Spatial_Info]=Shuffle_MI_Single_Session(Session_Spatial_Info, Session_RawData, Corrected_Occupancy_Map, RawXYData, ds, ap)
%% Shuffle single sessin calcium event times

% Plan A: Shuffle by shifting all calcium events by a random interval of time
Rand_Time=ds.Total_T_Recording*rand; % generate a random time between 0 & 600 sec

Session_RawData(:,1)=Session_RawData(:,1)+Rand_Time; % shuffle calcium event times by shifting all events by the same amount

for n=1:size(Session_RawData,1)
    if Session_RawData(n,1)>ds.Total_T_Recording
        Session_RawData(n,1)=Session_RawData(n,1)-ds.Total_T_Recording; % if event time goes beyond 600 then shift it back to beginning
    else
        Session_RawData(n,1)=Session_RawData(n,1); % if event time stays within 600 keep it same
    end
end

%% same procedures as a normal calculation of MI

[Session_RawData]=Find_XY_Position_Single_Session(Session_RawData, RawXYData, ds.Noldus_Sampling_Frequency);

for n=1:size(Session_RawData,1)
    if Session_RawData(n,6)<ap.Low_Mobility_Threshold || Session_RawData(n,6)>ap.High_Mobility_Threshold
        Session_RawData(n,1:6)=10^4;
    end
end
Session_RawData(Session_RawData==10^4)=[];
Session_RawData=reshape(Session_RawData,[],6);



% set XY boundaries of Ca event maps
n=(unique(Session_RawData(:,2)))';
MaxFiring=max(sum(Session_RawData(:,2)==n));

A=min(RawXYData(:,2),[],'omitnan');
B=max(RawXYData(:,2),[],'omitnan');
C=min(RawXYData(:,3),[],'omitnan');
D=max(RawXYData(:,3),[],'omitnan');

Map_LowX_Bound=fix(A)-ap.Arena_Add_Edge; % To fully depict the arena within a matrix without leaving out data points
Map_HighX_Bound=fix(B)+ap.Arena_Add_Edge;
Map_LowY_Bound=fix(C)-ap.Arena_Add_Edge;
Map_HighY_Bound=fix(D)+ap.Arena_Add_Edge;



% same thing generating Ca event maps
RawFiring3DMatrix=nan(MaxFiring, size(Session_RawData,2), max(unique(Session_RawData(:,2))));
for n=1:max(unique(Session_RawData(:,2)))
    y=sum(Session_RawData(:,2)==n);
    RawFiring3DMatrix(:,:,n)=reshape([Session_RawData(cat(1,Session_RawData(:,2))==n,:);NaN(MaxFiring-y,size(Session_RawData,2))],MaxFiring,size(Session_RawData,2));
end
RawFiring3DMatrix(RawFiring3DMatrix==0)=nan;


Ca_Event_Map=nan(size(Corrected_Occupancy_Map,1),size(Corrected_Occupancy_Map,2),size(RawFiring3DMatrix,3));
for k=1:size(RawFiring3DMatrix,3)
    Ca_Event_Map(:,:,k)=hist3(RawFiring3DMatrix(:,4:5,k),'ctrs',{Map_LowX_Bound:(Map_HighX_Bound-Map_LowX_Bound)/ap.Map_Division_For_MI:Map_HighX_Bound  Map_LowY_Bound:(Map_HighY_Bound-Map_LowY_Bound)/ap.Map_Division_For_MI:Map_HighY_Bound}); 
end
Ca_Event_Map=rot90(Ca_Event_Map,1);

Ca_Event_Rate_Map=nan(size(Corrected_Occupancy_Map,1),size(Corrected_Occupancy_Map,2),size(RawFiring3DMatrix,3));
for k=1:size(RawFiring3DMatrix,3)
    Ca_Event_Rate_Map(:,:,k)=(Ca_Event_Map(:,:,k))./Corrected_Occupancy_Map;
end
Ca_Event_Rate_Map(isnan(Ca_Event_Rate_Map)) = 0;

Normalized_Ca_Event_Rate_Map=nan(size(Corrected_Occupancy_Map,1),size(Corrected_Occupancy_Map,2),size(RawFiring3DMatrix,3));
for k=1:size(RawFiring3DMatrix,3)
    Normalized_Ca_Event_Rate_Map(:,:,k)=(Ca_Event_Rate_Map(:,:,k))./sum(Ca_Event_Rate_Map(:,:,k),'all');
end


Session_Shuffled_Spatial_Info=zeros(size(Ca_Event_Map,3),1);
for n = (unique(Session_RawData(:,2)))'
    [Session_Shuffled_Spatial_Info(n,1)] = calc_mutual_information(Normalized_Ca_Event_Rate_Map(:,:,n), Corrected_Occupancy_Map(:,:,1));
end
Session_Shuffled_Spatial_Info(Session_Shuffled_Spatial_Info==0)=nan; %get rid of useless info (if cell has MI of 0 then useless anyway)



% fix bugs: since this is shuffling, some cells might be excluded because
% their events are all shuffled to immobility periods, which makes it
% difficult for later concatenation of arrays in calculating MI p-values,
% so keep the shuffled data same size as the real data

d=size(Session_Spatial_Info,1)-size(Session_Shuffled_Spatial_Info,1);

if size(Session_Shuffled_Spatial_Info,1)>size(Session_Spatial_Info,1) % if for some reason the shuffled MI array size is bigger
    Session_Shuffled_Spatial_Info(size(Session_Spatial_Info,1)+1:end,1)=100000; % mark the last parts of the shuffled array as some number that is impossible for MI values in our situation 
end
Session_Shuffled_Spatial_Info(Session_Shuffled_Spatial_Info==100000)=[]; % eliminate the parts making the shuffled array bigger


if size(Session_Shuffled_Spatial_Info,1)<size(Session_Spatial_Info,1) % if for some reason the shuffled MI array is smaller than actual data
    Session_Shuffled_Spatial_Info(end+1:end+d,1)=nan; % add some NANs to the end of the shuffled array
end


end