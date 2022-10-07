function [Session1_Memory_Development_Memory_Index_Shuffle, Session1_Memory_Development_Memory_Index, Session2_Memory_Development_Memory_Index_Shuffle, Session2_Memory_Development_Memory_Index] = Memory_Development_Analysis(Behavior_Tracking1, Behavior_Tracking2, Cell_Events1, Cell_Events2, PlaceCell_Summary1, PlaceCell_Summary2, Num_Place_Cell, ds, ap)
%% Step 1: Process Session2 Reference Inscopix & Noldus Data

% Extract Session2 Inscopix & Noldus Data that are Used for Memory Formation Reference
[Behavior_Tracking2_Reference, Cell_Events2_Reference] = Segment_NoldusAndInscopix(Behavior_Tracking2, Cell_Events2, ap.Session2_T_Low_Cutoff, ap.Session2_T_High_Cutoff);

% Compute Occupancy Map for Session2 Noldus Reference Data
[Corrected_Occupancy_Map_Session2_Reference] = Compute_Occupancy_Map(Behavior_Tracking2_Reference, ds, ap);

% Compute Smoothed_Ca_Event_Rate_Map for Session2 Reference Inscopix Data
[Session2_Reference_Smoothed_Ca_Event_Rate_Map] = Compute_Smoothed_Ca_Event_Rate_Map(Behavior_Tracking2_Reference, Cell_Events2_Reference, Corrected_Occupancy_Map_Session2_Reference, ap);

%% Step 2: Process Session1 Reference Inscopix & Noldus Data

[Behavior_Tracking1_Reference, Cell_Events1_Reference] = Segment_NoldusAndInscopix(Behavior_Tracking1, Cell_Events1, ap.Session1_T_Low_Cutoff, ap.Session1_T_High_Cutoff);

[Corrected_Occupancy_Map_Session1_Reference] = Compute_Occupancy_Map(Behavior_Tracking1_Reference, ds, ap);

[Session1_Reference_Smoothed_Ca_Event_Rate_Map] = Compute_Smoothed_Ca_Event_Rate_Map(Behavior_Tracking1_Reference, Cell_Events1_Reference, Corrected_Occupancy_Map_Session1_Reference, ap);

%% Step 3: Segment Session1 Data and Calculate PV Correlation between Each Section & Session2 Reference

% Create a numeric matrix that defines the temporal segments
Time_Table = zeros(1, ap.Temporal_Division+1);
for n = 1:(size(Time_Table,2)-1)
    Time_Table(1,n+1) = Time_Table(1,n) + (ds.Total_T_Recording/ap.Temporal_Division);
end

% Segmented Inscopix & Noldus Data are stored within a structure
Session1_Segmented_Data = struct();
for ii = 1:(ap.Temporal_Division)
    [Session1_Segmented_Data(ii).Noldus, Session1_Segmented_Data(ii).Inscopix] = Segment_NoldusAndInscopix(Behavior_Tracking1, Cell_Events1, Time_Table(1,ii), Time_Table(1,ii+1));
end

% Calculate corrected occupancy map and organize into same structure
for jj = 1:size(Session1_Segmented_Data, 2)
    [Session1_Segmented_Data(jj).Corrected_Occupancy_Map] = Compute_Occupancy_Map(Session1_Segmented_Data(jj).Noldus, ds, ap);
end

% Calculate Smoothed_Ca_Event_Rate_Map and organize into same structure
for zz = 1:size(Session1_Segmented_Data, 2)
    [Session1_Segmented_Data(zz).Smoothed_Ca_Event_Rate_Map] = Compute_Smoothed_Ca_Event_Rate_Map(Session1_Segmented_Data(zz).Noldus, Session1_Segmented_Data(zz).Inscopix, Session1_Segmented_Data(zz).Corrected_Occupancy_Map, ap);
end

% Calculate Shuffled PV between each time segment and session1 reference
Session1_Memory_Development_Memory_Index_Shuffle = [];
for yy = 1:size(Session1_Segmented_Data, 2)
    [Session1_Segmented_Data(yy).Shuffled_PV_Corr, Session1_Memory_Development_Memory_Index_Shuffle(yy,1)] = Shuffle_PV(Num_Place_Cell, Session1_Segmented_Data(yy).Smoothed_Ca_Event_Rate_Map, Session1_Reference_Smoothed_Ca_Event_Rate_Map, PlaceCell_Summary1, PlaceCell_Summary2, ap);
end

% Calculate PV correlation between each time segment and session1 reference
Session1_Memory_Development_Memory_Index = [];
for qq  = 1:size(Session1_Segmented_Data, 2)
    [~, Session1_Memory_Development_Memory_Index(qq,1), ~, ~, ~]=PV_Correlation_Analysis(Num_Place_Cell, Session1_Segmented_Data(qq).Shuffled_PV_Corr, Session1_Segmented_Data(qq).Smoothed_Ca_Event_Rate_Map, Session1_Reference_Smoothed_Ca_Event_Rate_Map, PlaceCell_Summary1, PlaceCell_Summary2, ap);
end

%% Step 4: Segment Session2 Data and Calculate PV Correlation between Each Section & Session2 Reference

% Segmented Inscopix & Noldus Data are stored within a structure
Session2_Segmented_Data = struct();
for ii = 1:(ap.Temporal_Division)
    [Session2_Segmented_Data(ii).Noldus, Session2_Segmented_Data(ii).Inscopix] = Segment_NoldusAndInscopix(Behavior_Tracking2, Cell_Events2, Time_Table(1,ii), Time_Table(1,ii+1));
end

% Calculate corrected occupancy map and organize into same structure
for jj = 1:size(Session2_Segmented_Data, 2)
    [Session2_Segmented_Data(jj).Corrected_Occupancy_Map] = Compute_Occupancy_Map(Session2_Segmented_Data(jj).Noldus, ds, ap);
end

% Calculate Smoothed_Ca_Event_Rate_Map and organize into same structure
for zz = 1:size(Session2_Segmented_Data, 2)
    [Session2_Segmented_Data(zz).Smoothed_Ca_Event_Rate_Map] = Compute_Smoothed_Ca_Event_Rate_Map(Session2_Segmented_Data(zz).Noldus, Session2_Segmented_Data(zz).Inscopix, Session2_Segmented_Data(zz).Corrected_Occupancy_Map, ap);
end

% Calculate Shuffled PV between each time segment and session2 reference
Session2_Memory_Development_Memory_Index_Shuffle = [];
for yy = 1:size(Session2_Segmented_Data, 2)
    [Session2_Segmented_Data(yy).Shuffled_PV_Corr, Session2_Memory_Development_Memory_Index_Shuffle(yy,1)] = Shuffle_PV(Num_Place_Cell, Session2_Segmented_Data(yy).Smoothed_Ca_Event_Rate_Map, Session2_Reference_Smoothed_Ca_Event_Rate_Map, PlaceCell_Summary1, PlaceCell_Summary2, ap);
end

% Calculate PV correlation between each time segment and session2 reference
Session2_Memory_Development_Memory_Index = [];
for qq  = 1:size(Session2_Segmented_Data, 2)
    [~, Session2_Memory_Development_Memory_Index(qq,1), ~, ~, ~]=PV_Correlation_Analysis(Num_Place_Cell, Session2_Memory_Development_Memory_Index_Shuffle(qq,1), Session2_Segmented_Data(qq).Smoothed_Ca_Event_Rate_Map, Session2_Reference_Smoothed_Ca_Event_Rate_Map, PlaceCell_Summary1, PlaceCell_Summary2, ap);
end

end