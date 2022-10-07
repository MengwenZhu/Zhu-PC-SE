function Data_Summary = Mobility_and_ER_Analysis(Inscopix1, Inscopix2, Noldus1, Noldus2, Start_T1, Start_T2, PlaceCell_Summary1, PlaceCell_Summary2, ds, ap)
% This function gives more detailed analysis of the relationship between
% animal's behavioral mobility and calcium event rate (ER), and it also
% explores what happens during immobility phase of the experiment

%% Step 1: Organize Inscopix and Noldus Files & Extract Some Basic Parameters

% Extract mouse's positions as a function of time from Noldus files and
% name them as XY1 & XY2
Behavior_Tracking1=Noldus1(1:(ds.Total_T_Recording/(1/ds.Noldus_Sampling_Frequency))+1,2:4); 
Behavior_Tracking2=Noldus2(1:(ds.Total_T_Recording/(1/ds.Noldus_Sampling_Frequency))+1,2:4); 

% Shift the cell identities in Inscopix files up by one to avoid errors (sometimes cell names start with zero)
Cell_Events1 = Inscopix1;
Cell_Events1(:,2) = Cell_Events1(:,2)+1;
Cell_Events2 = Inscopix2;
Cell_Events2(:,2) = Cell_Events2(:,2)+1;

% Correct cell events time stamps to a number between 0sec and ds.Total_T_Recording
Cell_Events1(:,1) = Cell_Events1(:,1)-Start_T1;
Cell_Events2(:,1) = Cell_Events2(:,1)-Start_T2;

% adjust the amplitude of calcium events by a factor
Cell_Events1(:,3) = ap.Amp_factor_PV*Cell_Events1(:,3);
Cell_Events2(:,3) = ap.Amp_factor_PV*Cell_Events2(:,3);

% Check if any timestamp is outside of what we expect
if any(Cell_Events1(:,1)<0) || any(Cell_Events1(:,1)>ds.Total_T_Recording) || any(Cell_Events2(:,1)<0) || any(Cell_Events2(:,1)>ds.Total_T_Recording)
    disp('corrected timestamps are outside of expected time range')
else
    disp('corrected timestamps are within expected time range, good to go')
end

%% Step 2: Explore Relationship between ER & Place Cell (do place cells have different ER than non-place cells?)

% compute the firing rates of cells detected in each recording
Session1_Event_Rates=nan(max(unique(Cell_Events1(:,2))),1);
for n=(unique(Cell_Events1(:,2)))'
    Session1_Event_Rates(n,1)=sum(Cell_Events1(:,2)==n)/ds.Total_T_Recording; 
end

Session2_Event_Rates=nan(max(unique(Cell_Events2(:,2))),1);
for n=(unique(Cell_Events2(:,2)))'
    Session2_Event_Rates(n,1)=sum(Cell_Events2(:,2)==n)/ds.Total_T_Recording; 
end

% Preallocation
PlaceCell_session1_ER = nan(size(PlaceCell_Summary1, 1),1);
Non_PlaceCell_session1_ER = nan(size(PlaceCell_Summary1, 1),1);
PlaceCell_session2_ER = nan(size(PlaceCell_Summary2, 1),1);
Non_PlaceCell_session2_ER = nan(size(PlaceCell_Summary2, 1),1);
% Find out ER for place cells & non-place cells in both sessions
for n = 1:size(PlaceCell_Summary1, 1)
    if PlaceCell_Summary1(n,1) == 1
        PlaceCell_session1_ER(n,1) = Session1_Event_Rates(n,1);
    elseif PlaceCell_Summary1(n,1) == 0
        Non_PlaceCell_session1_ER(n,1) = Session1_Event_Rates(n,1);
    end
end
PlaceCell_session1_ER(isnan(PlaceCell_session1_ER)) = [];
Non_PlaceCell_session1_ER(isnan(Non_PlaceCell_session1_ER)) = [];

for n = 1:size(PlaceCell_Summary2, 1)
    if PlaceCell_Summary2(n,1) == 1
        PlaceCell_session2_ER(n,1) = Session2_Event_Rates(n,1);
    elseif PlaceCell_Summary2(n,1) == 0
        Non_PlaceCell_session2_ER(n,1) = Session2_Event_Rates(n,1);
    end
end
PlaceCell_session2_ER(isnan(PlaceCell_session2_ER)) = [];
Non_PlaceCell_session2_ER(isnan(Non_PlaceCell_session2_ER)) = [];

% summarize place cell & non place cell ER information into arrays and
% calculate the mean ERs for both groups of cells
PlaceCellER_Summary = cat(1, PlaceCell_session1_ER, PlaceCell_session2_ER);
Mean_PlaceCellER = mean(PlaceCellER_Summary, 'omitnan');
NonPlaceCellER_Summary = cat(1, Non_PlaceCell_session1_ER, Non_PlaceCell_session2_ER);
Mean_NonPlaceCellER = mean(NonPlaceCellER_Summary, 'omitnan');

% Optional: temporary visualization
% % figure
% % [PlaceCell_F, PlaceCell_X, ~, ~] = ecdf(PlaceCellER_Summary,'Function','cdf','Alpha',0.05);
% % [NonPlaceCell_F, NonPlaceCell_X, ~, ~] = ecdf(NonPlaceCellER_Summary,'Function','cdf','Alpha',0.05);
% % plot(PlaceCell_X, PlaceCell_F, 'r', 'LineWidth', 3)
% % hold on
% % plot(NonPlaceCell_X, NonPlaceCell_F, 'b', 'LineWidth', 3)
% % xlabel('Event Rate (spikes/sec)')
% % ylabel('Cumulative Fraction of Cells')
% % xlim([0,0.5])
% % legend('Place Cells', 'Non Place Cells')
% % set(gca,'FontSize',15)
% % hold off

%% Step 3: Process Noldus Behavioral Tracking Data

% Generate a filter: smoothing the speed of the mouse by convolution with a triangular filter 
smooth_win = triang(round(ds.Noldus_Sampling_Frequency * ap.Behavior_Filter_Size));
smooth_win = smooth_win / sum(smooth_win);
% Calculate the displacement of the mouse between each two adjust
% timestamps for session1
Displacement1=zeros(size(Behavior_Tracking1,1),1);
for n=1:size(Behavior_Tracking1,1)
    if n+1<size(Behavior_Tracking1,1)
        Displacement1(n,1)=sqrt((Behavior_Tracking1(n+1,2)-Behavior_Tracking1(n,2))^2+(Behavior_Tracking1(n+1,3)-Behavior_Tracking1(n,3))^2); % calculate displacement at each 0.04sec interval
    end
end
% Fix potential errors: sometimes displacement calculated is nan because of very occassinoal
% Noldus tracking errors
Displacement1(isnan(Displacement1))=0;
Displacement1(isinf(Displacement1))=0;
% calculate the speed during each 0.04 time interval
Speed1=Displacement1./(1/ds.Noldus_Sampling_Frequency); 
% Smooth the speed of animal during session1
Speed1 = conv(Speed1, smooth_win, 'same');


% Perform the same computation for session2
Displacement2=zeros(size(Behavior_Tracking2,1),1);
for n=1:size(Behavior_Tracking2,1)
    if n+1<size(Behavior_Tracking2,1)
        Displacement2(n,1)=sqrt((Behavior_Tracking2(n+1,2)-Behavior_Tracking2(n,2))^2+(Behavior_Tracking2(n+1,3)-Behavior_Tracking2(n,3))^2); 
    end
end
Displacement2(isnan(Displacement2))=0;
Displacement2(isinf(Displacement2))=0;
Speed2=Displacement2./(1/ds.Noldus_Sampling_Frequency);
Speed2 = conv(Speed2, smooth_win, 'same');

% Compute the fraction of time that mouse actively explores for session1&2
Mobility_Session1=(sum(Speed1>ap.Low_Mobility_Threshold))/(size(Speed1,1)); 
Mobility_Session2=(sum(Speed2>ap.Low_Mobility_Threshold))/(size(Speed2,1)); 

% Finishing up this section by concatenating speed info of the animal with
% the original behavioral tracking matrices
Behavior_Tracking1=cat(2, Behavior_Tracking1, Speed1); 
Behavior_Tracking2=cat(2, Behavior_Tracking2, Speed2); 

%% Step 4: Explore the relationship between ER and mobility

% Find position of mouse for each calcium event during session1
% Notes about the organization of Cell_Events1&2 matrices after this step:
% Column 1 = timestamps of calcium event
% Column 2 = cell identity
% Column 3 = event amplitude
% Column 4&5 = X&Y position of the mouse at the event
% Column 6 = speed of the mouse at the event
[Cell_Events1]=Find_XY_Position_Single_Session(Cell_Events1, Behavior_Tracking1, ds.Noldus_Sampling_Frequency);
[Cell_Events2]=Find_XY_Position_Single_Session(Cell_Events2, Behavior_Tracking2, ds.Noldus_Sampling_Frequency);

% find the distribution of speeds during calcium events
Session1_Speed_During_Event = Cell_Events1(:,6);
Session2_Speed_During_Event = Cell_Events2(:,6);

% Optional: temporary visualization
% % [Speed1_F, Speed1_X, ~, ~] = ecdf(Session1_Speed_During_Event,'Function','cdf','Alpha',0.05);
% % [Speed2_F, Speed2_X, ~, ~] = ecdf(Session2_Speed_During_Event,'Function','cdf','Alpha',0.05);
% % figure
% % plot(Speed1_X, Speed1_F, 'r', 'LineWidth', 3)
% % hold on
% % plot(Speed2_X, Speed2_F, 'b', 'LineWidth', 3)
% % xlabel('Animal Speed (cm/sec)')
% % ylabel('Cumulative Fraction of Events')
% % ylim([0,1])
% % xlim([0,30])
% % legend('Session1', 'Session2')
% % set(gca,'FontSize',15)
% % hold off

% calculate the mean ER of each cell during mobility and imobility
% also calculate the overall mean ER during mobility and imobility

% compute for session1 first
session1_cell_ID = unique(Cell_Events1(:,2));
Session1_Cell_ER_During_Imobility = nan(max(session1_cell_ID),1);
Session1_Cell_ER_During_Mobility = nan(max(session1_cell_ID),1);
for n = 1:max(session1_cell_ID)
    Session1_Cell_ER_During_Imobility(n,1) = (sum(Cell_Events1(:,2) == n & Cell_Events1(:,6) < ap.Low_Mobility_Threshold))/((1-Mobility_Session1)*ds.Total_T_Recording);
    Session1_Cell_ER_During_Mobility(n,1) = (sum(Cell_Events1(:,2) == n & Cell_Events1(:,6) > ap.Low_Mobility_Threshold))/(Mobility_Session1*ds.Total_T_Recording);
end
for n = 1:size(Session1_Cell_ER_During_Mobility, 1)
    if Session1_Cell_ER_During_Mobility(n,1) == 0 && Session1_Cell_ER_During_Imobility(n,1) == 0
        Session1_Cell_ER_During_Mobility(n,1) = nan;
        Session1_Cell_ER_During_Imobility(n,1) = nan;
    end
end
Session1_Mean_ER_During_Mobility = mean(Session1_Cell_ER_During_Mobility, 'omitnan');
Session1_Mean_ER_During_Imobility = mean(Session1_Cell_ER_During_Imobility, 'omitnan');


% compute for session2
Session2_cell_ID = unique(Cell_Events2(:,2));
Session2_Cell_ER_During_Imobility = nan(max(Session2_cell_ID),1);
Session2_Cell_ER_During_Mobility = nan(max(Session2_cell_ID),1);
for n = 1:max(Session2_cell_ID)
    Session2_Cell_ER_During_Imobility(n,1) = (sum(Cell_Events2(:,2) == n & Cell_Events2(:,6) < ap.Low_Mobility_Threshold))/((1-Mobility_Session2)*ds.Total_T_Recording);
    Session2_Cell_ER_During_Mobility(n,1) = (sum(Cell_Events2(:,2) == n & Cell_Events2(:,6) > ap.Low_Mobility_Threshold))/(Mobility_Session2*ds.Total_T_Recording);
end
for n = 1:size(Session2_Cell_ER_During_Mobility, 1)
    if Session2_Cell_ER_During_Mobility(n,1) == 0 && Session2_Cell_ER_During_Imobility(n,1) == 0
        Session2_Cell_ER_During_Mobility(n,1) = nan;
        Session2_Cell_ER_During_Imobility(n,1) = nan;
    end
end
Session2_Mean_ER_During_Mobility = mean(Session2_Cell_ER_During_Mobility, 'omitnan');
Session2_Mean_ER_During_Imobility = mean(Session2_Cell_ER_During_Imobility, 'omitnan');

%% Put outputs into a Data_Summary structure

Data_Summary = struct(...
    'animalName', [],...
    'genotype', [],...
    'dateToFind', [],...    
    'drug', [],...
    'session1_dose', [],...
    'session2_dose', [],...
    'session1_context', [],...
    'session2_context', [],...
    'exptParadigm', [],...
    'PlaceCellER_Summary', PlaceCellER_Summary,...
    'NonPlaceCellER_Summary', NonPlaceCellER_Summary,...
    'Mean_PlaceCellER', Mean_PlaceCellER,...
    'Mean_NonPlaceCellER', Mean_NonPlaceCellER,...
    'Session1_Speed_During_Event', Session1_Speed_During_Event,...
    'Session2_Speed_During_Event', Session2_Speed_During_Event,...
    'Session1_Cell_ER_During_Mobility', Session1_Cell_ER_During_Mobility,...
    'Session1_Cell_ER_During_Imobility', Session1_Cell_ER_During_Imobility,...
    'Session2_Cell_ER_During_Mobility', Session2_Cell_ER_During_Mobility,...
    'Session2_Cell_ER_During_Imobility', Session2_Cell_ER_During_Imobility,...
    'Session1_Mean_ER_During_Mobility', Session1_Mean_ER_During_Mobility,...
    'Session1_Mean_ER_During_Imobility', Session1_Mean_ER_During_Imobility,...
    'Session2_Mean_ER_During_Mobility', Session2_Mean_ER_During_Mobility,...
    'Session2_Mean_ER_During_Imobility', Session2_Mean_ER_During_Imobility);

end