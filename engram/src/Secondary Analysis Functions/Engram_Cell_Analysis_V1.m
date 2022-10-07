function Data_Summary = Engram_Cell_Analysis_V1(Inscopix1, Inscopix2, Noldus1, Noldus2, Start_T1, Start_T2, PlaceCell_Summary1, PlaceCell_Summary2, ds, ap)
% IMPORTANT NOTES:
% Since primary analysis is automated, all detected cellsets from each
% animal will be comprehensively registered, events detected, and events
% detected from each recording will be exported as a csv file that has
% global cell name.

% Also, place cell characteristics (mutual information shuffle tests) will be performed
% in a separate function, which has the output Placeness_MetaData, 
% stored as a structure that is located within metadata
% directory. They will be used as a reference for determining whether a
% cell is a place cell or not.

% There are several major changes to this function compared to before

% Input of this function:
% 1. Inscopix1 & Inscopix2 are the pair of event csv files
% 2. Noldus1 & Noldus2 are the pair of behavioral tracking excel files
% 3. Start_T1 & Start_T2 are the starting times of the recording pair (posix time)
% 4. PlaceCell_Summary1 & PlaceCell_Summary2 are two numeric arrays with 0, 1, and NaN that tells you which cell is place cell 
% 5. ds & ap are from function getdsap

% Prerequisites:
% To implement both primary and secondary batch analysis, as well as graphing, you will need the following MATLAB add-on's:
% (P.S.: You will also need IDPS 1.6 or higher installed and add path to Data Processing folder within Inscopix folder)
% Image processing
% Optimization
% Signal Processing
% Statistics and Machine Learning
% Curve Fitting
% Parallel Computing
% Financial Toolbox
% Mapping Toolbox

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

% Extract the total number of cells active in each recording
Num_Cell_Active_Session1 = size(unique(Cell_Events1(:,2)),1);
Num_Cell_Active_Session2 = size(unique(Cell_Events2(:,2)),1);

% Extract the cell identities of cells active during both recordings &
% number of cells active during both recordings
Coactive_Cell_Identity = intersect(unique(Cell_Events1(:,2)), unique(Cell_Events2(:,2)));
Num_Cell_Active_During_Both_Session = size(Coactive_Cell_Identity,1);

% Extract Total Number of Cells Ever Detected During Both Recordings
Total_Cell_Active = Num_Cell_Active_Session1 + Num_Cell_Active_Session2 - Num_Cell_Active_During_Both_Session;

% Extract some parameters
Coactive_Cell_Proportion = Num_Cell_Active_During_Both_Session/Total_Cell_Active;
Session1_ONLY_Cell_Proportion = (Num_Cell_Active_Session1-Num_Cell_Active_During_Both_Session)/Total_Cell_Active;
Session2_ONLY_Cell_Proportion = (Num_Cell_Active_Session2-Num_Cell_Active_During_Both_Session)/Total_Cell_Active;

% Compute probability that an active cell in session1 remains active during session2
Cell_Recurring_Probability = Num_Cell_Active_During_Both_Session/Num_Cell_Active_Session1;

%% Step 2: Extract Calcium Event Rate Parameters

% compute the firing rates of cells detected in each recording
Session1_Event_Rates=nan(max(unique(Cell_Events1(:,2))),1);
for n=(unique(Cell_Events1(:,2)))'
    Session1_Event_Rates(n,1)=sum(Cell_Events1(:,2)==n)/ds.Total_T_Recording; 
end

Session2_Event_Rates=nan(max(unique(Cell_Events2(:,2))),1);
for n=(unique(Cell_Events2(:,2)))'
    Session2_Event_Rates(n,1)=sum(Cell_Events2(:,2)==n)/ds.Total_T_Recording; 
end

% compute average firing rate for each recording session
Session1_Mean_Ca_Event_Rate=mean(Session1_Event_Rates,'omitnan'); 
Session2_Mean_Ca_Event_Rate=mean(Session2_Event_Rates,'omitnan');

% Now calculate the proportion of cells within each pre-defined firing rate range for later plotting & summarizing the data
BinCounts_Session1_Event_Rate = histcounts(Session1_Event_Rates(:),(ap.Ca_Rate_Lower_Bound:ap.Ca_Rate_Interval:ap.Ca_Rate_Upper_Bound));
BinCounts_Session1_Event_Rate_Proportion = BinCounts_Session1_Event_Rate/sum(BinCounts_Session1_Event_Rate); 

BinCounts_Session2_Event_Rate = histcounts(Session2_Event_Rates(:),(ap.Ca_Rate_Lower_Bound:ap.Ca_Rate_Interval:ap.Ca_Rate_Upper_Bound));
BinCounts_Session2_Event_Rate_Proportion = BinCounts_Session2_Event_Rate/sum(BinCounts_Session2_Event_Rate); 

% % Optional: for temporary visualization
% subplot(1,2,1);
% bar((ap.Ca_Rate_Lower_Bound+0.5*ap.Ca_Rate_Interval:ap.Ca_Rate_Interval:ap.Ca_Rate_Upper_Bound),BinCounts_Session1_Event_Rate_Proportion,1,'b','FaceAlpha',0.6);
% ylim([0 0.3]);
% xlabel('calcium event rate', 'FontSize', 11);
% ylabel('proportion of cells', 'FontSize', 11);
% title('AM');
% 
% subplot(1,2,2);
% bar((ap.Ca_Rate_Lower_Bound+0.5*ap.Ca_Rate_Interval:ap.Ca_Rate_Interval:ap.Ca_Rate_Upper_Bound),BinCounts_Session2_Event_Rate_Proportion,1,'b','FaceAlpha',0.6);
% ylim([0 0.3]);
% xlabel('calcium event rate', 'FontSize', 11);
% ylabel('proportion of cells', 'FontSize', 11);
% title('PM');

% Additional analysis: Event Rate Correlation (ER) Analysis
% Prepare the data
a = size(Session1_Event_Rates,1);
b = size(Session2_Event_Rates,1);

if a - b == 1
    Session2_Event_Rates = cat(1, Session2_Event_Rates, 0);
elseif b - a == 1
    Session1_Event_Rates = cat(1, Session1_Event_Rates, 0);
elseif abs(a-b) ~= 1 && a>b
    Session2_Event_Rates = cat(1, Session2_Event_Rates, zeros(a-b, 1));
elseif abs(a-b) ~= 1 && b>a
    Session1_Event_Rates = cat(1, Session1_Event_Rates, zeros(b-a, 1));
else
    disp('both event rate vectors have same size, good to go')
end

for n = 1:size(Session1_Event_Rates,1)
    if isnan(Session1_Event_Rates(n,1)) || isnan(Session2_Event_Rates(n,1))
        Session1_Event_Rates(n,1) = nan;
        Session2_Event_Rates(n,1) = nan;
    elseif Session1_Event_Rates(n,1) == 0 || Session2_Event_Rates(n,1) == 0
        Session1_Event_Rates(n,1) = nan;
        Session2_Event_Rates(n,1) = nan;
    end
end
Session1_Event_Rates(isnan(Session1_Event_Rates)) = [];
Session2_Event_Rates(isnan(Session2_Event_Rates)) = [];

% calculate ER correlation
Event_Rate_Correlation = corr(Session1_Event_Rates(:), Session2_Event_Rates(:));

%% Step 2.1 Compute Calcium Event Amplitude (after linear adjustment) Distribution for Each Recording

% get adjusted calcium event amplitude distribution for each session
Session1_Event_Amps = Cell_Events1(:,3);
Session2_Event_Amps = Cell_Events2(:,3);

% get adjusted calcium event amplitude mean for each session
Session1_Event_Amps_Mean = mean(Session1_Event_Amps,'all','omitnan');
Session2_Event_Amps_Mean = mean(Session2_Event_Amps,'all','omitnan'); 

% % Optional: for temporary visualization
% subplot(1,2,1);
% bar((ap.Ca_Amp_Lower_Bound+0.5*ap.Ca_Amp_Interval:ap.Ca_Amp_Interval:ap.Ca_Amp_Upper_Bound),BinCounts_Session1_Event_Amp_Proportion,1,'b','FaceAlpha',0.6);
% ylim([0 0.5]);
% xlabel('calcium event amplitude', 'FontSize', 11);
% ylabel('proportion of cells', 'FontSize', 11);
% title('AM');
% 
% subplot(1,2,2);
% bar((ap.Ca_Amp_Lower_Bound+0.5*ap.Ca_Amp_Interval:ap.Ca_Amp_Interval:ap.Ca_Amp_Upper_Bound),BinCounts_Session2_Event_Amp_Proportion,1,'b','FaceAlpha',0.6);
% ylim([0 0.5]);
% xlabel('calcium event amplitude', 'FontSize', 11);
% ylabel('proportion of cells', 'FontSize', 11);
% title('PM');


%% Optional Additional Step: Choose a Portion Inscopix & Noldus Data to Include in the Following Analysis

[Behavior_Tracking1, Cell_Events1] = Segment_NoldusAndInscopix(Behavior_Tracking1, Cell_Events1, ap.Engram_Start_T, ap.Engram_End_T);
[Behavior_Tracking2, Cell_Events2] = Segment_NoldusAndInscopix(Behavior_Tracking2, Cell_Events2, ap.Engram_Start_T, ap.Engram_End_T);

%% Step 3 (Optional): Linearly adjust each Ca2+ event according to its peak amplitude, measured in dF/noise (refer to IDPS documentation for further details)

% Adding this step will decrease processing speed by quite a lot
% Linear amplitude correction for Cell_Events1

% perform some reorganizations
for n=1:size(Cell_Events1,1)
    if Cell_Events1(n,3)<=1
        Cell_Events1(n,3)=1; 
    elseif Cell_Events1(n,3)>1 && rem(Cell_Events1(n,3),1)<0.5
       Cell_Events1(n,3)=Cell_Events1(n,3)-rem(Cell_Events1(n,3),1); 
    else
       Cell_Events1(n,3)=Cell_Events1(n,3)+(1-rem(Cell_Events1(n,3),1)); 
    end
end

R1=reshape(Cell_Events1(:,1),[1,1,size(Cell_Events1,1)]);
R2=reshape(Cell_Events1(:,2),[1,1,size(Cell_Events1,1)]);
R3=reshape(Cell_Events1(:,3),[1,1,size(Cell_Events1,1)]);

CatR=cat(2,R1,R2,R3);

CorrectedCat=nan(max(Cell_Events1(:,3),[],'all'), size(CatR,2), size(CatR,3));
for n=1:size(CatR,3)
    CorrectedCat(1:CatR(1,3,n),:,n)=repmat(CatR(1,1:3,n),[CatR(1,3,n),1,1]); 
end

% Add this back
for n=2:size(CorrectedCat,1)
    CorrectedCat(n,1,:)=CorrectedCat(n,1,:)+(1/ds.Noldus_Sampling_Frequency)*(n-1); 
end

C1=reshape(CorrectedCat(:,1,:),[size(CorrectedCat,1)*size(CorrectedCat,3),1]);
C2=reshape(CorrectedCat(:,2,:),[size(CorrectedCat,1)*size(CorrectedCat,3),1]);
C3=reshape(CorrectedCat(:,3,:),[size(CorrectedCat,1)*size(CorrectedCat,3),1]);

C1(isnan(C1))=[];
C2(isnan(C2))=[];
C3(isnan(C3))=[];

Cell_Events1=cat(2,C1,C2,C3); 

% fix bugs of going beyond total recording time:
for n=1:size(Cell_Events1,1)
    if Cell_Events1(n,1)>ds.Total_T_Recording
        Cell_Events1(n,:)=nan;
    end
end
Cell_Events1(isnan(Cell_Events1))=[];
Cell_Events1=reshape(Cell_Events1,[],3);

% Same computation for session 2
for n=1:size(Cell_Events2,1)
    if Cell_Events2(n,3)<=1
        Cell_Events2(n,3)=1; 
    elseif Cell_Events2(n,3)>1 && rem(Cell_Events2(n,3),1)<0.5
       Cell_Events2(n,3)=Cell_Events2(n,3)-rem(Cell_Events2(n,3),1); 
    else
       Cell_Events2(n,3)=Cell_Events2(n,3)+(1-rem(Cell_Events2(n,3),1)); 
    end
end


R1=reshape(Cell_Events2(:,1),[1,1,size(Cell_Events2,1)]);
R2=reshape(Cell_Events2(:,2),[1,1,size(Cell_Events2,1)]);
R3=reshape(Cell_Events2(:,3),[1,1,size(Cell_Events2,1)]);

CatR=cat(2,R1,R2,R3);

CorrectedCat=nan(max(Cell_Events2(:,3),[],'all'), size(CatR,2), size(CatR,3));
for n=1:size(CatR,3)
    CorrectedCat(1:CatR(1,3,n),:,n)=repmat(CatR(1,1:3,n),[CatR(1,3,n),1,1]); 
end

% Add this back
for n=2:size(CorrectedCat,1)
    CorrectedCat(n,1,:)=CorrectedCat(n,1,:)+(1/ds.Noldus_Sampling_Frequency)*(n-1); 
end

C1=reshape(CorrectedCat(:,1,:),[size(CorrectedCat,1)*size(CorrectedCat,3),1]);
C2=reshape(CorrectedCat(:,2,:),[size(CorrectedCat,1)*size(CorrectedCat,3),1]);
C3=reshape(CorrectedCat(:,3,:),[size(CorrectedCat,1)*size(CorrectedCat,3),1]);

C1(isnan(C1))=[];
C2(isnan(C2))=[];
C3(isnan(C3))=[];

Cell_Events2=cat(2,C1,C2,C3);

% fix bugs
for n=1:size(Cell_Events2,1)
    if Cell_Events2(n,1)>ds.Total_T_Recording
        Cell_Events2(n,:)=nan;
    end
end
Cell_Events2(isnan(Cell_Events2))=[];
Cell_Events2=reshape(Cell_Events2,[],3);

%% Additional Optional Step: Exclude Events with Amplitude under a Certain Threshold

for n = 1:size(Cell_Events1,1)
    if Cell_Events1(n,3)<=ap.Amp_Exclusion_Threshold
        Cell_Events1(n,:)=nan;
    end
end
Cell_Events1(isnan(Cell_Events1)) = [];
Cell_Events1 = reshape(Cell_Events1, [], 3);

for n = 1:size(Cell_Events2,1)
    if Cell_Events2(n,3)<=ap.Amp_Exclusion_Threshold
        Cell_Events2(n,:)=nan;
    end
end
Cell_Events2(isnan(Cell_Events2)) = [];
Cell_Events2 = reshape(Cell_Events2, [], 3);

%% Step 4: Process Noldus Behavioral Tracking Data

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

% give an indication of whether the recording should be included in
% statistical analysis (mobility < certain value should not be included
% within the summary of place cell & engram statistics)
if Mobility_Session1 >= 0.15
    Mobility_Session1_Pass = 1;
else
    Mobility_Session1_Pass = 0;
end

if Mobility_Session2 >= 0.15
    Mobility_Session2_Pass = 1;
else
    Mobility_Session2_Pass = 0;
end

% Finishing up this section by concatenating speed info of the animal with
% the original behavioral tracking matrices
Behavior_Tracking1=cat(2, Behavior_Tracking1, Speed1); 
Behavior_Tracking2=cat(2, Behavior_Tracking2, Speed2); 

% Optional graph: plot mouse speed over time
% % plot(Behavior_Tracking1(:,1),Behavior_Tracking1(:,4),'-b','lineWidth',1.5)
% % xlim([0,ds.Total_T_Recording]);
% % ylim([0,60]);
% % title('Session 1&2: mouse velocity vs. time','FontSize',14)
% % xlabel('time (sec)','FontSize',14)
% % ylabel('speed (cm/sec)','FontSize',14)
% % hold on;
% % plot(Behavior_Tracking2(:,1),Behavior_Tracking2(:,4),'-g','lineWidth',1.5)
% % xlim([0,ds.Total_T_Recording]);
% % ylim([0,60]);
% % legend('session 1','session 2')
% % hold off;

%% Step 5: Generate Occupancy_Map for Both Sessions from Noldus Behavioral Tracking Data

% Calculate & extend occupancy map boundaries
A=min(Behavior_Tracking1(:,2),[],'omitnan');
B=max(Behavior_Tracking1(:,2),[],'omitnan');
C=min(Behavior_Tracking1(:,3),[],'omitnan');
D=max(Behavior_Tracking1(:,3),[],'omitnan');

% extend the arena size to make sure all behavioral tracking points fall
% within occupancy matrix
Map1_LowX_Bound=fix(A)-ap.Arena_Add_Edge; 
Map1_HighX_Bound=fix(B)+ap.Arena_Add_Edge;
Map1_LowY_Bound=fix(C)-ap.Arena_Add_Edge;
Map1_HighY_Bound=fix(D)+ap.Arena_Add_Edge;


% Calculate & extend occupancy map boundaries
D=min(Behavior_Tracking1(:,2),[],'omitnan');
E=max(Behavior_Tracking1(:,2),[],'omitnan');
F=min(Behavior_Tracking1(:,3),[],'omitnan');
G=max(Behavior_Tracking1(:,3),[],'omitnan');

% extend the arena size to make sure all behavioral tracking points fall
% within occupancy matrix
Map2_LowX_Bound=fix(D)-ap.Arena_Add_Edge; 
Map2_HighX_Bound=fix(E)+ap.Arena_Add_Edge;
Map2_LowY_Bound=fix(F)-ap.Arena_Add_Edge;
Map2_HighY_Bound=fix(G)+ap.Arena_Add_Edge;


[Corrected_Occupancy_Map1] = Compute_Occupancy_Map(Behavior_Tracking1, ds, ap);
[Corrected_Occupancy_Map2] = Compute_Occupancy_Map(Behavior_Tracking2, ds, ap);

%% Step 6: Compute the XY-positions of the Mouse Given a Calcium Event Timestamp & Exclude Events During Immobility

% Find position of mouse for each calcium event during session1
% Notes about the organization of Cell_Events1&2 matrices after this step:
% Column 1 = timestamps of calcium event
% Column 2 = cell identity
% Column 3 = event amplitude
% Column 4&5 = X&Y position of the mouse at the event
% Column 6 = speed of the mouse at the event
[Cell_Events1]=Find_XY_Position_Single_Session(Cell_Events1, Behavior_Tracking1, ds.Noldus_Sampling_Frequency);
% Exclude events during immobility
for n=1:size(Cell_Events1,1)
    if Cell_Events1(n,6)>ap.High_Mobility_Threshold || Cell_Events1(n,6)<ap.Low_Mobility_Threshold
        Cell_Events1(n,1:6)=nan(1,6);
    end
end
% Clean up
Cell_Events1(any(isnan(Cell_Events1),2),:) = []; 


% Same computation for session2
[Cell_Events2]=Find_XY_Position_Single_Session(Cell_Events2, Behavior_Tracking2, ds.Noldus_Sampling_Frequency);
for n=1:size(Cell_Events2,1)
    if Cell_Events2(n,6)>ap.High_Mobility_Threshold || Cell_Events2(n,6)<ap.Low_Mobility_Threshold
        Cell_Events2(n,1:6)=nan(1,6); 
    end
end
Cell_Events2(any(isnan(Cell_Events2),2),:) = []; 

%% Step 7: Extract Cells (and Their Events) That Are Active During Both Sessions From Each Cell_Events Matrix

% Extract such info for session1
Cell_Active_Both_Session_Events1 = nan(size(Cell_Events1,1), size(Cell_Events1,2));
for ii = 1:size(Cell_Events1,1)
    if ismember(Cell_Events1(ii,2), Coactive_Cell_Identity(:))
        Cell_Active_Both_Session_Events1 (ii,:) = Cell_Events1 (ii,:);
    end
end
% Clean up
Cell_Active_Both_Session_Events1(isnan(Cell_Active_Both_Session_Events1)) = [];
Cell_Active_Both_Session_Events1 = reshape(Cell_Active_Both_Session_Events1, [], size(Cell_Events1,2));


% Perform the same computation for session2
Cell_Active_Both_Session_Events2 = nan(size(Cell_Events2,1), size(Cell_Events2,2));
for ii = 1:size(Cell_Events2,1)
    if ismember(Cell_Events2(ii,2), Coactive_Cell_Identity(:))
        Cell_Active_Both_Session_Events2 (ii,:) = Cell_Events2 (ii,:);
    end
end
Cell_Active_Both_Session_Events2(isnan(Cell_Active_Both_Session_Events2)) = [];
Cell_Active_Both_Session_Events2 = reshape(Cell_Active_Both_Session_Events2, [], size(Cell_Events2,2));

%% Step 8: Compute Ca_Event_Map For Both Sessions

% Calculate the maximum number of calcium events observed from one cell that appeared during session1
% e.g. If Cell #20 fires 300 times during session1 recording, and this value is more than any other cell
% This computation is for organizing data in the next steps
n=(unique(Cell_Active_Both_Session_Events1(:,2)))';
Session1_Cell_Fire_Most_Num_Event=max(sum(Cell_Active_Both_Session_Events1(:,2)==n)); 

% Reorganize data for generating event maps for each cell
Cell_Active_Both_Session_Reorganize1=nan(Session1_Cell_Fire_Most_Num_Event, size(Cell_Active_Both_Session_Events1,2), size(unique(Cell_Active_Both_Session_Events1(:,2)),1));
for n=(unique(Cell_Active_Both_Session_Events1(:,2)))'
    y=sum(Cell_Active_Both_Session_Events1(:,2)==n);
    Cell_Active_Both_Session_Reorganize1(:,:,n)=reshape([Cell_Active_Both_Session_Events1(cat(1,Cell_Active_Both_Session_Events1(:,2))==n,:);NaN(Session1_Cell_Fire_Most_Num_Event-y,size(Cell_Active_Both_Session_Events1,2))],Session1_Cell_Fire_Most_Num_Event,size(Cell_Active_Both_Session_Events1,2)); 
end
Cell_Active_Both_Session_Reorganize1(Cell_Active_Both_Session_Reorganize1==0)=nan;

% Compute Ca Event Map, which is how many calcium spikes happen in a
% particular spatial bin
Cell_Active_Both_Session_Ca_Event_Map1=nan(size(Corrected_Occupancy_Map1,1),size(Corrected_Occupancy_Map1,2),size(Cell_Active_Both_Session_Reorganize1,3));
for k=1:size(Cell_Active_Both_Session_Reorganize1,3)
    Cell_Active_Both_Session_Ca_Event_Map1(:,:,k)=hist3(Cell_Active_Both_Session_Reorganize1(:,4:5,k),'ctrs',{Map1_LowX_Bound:(Map1_HighX_Bound-Map1_LowX_Bound)/ap.Map_Division_For_PV:Map1_HighX_Bound  Map1_LowY_Bound:(Map1_HighY_Bound-Map1_LowY_Bound)/ap.Map_Division_For_PV:Map1_HighY_Bound});
end
% Rotate map so that it matches the actual experimental setup
Cell_Active_Both_Session_Ca_Event_Map1=rot90(Cell_Active_Both_Session_Ca_Event_Map1,1); 

% Fix bugs: make sure non-existent cells are represented by a matrix of zeros
for n=1:size(Cell_Active_Both_Session_Ca_Event_Map1,3)
    if ismember(n, Coactive_Cell_Identity)
        Cell_Active_Both_Session_Ca_Event_Map1(:,:,n)=Cell_Active_Both_Session_Ca_Event_Map1(:,:,n);
    else 
        Cell_Active_Both_Session_Ca_Event_Map1(:,:,n)=zeros(size(Cell_Active_Both_Session_Ca_Event_Map1,1), size(Cell_Active_Both_Session_Ca_Event_Map1,2));
    end
end




% Perform exactly the same computation for session2
n=(unique(Cell_Active_Both_Session_Events2(:,2)))';
Session2_Cell_Fire_Most_Num_Event=max(sum(Cell_Active_Both_Session_Events2(:,2)==n)); 

Cell_Active_Both_Session_Reorganize2=nan(Session2_Cell_Fire_Most_Num_Event, size(Cell_Active_Both_Session_Events2,2), size(unique(Cell_Active_Both_Session_Events2(:,2)),1));
for n=(unique(Cell_Active_Both_Session_Events2(:,2)))'
    y=sum(Cell_Active_Both_Session_Events2(:,2)==n);
    Cell_Active_Both_Session_Reorganize2(:,:,n)=reshape([Cell_Active_Both_Session_Events2(cat(1,Cell_Active_Both_Session_Events2(:,2))==n,:);NaN(Session2_Cell_Fire_Most_Num_Event-y,size(Cell_Active_Both_Session_Events2,2))],Session2_Cell_Fire_Most_Num_Event,size(Cell_Active_Both_Session_Events2,2));
end

Cell_Active_Both_Session_Ca_Event_Map2=nan(size(Corrected_Occupancy_Map2,1),size(Corrected_Occupancy_Map2,2),size(Cell_Active_Both_Session_Reorganize2,3));
for k=1:size(Cell_Active_Both_Session_Reorganize2,3)
    Cell_Active_Both_Session_Ca_Event_Map2(:,:,k)=hist3(Cell_Active_Both_Session_Reorganize2(:,4:5,k),'ctrs',{Map2_LowX_Bound:(Map2_HighX_Bound-Map2_LowX_Bound)/ap.Map_Division_For_PV:Map2_HighX_Bound  Map2_LowY_Bound:(Map2_HighY_Bound-Map2_LowY_Bound)/ap.Map_Division_For_PV:Map2_HighY_Bound}); 
end
Cell_Active_Both_Session_Ca_Event_Map2=rot90(Cell_Active_Both_Session_Ca_Event_Map2,1);

for n=1:size(Cell_Active_Both_Session_Ca_Event_Map2,3)
    if ismember(n, Coactive_Cell_Identity)
        Cell_Active_Both_Session_Ca_Event_Map2(:,:,n)=Cell_Active_Both_Session_Ca_Event_Map2(:,:,n);
    else 
        Cell_Active_Both_Session_Ca_Event_Map2(:,:,n)=zeros(size(Cell_Active_Both_Session_Ca_Event_Map2,1), size(Cell_Active_Both_Session_Ca_Event_Map2,2));
    end
end


%% Step 9: Compute Ca_Event_Rate_Map For Both Sessions

% Compute Ca_Event_Rate_Map for session1
Cell_Active_Both_Session_Ca_Event_Rate_Map1=Cell_Active_Both_Session_Ca_Event_Map1./Corrected_Occupancy_Map1; 
Cell_Active_Both_Session_Ca_Event_Rate_Map1(isnan(Cell_Active_Both_Session_Ca_Event_Rate_Map1))=0; 
Cell_Active_Both_Session_Ca_Event_Rate_Map1=reshape(Cell_Active_Both_Session_Ca_Event_Rate_Map1,size(Corrected_Occupancy_Map1,1),size(Corrected_Occupancy_Map1,2),[]); 

% Compute Ca_Event_Rate_Map for session2
Cell_Active_Both_Session_Ca_Event_Rate_Map2=Cell_Active_Both_Session_Ca_Event_Map2./Corrected_Occupancy_Map2;
Cell_Active_Both_Session_Ca_Event_Rate_Map2(isnan(Cell_Active_Both_Session_Ca_Event_Rate_Map2))=0;
Cell_Active_Both_Session_Ca_Event_Rate_Map2=reshape(Cell_Active_Both_Session_Ca_Event_Rate_Map2,size(Corrected_Occupancy_Map2,1),size(Corrected_Occupancy_Map2,2),[]);

%Fix bugs: non-existent cells are represented by zero matrices 
for n=1:size(Cell_Active_Both_Session_Ca_Event_Rate_Map1,3)
    if ismember(n,unique(Cell_Active_Both_Session_Events1(:,2)))
        Cell_Active_Both_Session_Ca_Event_Rate_Map1(:,:,n)=Cell_Active_Both_Session_Ca_Event_Rate_Map1(:,:,n);
    else
        Cell_Active_Both_Session_Ca_Event_Rate_Map1(:,:,n)=zeros(size(Cell_Active_Both_Session_Ca_Event_Rate_Map1,1),size(Cell_Active_Both_Session_Ca_Event_Rate_Map1,2));
    end
end

for n=1:size(Cell_Active_Both_Session_Ca_Event_Rate_Map2,3)
    if ismember(n,unique(Cell_Active_Both_Session_Events2(:,2)))
        Cell_Active_Both_Session_Ca_Event_Rate_Map2(:,:,n)=Cell_Active_Both_Session_Ca_Event_Rate_Map2(:,:,n);
    else
        Cell_Active_Both_Session_Ca_Event_Rate_Map2(:,:,n)=zeros(size(Cell_Active_Both_Session_Ca_Event_Rate_Map2,1),size(Cell_Active_Both_Session_Ca_Event_Rate_Map2,2));
    end
end

% Fix bugs: make sure that sizes of calcium rate map matrices are the same
d1=size(Cell_Active_Both_Session_Ca_Event_Rate_Map1,3)-size(Cell_Active_Both_Session_Ca_Event_Rate_Map2,3);
d2=size(Cell_Active_Both_Session_Ca_Event_Rate_Map2,3)-size(Cell_Active_Both_Session_Ca_Event_Rate_Map1,3);

if size(Cell_Active_Both_Session_Ca_Event_Rate_Map1,3)<size(Cell_Active_Both_Session_Ca_Event_Rate_Map2,3) && d2>1
    Cell_Active_Both_Session_Ca_Event_Rate_Map1(1:size(Cell_Active_Both_Session_Ca_Event_Rate_Map1,1),1:size(Cell_Active_Both_Session_Ca_Event_Rate_Map1,2),(end+1):size(Cell_Active_Both_Session_Ca_Event_Rate_Map2,3))=zeros(size(Cell_Active_Both_Session_Ca_Event_Rate_Map1,1),size(Cell_Active_Both_Session_Ca_Event_Rate_Map1,2),size(Cell_Active_Both_Session_Ca_Event_Rate_Map2,3)-size(Cell_Active_Both_Session_Ca_Event_Rate_Map1,3)); 
elseif size(Cell_Active_Both_Session_Ca_Event_Rate_Map1,3)<size(Cell_Active_Both_Session_Ca_Event_Rate_Map2,3) && d2==1
    Cell_Active_Both_Session_Ca_Event_Rate_Map1(1:size(Cell_Active_Both_Session_Ca_Event_Rate_Map1,1),1:size(Cell_Active_Both_Session_Ca_Event_Rate_Map1,2),(end+1))=zeros(size(Cell_Active_Both_Session_Ca_Event_Rate_Map1,1),size(Cell_Active_Both_Session_Ca_Event_Rate_Map1,2),size(Cell_Active_Both_Session_Ca_Event_Rate_Map2,3)-size(Cell_Active_Both_Session_Ca_Event_Rate_Map1,3)); 
end

if size(Cell_Active_Both_Session_Ca_Event_Rate_Map1,3)>size(Cell_Active_Both_Session_Ca_Event_Rate_Map2,3) && d1>1
    Cell_Active_Both_Session_Ca_Event_Rate_Map2(1:size(Cell_Active_Both_Session_Ca_Event_Rate_Map2,1),1:size(Cell_Active_Both_Session_Ca_Event_Rate_Map2,2),(end+1):size(Cell_Active_Both_Session_Ca_Event_Rate_Map1,3))=zeros(size(Cell_Active_Both_Session_Ca_Event_Rate_Map2,1),size(Cell_Active_Both_Session_Ca_Event_Rate_Map2,2),size(Cell_Active_Both_Session_Ca_Event_Rate_Map1,3)-size(Cell_Active_Both_Session_Ca_Event_Rate_Map2,3)); 
elseif size(Cell_Active_Both_Session_Ca_Event_Rate_Map1,3)>size(Cell_Active_Both_Session_Ca_Event_Rate_Map2,3) && d1==1
    Cell_Active_Both_Session_Ca_Event_Rate_Map2(1:size(Cell_Active_Both_Session_Ca_Event_Rate_Map2,1),1:size(Cell_Active_Both_Session_Ca_Event_Rate_Map2,2),(end+1))=zeros(size(Cell_Active_Both_Session_Ca_Event_Rate_Map2,1),size(Cell_Active_Both_Session_Ca_Event_Rate_Map2,2),size(Cell_Active_Both_Session_Ca_Event_Rate_Map1,3)-size(Cell_Active_Both_Session_Ca_Event_Rate_Map2,3)); 
end

%% Step 10: Compute Gaussian_Smoothed_Ca_Event_Rate_Maps For Both Sessions

Cell_Active_Both_Session_Smoothed_Ca_Event_Rate_Map1=zeros(size(Cell_Active_Both_Session_Ca_Event_Rate_Map1,1), size(Cell_Active_Both_Session_Ca_Event_Rate_Map1,2), size(Cell_Active_Both_Session_Ca_Event_Rate_Map1,3));
for n=1:size(Cell_Active_Both_Session_Ca_Event_Rate_Map1,3)
    Cell_Active_Both_Session_Smoothed_Ca_Event_Rate_Map1(:,:,n)=filter2(ap.Gaussian_Filter,Cell_Active_Both_Session_Ca_Event_Rate_Map1(:,:,n)); 
end
% fix bugs & clean up
Cell_Active_Both_Session_Smoothed_Ca_Event_Rate_Map1(isinf(Cell_Active_Both_Session_Smoothed_Ca_Event_Rate_Map1))=0;
Cell_Active_Both_Session_Smoothed_Ca_Event_Rate_Map1(isnan(Cell_Active_Both_Session_Smoothed_Ca_Event_Rate_Map1))=0;

Cell_Active_Both_Session_Smoothed_Ca_Event_Rate_Map2=zeros(size(Cell_Active_Both_Session_Ca_Event_Rate_Map2,1), size(Cell_Active_Both_Session_Ca_Event_Rate_Map2,2), size(Cell_Active_Both_Session_Ca_Event_Rate_Map2,3));
for n=1:size(Cell_Active_Both_Session_Ca_Event_Rate_Map2,3)
    Cell_Active_Both_Session_Smoothed_Ca_Event_Rate_Map2(:,:,n)=filter2(ap.Gaussian_Filter,Cell_Active_Both_Session_Ca_Event_Rate_Map2(:,:,n));
end
% fix bugs & clean up
Cell_Active_Both_Session_Smoothed_Ca_Event_Rate_Map2(isinf(Cell_Active_Both_Session_Smoothed_Ca_Event_Rate_Map2))=0;
Cell_Active_Both_Session_Smoothed_Ca_Event_Rate_Map2(isnan(Cell_Active_Both_Session_Smoothed_Ca_Event_Rate_Map2))=0;


%% Step 11: Generate Ca_Event_Rate_Maps That are Used Just For Graphing 

[Smoothed_Ca_Rate_Map1_ForGraph, Smoothed_Ca_Rate_Map2_ForGraph]=Generate_Smoothed_Ca_Event_Maps_For_Graphing(Cell_Active_Both_Session_Events1, Cell_Active_Both_Session_Events2, Behavior_Tracking1, Behavior_Tracking2, ds, ap);

% Optional Graphing to visualize place field of cells, if any
% Note, this cell has to be active during both sessions, or else you cannot graph it correctly
% % n=;
% % subplot(2,2,1);
% % plot(Behavior_Tracking1(:,2),Behavior_Tracking1(:,3));
% % xlim([min(Behavior_Tracking1(:,2))-1,max(Behavior_Tracking1(:,2)+1)]);
% % ylim([min(Behavior_Tracking1(:,3))-1,max(Behavior_Tracking1(:,3)+1)]);
% % hold on
% % scatter(Cell_Active_Both_Session_Reorganize1(:,4,n), Cell_Active_Both_Session_Reorganize1(:,5,n), 70, 'filled'); % plot calcium event dots on AM mouse trace (might have to adjust dot size manually)
% % set(gca,'XColor','none');
% % set(gca,'YColor','none');
% % hold off
% % 
% % subplot(2,2,2);
% % plot(Behavior_Tracking2(:,2),Behavior_Tracking2(:,3));
% % xlim([min(Behavior_Tracking2(:,2))-1,max(Behavior_Tracking2(:,2)+1)]);
% % ylim([min(Behavior_Tracking2(:,3))-1,max(Behavior_Tracking2(:,3)+1)]);
% % hold on
% % scatter(Cell_Active_Both_Session_Reorganize2(:,4,n), Cell_Active_Both_Session_Reorganize2(:,5,n), 70, 'filled'); % plot calcium event dots on PM mouse trace
% % set(gca,'XColor','none');
% % set(gca,'YColor','none');
% % hold off
% % 
% % subplot(2,2,3);
% % heatmap(Smoothed_Ca_Rate_Map1_ForGraph(:,:,n),'CellLabelColor','none'); % plot intensity-corrected Ca rate heat maps for AM session
% % colormap jet; 
% % Ax = gca;
% % Ax.XDisplayLabels = nan(size(Ax.XDisplayData));
% % Ax.YDisplayLabels = nan(size(Ax.YDisplayData));
% % 
% % subplot(2,2,4);
% % heatmap(Smoothed_Ca_Rate_Map2_ForGraph(:,:,n),'CellLabelColor','none'); % plot intensity-corrected Ca rate heat maps for PM session
% % colormap jet;
% % Ax = gca;
% % Ax.XDisplayLabels = nan(size(Ax.XDisplayData));
% % Ax.YDisplayLabels = nan(size(Ax.YDisplayData));


%% Step 12: Calculate the Number of Place Cells Observed During Two Sessions

% Make sure that place cell matrices have the same size
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

% Compute # of Coactive Place Cells in Both Sessions 
% Criteria: it has to appear during both sessions, and it has to at least be a place
% cell in either one of the two sessions
K = zeros(size(PlaceCell_Summary1,1),1);
for n = 1:size(PlaceCell_Summary1,1)
    if ~isnan(PlaceCell_Summary1(n,1)) && ~isnan(PlaceCell_Summary2(n,1))
        K(n,1) = (PlaceCell_Summary1(n,1) == 1 || PlaceCell_Summary2(n,1) == 1);
    end
end
Coactive_Num_Place_Cell = sum(K);

% Compute total # of place cells in both sessions
% Criteria: as long as it is a place cell in either session (regardless of whether if reappears or not)
Q = zeros(size(PlaceCell_Summary1,1),1);
for n = 1:size(PlaceCell_Summary1,1)
        Q(n,1) = (PlaceCell_Summary1(n,1) == 1 || PlaceCell_Summary2(n,1) == 1);
end
Total_Num_Place_Cell = sum(Q);

%% Step 13: Compute Rate Map Correlation & Run Corresponding Shuffle Test 

% Calculate shuffled rate map correlation (null distribution)
Shuffled_Rate_Map_Corr_values=nan(size(Cell_Active_Both_Session_Smoothed_Ca_Event_Rate_Map1,3),1,ap.Ca_Rate_Map_Num_Shuffle);
for n=1:ap.Ca_Rate_Map_Num_Shuffle
    [Shuffled_Rate_Map_Corr_values(:,1,n)]=Shuffle_Smoothed_Rate_Map_Correlation(Coactive_Num_Place_Cell, Cell_Active_Both_Session_Smoothed_Ca_Event_Rate_Map1, Cell_Active_Both_Session_Smoothed_Ca_Event_Rate_Map2, PlaceCell_Summary1, PlaceCell_Summary2, ap);
end
Mean_Shuffled_Rate_Map_Corr = nan(size(Shuffled_Rate_Map_Corr_values,1),1);
for n = 1:size(Shuffled_Rate_Map_Corr_values,1)
    Mean_Shuffled_Rate_Map_Corr(n,1) = mean(Shuffled_Rate_Map_Corr_values(n,1,:));
end
    
Shuffled_RateMap_Corr_Memory_Index = mean(Mean_Shuffled_Rate_Map_Corr, 'omitnan');

% Calculate RM correlation distrbution for the actual data
[Rate_Map_Corr, RateMap_Corr_Memory_Index, RateMap_Corr_Coherent_Rot, Rate_Map_Correlation_Proportion, Shuffled_Rate_Map_Correlation_Proportion]=RM_Correlation_Analysis(Coactive_Num_Place_Cell, Shuffled_Rate_Map_Corr_values, Cell_Active_Both_Session_Smoothed_Ca_Event_Rate_Map1, Cell_Active_Both_Session_Smoothed_Ca_Event_Rate_Map2, PlaceCell_Summary1, PlaceCell_Summary2, ap);

%% Step 15: Compute Population-Vector (PV) Correlation and Perform PV Shuffle Test

% Calculate shuffled PV null distribution
Shuffled_PV_Dist = nan(1, size(Cell_Active_Both_Session_Smoothed_Ca_Event_Rate_Map1,1)*size(Cell_Active_Both_Session_Smoothed_Ca_Event_Rate_Map1,2), ap.PV_Corr_Num_Shuffle);
Shuffled_PV_Mean = nan(ap.PV_Corr_Num_Shuffle, 1);
for n=1:ap.PV_Corr_Num_Shuffle
    [Shuffled_PV_Dist(1,:,n), Shuffled_PV_Mean(n,1)] = Shuffle_PV(Coactive_Num_Place_Cell, Cell_Active_Both_Session_Smoothed_Ca_Event_Rate_Map1, Cell_Active_Both_Session_Smoothed_Ca_Event_Rate_Map2, PlaceCell_Summary1, PlaceCell_Summary2, ap);
end
Shuffled_PV_Corr_Memory_Index = mean(Shuffled_PV_Mean, 'omitnan');

Shuffled_PV_Mean_Dist = nan(1, size(Shuffled_PV_Dist,2));
for n = 1:size(Shuffled_PV_Dist,2)
    Shuffled_PV_Mean_Dist (1,n) = mean(Shuffled_PV_Dist(1,n,:), 'omitnan');
end

% Calculate PV correlation distrbution for the actual data
[PV_Corr_Values, PV_Corr_Memory_Index, PV_Corr_Coherent_Rot, PV_Proportion, Shuffled_PV_Proportion] = PV_Correlation_Analysis(Coactive_Num_Place_Cell, Shuffled_PV_Dist, Cell_Active_Both_Session_Smoothed_Ca_Event_Rate_Map1, Cell_Active_Both_Session_Smoothed_Ca_Event_Rate_Map2, PlaceCell_Summary1, PlaceCell_Summary2, ap);

% Optional Graphing
% % bar((ap.PV_Corr_Lower_Bound+0.5*ap.PV_Corr_Interval:ap.PV_Corr_Interval:ap.PV_Corr_Upper_Bound),PV_Proportion,1,'b','FaceAlpha',0.6);
% % ylim([0 0.4]);
% % xlabel('PV-correlation values', 'FontSize', 11);
% % ylabel('proportion of spatial bins', 'FontSize', 11);
% % hold on;
% % plot((ap.PV_Corr_Lower_Bound+ap.PV_Corr_Interval:ap.PV_Corr_Interval:ap.PV_Corr_Upper_Bound),Shuffled_PV_Proportion,'--k','lineWidth',1.5);
% % hold off;
% % legend('Proportion of spatial bins','Shuffle');

%% Step 16: Perform Memory Development Analysis (Optional)

% The basic concept is the following: we assume that mice start to develop stable
% spatial memory over the first 2-3min when presented with a new context,
% and they are able to retrieve that memory during a second exposure.
% Therefore, memory formation is a time-dependent function that could be
% potentially examined with this current analysis.

% Memory should be formed last 5min of each recording, and it will be our
% reference memory that is used to compare with each segment of data 
% (currently proposing each 2min)

% The algorithm used for doing such comparison will be PV correlation.
% Therefore, you could plot PV memory index as a function of time, and
% maybe calculate the rate of change of PV by fitting the points to a
% curve. 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Perform the analysis
% [Session1_Memory_Development_MI_Shuffle, Session1_Memory_Development_MI, Session2_Memory_Development_MI_Shuffle, Session2_Memory_Development_MI] = Memory_Development_Analysis(Behavior_Tracking1, Behavior_Tracking2, Cell_Events1, Cell_Events2, PlaceCell_Summary1, PlaceCell_Summary2, Coactive_Num_Place_Cell, ds, ap);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Optional Graphing
% % figure 
% % subplot(1,2,1)
% % rec_time_in_min = ds.Total_T_Recording/60;
% % T_segment = rec_time_in_min/ap.Temporal_Division;
% % Time_Table = zeros(1, (rec_time_in_min/T_segment));
% % Time_Table(1,1) = T_segment;
% % for n = 1:(size(Time_Table,2)-1)
% %     Time_Table(n+1) = Time_Table(n) + T_segment;
% % end
% % scatter(Time_Table, Session1_Memory_Development_MI, 70, 'filled', 'g');
% % xlim([0 rec_time_in_min])
% % ylim([-0.1 1])
% % hold on
% % Fitted_Curve1 = fit(Time_Table(:), Session1_Memory_Development_MI(:), 'poly2');
% % plot(Fitted_Curve1, Time_Table, Session1_Memory_Development_MI_Shuffle, '--*k');
% % xlabel('Time(min)')
% % ylabel('PV Memory Index')
% % title('Session1 Memory Development & Shuffle')
% % legend('off')
% % set(gca,'FontSize',12);
% % hold off
% % 
% % subplot(1,2,2)
% % scatter(Time_Table, Session2_Memory_Development_MI, 70, 'filled', 'g');
% % xlim([0 rec_time_in_min])
% % ylim([-0.1 1])
% % hold on
% % Fitted_Curve2 = fit(Time_Table(:), Session2_Memory_Development_MI(:), 'poly2');
% % plot(Fitted_Curve2, Time_Table, Session2_Memory_Development_MI_Shuffle, '--*k');
% % xlabel('Time(min)')
% % ylabel('PV Memory Index')
% % title('Session2 Memory Development & Shuffle')
% % set(gca,'FontSize',12);
% % legend('off')
% % hold off

%% Finally, after all analysis is done, bundle results variables into struct Data_Summary
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
    'Num_Cell_Active_Session1', Num_Cell_Active_Session1,...
    'Num_Cell_Active_Session2', Num_Cell_Active_Session2,...
    'Num_Cell_Active_During_Both_Session', Num_Cell_Active_During_Both_Session,...
    'Total_Cell_Active', Total_Cell_Active,...
    'Coactive_Cell_Proportion', Coactive_Cell_Proportion,...
    'Session1_ONLY_Cell_Proportion', Session1_ONLY_Cell_Proportion,...
    'Session2_ONLY_Cell_Proportion', Session2_ONLY_Cell_Proportion,...
    'Cell_Recurring_Probability', Cell_Recurring_Probability,...
    'Session1_Event_Rates', Session1_Event_Rates,...
    'Session2_Event_Rates', Session2_Event_Rates,...
    'Session1_Mean_Ca_Event_Rate', Session1_Mean_Ca_Event_Rate,...
    'Session2_Mean_Ca_Event_Rate', Session2_Mean_Ca_Event_Rate,...
    'BinCounts_Session1_Event_Rate_Proportion', BinCounts_Session1_Event_Rate_Proportion,...
    'BinCounts_Session2_Event_Rate_Proportion', BinCounts_Session2_Event_Rate_Proportion,...
    'Event_Rate_Correlation', Event_Rate_Correlation,...
    'Session1_Event_Amps', Session1_Event_Amps,...
    'Session2_Event_Amps', Session2_Event_Amps,...
    'Session1_Event_Amps_Mean', Session1_Event_Amps_Mean,...
    'Session2_Event_Amps_Mean', Session2_Event_Amps_Mean,...
    'Mobility_Session1', Mobility_Session1,...
    'Mobility_Session2', Mobility_Session2,...
    'Cell_Active_Both_Session_Reorganize1', Cell_Active_Both_Session_Reorganize1,...
    'Cell_Active_Both_Session_Reorganize2', Cell_Active_Both_Session_Reorganize2,...
    'Behavior_Tracking1', Behavior_Tracking1,...
    'Behavior_Tracking2', Behavior_Tracking2,...
    'Total_Num_Place_Cell', Total_Num_Place_Cell,...
    'Coactive_Num_Place_Cell', Coactive_Num_Place_Cell,...
    'Smoothed_Ca_Rate_Map1_ForGraph', Smoothed_Ca_Rate_Map1_ForGraph,...
    'Smoothed_Ca_Rate_Map2_ForGraph', Smoothed_Ca_Rate_Map2_ForGraph,...
    'Shuffled_Rate_Map_Corr_values', Shuffled_Rate_Map_Corr_values,...
    'Mean_Shuffled_Rate_Map_Corr', Mean_Shuffled_Rate_Map_Corr,...
    'Rate_Map_Corr', Rate_Map_Corr,...
    'Shuffled_RateMap_Corr_Memory_Index', Shuffled_RateMap_Corr_Memory_Index,...
    'RateMap_Corr_Memory_Index', RateMap_Corr_Memory_Index,...
    'RateMap_Corr_Coherent_Rot', RateMap_Corr_Coherent_Rot,...
    'Shuffled_RateMapCorr_Proportion', Shuffled_Rate_Map_Correlation_Proportion,...
    'RateMapCorr_Proportion', Rate_Map_Correlation_Proportion,...
    'Shuffled_PV_Dist', Shuffled_PV_Dist,...
    'Shuffled_PV_Mean_Dist', Shuffled_PV_Mean_Dist,...
    'PV_Corr_Values', PV_Corr_Values,...
    'PV_Corr_Memory_Index', PV_Corr_Memory_Index,...
    'Shuffled_PV_Corr_Memory_Index', Shuffled_PV_Corr_Memory_Index,...
    'PV_Corr_Coherent_Rot', PV_Corr_Coherent_Rot,...
    'PV_Proportion', PV_Proportion,...
    'Shuffled_PV_Proportion', Shuffled_PV_Proportion,...
    'Mobility_Session1_Pass', Mobility_Session1_Pass,...
    'Mobility_Session2_Pass', Mobility_Session2_Pass);

% %     'Session1_Memory_Development_MI_Shuffle', Session1_Memory_Development_MI_Shuffle,...
% %     'Session1_Memory_Development_MI', Session1_Memory_Development_MI,...
% %     'Session2_Memory_Development_MI_Shuffle', Session2_Memory_Development_MI_Shuffle,...
% %     'Session2_Memory_Development_MI', Session2_Memory_Development_MI...


end
