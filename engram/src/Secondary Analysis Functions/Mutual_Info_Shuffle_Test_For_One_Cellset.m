function Data_Summary = Mutual_Info_Shuffle_Test_For_One_Cellset(Inscopix, Noldus, Start_T, META, ds, ap)
%% Organize the raw input data

NoldusData = Noldus(1:(ds.Total_T_Recording/(1/ds.Noldus_Sampling_Frequency))+1,2:4); % extract useful Noldus data 
InscopixData = Inscopix; 
InscopixData(:,2) = Inscopix(:,2)+1;  % make cell names from 1 to N
InscopixData(:,1) = InscopixData(:,1)-Start_T; % correct the timestamps of the raw data

%% Calcium event amplitude correction

% adjust the amplitude of calcium events
InscopixData(:,3) = ap.Amp_factor_placeness*InscopixData(:,3);
% Additional Optional Step: Exclude Events with Amplitude under a Certain Threshold
for n = 1:size(InscopixData,1)
    if InscopixData(n,3)<=ap.Amp_Exclusion_Threshold
        InscopixData(n,:)=nan;
    end
end
InscopixData(isnan(InscopixData)) = [];
InscopixData = reshape(InscopixData, [], 3);



for n=1:size(InscopixData,1)
    if InscopixData(n,3)<=1 
        InscopixData(n,3)=1; 
    elseif InscopixData(n,3)>1 && rem(InscopixData(n,3),1)<0.5
       InscopixData(n,3)=InscopixData(n,3)-rem(InscopixData(n,3),1); 
    else
       InscopixData(n,3)=InscopixData(n,3)+(1-rem(InscopixData(n,3),1)); 
    end
end

R1=reshape(InscopixData(:,1),[1,1,size(InscopixData,1)]);
R2=reshape(InscopixData(:,2),[1,1,size(InscopixData,1)]);
R3=reshape(InscopixData(:,3),[1,1,size(InscopixData,1)]);

CatR=cat(2,R1,R2,R3);

CorrectedCat=nan(max(InscopixData(:,3),[],'all'), size(CatR,2), size(CatR,3));
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

InscopixData=cat(2,C1,C2,C3); 

% fix bugs of going beyond recording time:
for n=1:size(InscopixData,1)
    if InscopixData(n,1)>ds.Total_T_Recording
        InscopixData(n,:)=nan;
    end
end
InscopixData(isnan(InscopixData))=[];
InscopixData=reshape(InscopixData,[],3);

%% Process Noldus behavioral data

% Smoothing by convolution with a triangular filter 
smooth_win = triang(round(ds.Noldus_Sampling_Frequency * ap.Behavior_Filter_Size));
smooth_win = smooth_win / sum(smooth_win);

Displacement=zeros(size(NoldusData,1),1);
for n=1:size(NoldusData,1)
    if n+1<size(NoldusData,1)
        Displacement(n,1)=sqrt((NoldusData(n+1,2)-NoldusData(n,2))^2+(NoldusData(n+1,3)-NoldusData(n,3))^2); % calculate displacement at each 0.04sec interval
    end
end
% fix bugs: sometimes displacement calculate is nan because very occassinoal
% Noldus tracking errors
Displacement(isnan(Displacement))=0;
Displacement(isinf(Displacement))=0;
Speed=Displacement./(1/ds.Noldus_Sampling_Frequency); % calculate the speed during each 0.04 time interval
% smooth trace
Speed = conv(Speed, smooth_win, 'same');

% Compute the fraction of time that mouse actively explores for session1&2
Mobility=(sum(Speed>ap.Low_Mobility_Threshold))/(size(Speed,1)); 

% give an indication of whether the recording should be included in
% statistical analysis (mobility < certain value should not be included
% within the summary of place cell & engram statistics)
if Mobility >= 0.15
    Mobility_Pass = 1;
else
    Mobility_Pass = 0;
end

NoldusData=cat(2, NoldusData, Speed); % finish up analysis by concatenating this info to Noldus files

%% Generate Corrected_Occupancy_Map 

% AM session map boundaries
A=min(NoldusData(:,2),[],'omitnan'); % calculate map boundaries
B=max(NoldusData(:,2),[],'omitnan');
C=min(NoldusData(:,3),[],'omitnan');
D=max(NoldusData(:,3),[],'omitnan');

Map_LowX_Bound=fix(A)-ap.Arena_Add_Edge; % To fully depict the arena within a matrix without leaving out data points
Map_HighX_Bound=fix(B)+ap.Arena_Add_Edge;
Map_LowY_Bound=fix(C)-ap.Arena_Add_Edge;
Map_HighY_Bound=fix(D)+ap.Arena_Add_Edge;

Occupancy_Map=hist3(NoldusData(:,2:3),'ctrs',{Map_LowX_Bound:(Map_HighX_Bound-Map_LowX_Bound)/ap.Map_Division_For_MI:Map_HighX_Bound  Map_LowY_Bound:(Map_HighY_Bound-Map_LowY_Bound)/ap.Map_Division_For_MI:Map_HighY_Bound}); 
Occupancy_Map=rot90(Occupancy_Map,1); % Rotation correction (to match the experimental setup)
Corrected_Occupancy_Map=Occupancy_Map.*(1/ds.Noldus_Sampling_Frequency); % convert the values in matrix to be measured in seconds (each position correspond to 0.04sec)

%% Find the XY-positions of the animal given a calcium event 

% for AM session find XY positions of calcium events
[InscopixData]=Find_XY_Position_Single_Session(InscopixData, NoldusData, ds.Noldus_Sampling_Frequency);

% exclude events during immobility as well
for n=1:size(InscopixData,1)
    if InscopixData(n,6)>ap.High_Mobility_Threshold || InscopixData(n,6)<ap.Low_Mobility_Threshold
        InscopixData(n,1:6)=nan(1,6); % exclude events during imobility
    end
end
InscopixData(any(isnan(InscopixData),2),:) = []; % remove nan

%% Generate Ca_Event_Map & Ca_Event_Rate_Map & Normalized_Ca_Event_Rate_Map

n=(unique(InscopixData(:,2)))';
MaxFiring=max(sum(InscopixData(:,2)==n)); 

% Reorganize data in to matrix for generating event maps for each neuron
RawFiring3DMatrix=nan(MaxFiring, size(InscopixData,2), size(unique(InscopixData(:,2)),1));
for n=(unique(InscopixData(:,2)))'
    y=sum(InscopixData(:,2)==n);
    RawFiring3DMatrix(:,:,n)=reshape([InscopixData(cat(1,InscopixData(:,2))==n,:);NaN(MaxFiring-y,size(InscopixData,2))],MaxFiring,size(InscopixData,2)); % reorganize the raw firing data for analysis
end
RawFiring3DMatrix(RawFiring3DMatrix==0)=nan;


% Ca event map is the number of Ca events happened within each spatial bin of the arena, for each neuron
Ca_Event_Map=nan(size(Occupancy_Map,1),size(Occupancy_Map,2),size(RawFiring3DMatrix,3));
for k=1:size(RawFiring3DMatrix,3)
    Ca_Event_Map(:,:,k)=hist3(RawFiring3DMatrix(:,4:5,k),'ctrs',{Map_LowX_Bound:(Map_HighX_Bound-Map_LowX_Bound)/ap.Map_Division_For_MI:Map_HighX_Bound  Map_LowY_Bound:(Map_HighY_Bound-Map_LowY_Bound)/ap.Map_Division_For_MI:Map_HighY_Bound});
end
Ca_Event_Map=rot90(Ca_Event_Map,1); %Rotational Correction


% Ca Event Rate Map (events/sec)
Ca_Event_Rate_Map=nan(size(Occupancy_Map,1),size(Occupancy_Map,2),size(RawFiring3DMatrix,3));
for k=1:size(RawFiring3DMatrix,3)
    Ca_Event_Rate_Map(:,:,k)=(Ca_Event_Map(:,:,k))./Occupancy_Map;
end
Ca_Event_Rate_Map(isnan(Ca_Event_Rate_Map)) = 0;


% Normalized Ca Event Rate Map (prob events/sec)
Normalized_Ca_Event_Rate_Map=nan(size(Occupancy_Map,1),size(Occupancy_Map,2),size(RawFiring3DMatrix,3));
for k=1:size(RawFiring3DMatrix,3)
    Normalized_Ca_Event_Rate_Map(:,:,k)=(Ca_Event_Rate_Map(:,:,k))./sum(Ca_Event_Rate_Map(:,:,k),'all');
end

%% Calculate Mutual Information 

% this needs the 'calc_mutual_information' function
% calculate MI values for the AM session cells
MI_values=zeros(size(Ca_Event_Map,3),1);
for n = (unique(InscopixData(:,2)))'
    [MI_values(n,1)]= calc_mutual_information(Normalized_Ca_Event_Rate_Map(:,:,n), Corrected_Occupancy_Map);
end
% get rid of useless info (if cell has any firing it cannot have MI == 0, always a positive number)
MI_values(MI_values==0)=nan; 

%% Shuffle Mutual Information

tic
MI_Shuffled_values=nan(size(MI_values,1),1,ap.Num_MI_Shuffle);

% shuffle MI values for AM session cells
for n=1:ap.Num_MI_Shuffle
    [MI_Shuffled_values(:,1,n)]=Shuffle_MI_Single_Session(MI_values, InscopixData, Corrected_Occupancy_Map, NoldusData, ds, ap);
end

toc

%% Get MI distribution calculated as proportion of cells at each interval of MI values (for later plotting)

% replaced histogram with histcounts to speed up code
BinCounts_h=histcounts(MI_Shuffled_values, (ap.MI_Lower_Bound:ap.MI_Interval:ap.MI_Upper_Bound));
BinCounts_RateMapCorr_Shuffle_Mean=BinCounts_h/(size(MI_Shuffled_values,3));
ShuffledMI_Proportion=BinCounts_RateMapCorr_Shuffle_Mean/sum(BinCounts_RateMapCorr_Shuffle_Mean); % get proportion of cells that fall within each bin of MI values for null distribution

% replaced histogram with histcounts as well
BinCounts_MI=histcounts(MI_values(:), (ap.MI_Lower_Bound:ap.MI_Interval:ap.MI_Upper_Bound));
MI_Proportion=BinCounts_MI/sum(BinCounts_MI); % get proportion of cells that fall within each bin of MI values for real data

% make graphs
% % bar((ap.MI_Lower_Bound+0.5*ap.MI_Interval:ap.MI_Interval:ap.MI_Upper_Bound),MI_Proportion,1,'b','FaceAlpha',0.6);
% % ylim([0 0.5]);
% % xlabel('Mutual Information Values', 'FontSize', 11);
% % ylabel('proportion of cells', 'FontSize', 11);
% % hold on;
% % plot((ap.MI_Lower_Bound+ap.MI_Interval:ap.MI_Interval:ap.MI_Upper_Bound),ShuffledMI_Proportion,'--k','lineWidth',1.5);
% % hold off;


%% Calculate MI-shuffle test p-values

MICat=cat(3,MI_values,MI_Shuffled_values); 

MI_p_values=nan(size(MICat,1),1);
for n = (unique(InscopixData(:,2)))'
      MI_p_values(n,1)=(sum(MICat(n,1,2:end)>MICat(n,1,1)))/(size(MICat,3)-1);
end

%% Plot the distribution of MI shuffle test p-values (for later making figures)

% replaced histogram with hitcounts
BinCounts_MI_pvalues=histcounts(MI_p_values, (ap.MI_pvalue_Lower_Bound:ap.MI_pvalue_Interval:ap.MI_pvalue_Upper_Bound));
MI_pvalues_proportion=BinCounts_MI_pvalues/sum(BinCounts_MI_pvalues);
% find the proportion of cells that have p(MI) < 0.05
PlaceCell_Proportion=MI_pvalues_proportion(1,1);
% make graphs
% % bar((ap.MI_pvalue_Lower_Bound+0.5*ap.MI_pvalue_Interval:ap.MI_pvalue_Interval:ap.MI_pvalue_Upper_Bound),MI_pvalues_proportion,1,'b','FaceAlpha',0.6);
% % ylim([0 0.3]);
% % xlabel('mutual info shuffle test p-values', 'FontSize', 11);
% % ylabel('proportion of cells', 'FontSize', 11);

%% Calculate p-values for cell events that represent true randomness (totally random firing pattern vs. null MI distribution)
% This is done by spatially shuffle the normalized calcium rate map
% failed, didn't work, might have to figure out in the future

% % K = Normalized_Ca_Event_Rate_Map;
% % 
% % % spatially shuffle the rate maps
% % Spatially_Shuffled_Rate_Map = zeros(1, size(Normalized_Ca_Event_Rate_Map,1)*size(Normalized_Ca_Event_Rate_Map,2), size(Normalized_Ca_Event_Rate_Map,3));
% % for n = 1:size(K,3)
% %     Spatially_Shuffled_Rate_Map(:,:,n) = Matrix_Perm_Shuffling(K(:,:,n));
% % end
% % 
% % % calculate spatially shuffled rate map MI values
% % True_Shuffled_MI_Values = nan(size(Normalized_Ca_Event_Rate_Map,3),1);
% % for n = 1:size(Spatially_Shuffled_Rate_Map,3)
% %     [True_Shuffled_MI_Values(n,1)] = calc_mutual_information(Spatially_Shuffled_Rate_Map(:,:,n),Corrected_Occupancy_Map(:));
% % end
% % 
% % % calculate spatially shuffled MI p-values
% % MICat=cat(3,True_Shuffled_MI_Values,MI_Shuffled_values); 
% % Shuffled_MI_p_values=nan(size(MICat,1),1);
% % for n = (unique(InscopixData(:,2)))'
% %       Shuffled_MI_p_values(n,1)=(sum(MICat(n,1,2:end)>MICat(n,1,1)))/(size(MICat,3)-1);
% % end

%% Make a matrix that tells whether a cell is a place cell or not (1=place cell, 0=non-place cell, nan=cell not detected)
% This is the most essential part of this function, which will be used for
% later PV and Rate_Map_Corr analysis to select only place cells to
% participate (and potentially some other analyses)

Place_Cell_Summary = nan(size(MI_p_values,1),1);
for n = 1:size(MI_p_values,1)
    if ~isnan(MI_p_values(n,1)) && MI_p_values(n,1)<ap.MI_p_value_Threshold
        Place_Cell_Summary(n,1) = 1;
    elseif ~isnan(MI_p_values(n,1)) && MI_p_values(n,1)>=ap.MI_p_value_Threshold
        Place_Cell_Summary(n,1) = 0;
    else
        Place_Cell_Summary(n,1) = nan;
    end
end


%% Organize all outputs into a structure

Data_Summary = struct(...
    'animalName', META.animalName, ...
    'genotype', META.genotype, ...
    'exptDate1', META.exptDate1, ...
    'exptDate2', META.exptDate2, ...
    'session', [], ...
    'drug', META.drug, ...
    'session1_dose', META.session1_dose, ...
    'session2_dose', META.session2_dose, ...
    'session1_context', META.session1_context,...
    'session2_context', META.session2_context,...
    'exptParadigm', META.exptParadigm, ...
    'Place_Cell_Summary', Place_Cell_Summary, ...
    'MI_values', MI_values, ...
    'MI_Shuffled_values', MI_Shuffled_values, ...
    'ShuffledMI_Proportion', ShuffledMI_Proportion, ...
    'MI_Proportion', MI_Proportion, ...
    'MI_p_values', MI_p_values, ...
    'MI_pvalues_proportion', MI_pvalues_proportion,...
    'PlaceCell_Proportion', PlaceCell_Proportion,...
    'Mobility', Mobility,...
    'Mobility_Pass', Mobility_Pass);

end