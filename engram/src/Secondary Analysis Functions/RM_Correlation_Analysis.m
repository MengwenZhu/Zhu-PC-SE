function [Rate_Map_Corr, RateMap_Corr_Memory_Index, RateMap_Corr_Coherent_Rot, Rate_Map_Correlation_Proportion, Shuffled_Rate_Map_Correlation_Proportion]=RM_Correlation_Analysis(Coactive_Num_Place_Cell, Shuffled_Rate_Map_Corr_values, Cell_Active_Both_Session_Smoothed_Ca_Event_Rate_Map1, Cell_Active_Both_Session_Smoothed_Ca_Event_Rate_Map2, PlaceCell_Summary1, PlaceCell_Summary2, ap)

% Reject non-place cells if and only if the total number of place cells
% observed during both sessions exceeds a certain number (if there are too 
% few place cells we will use all of them for PV and RM correlation analysis)
if Coactive_Num_Place_Cell > ap.Num_Place_Cell_Threshold
    for n=1:size(Cell_Active_Both_Session_Smoothed_Ca_Event_Rate_Map1,3)
        if PlaceCell_Summary1(n,1) == 0 && PlaceCell_Summary2(n,1) == 0
            Cell_Active_Both_Session_Smoothed_Ca_Event_Rate_Map1(:,:,n)=0;
        end
    end
    for n=1:size(Cell_Active_Both_Session_Smoothed_Ca_Event_Rate_Map2,3)
        if PlaceCell_Summary1(n,1) == 0 && PlaceCell_Summary2(n,1) == 0
            Cell_Active_Both_Session_Smoothed_Ca_Event_Rate_Map2(:,:,n)=0;
        end
    end
end

% Calculate RM Correlation for Difference Rotation Angles 
% Then compute the RM Correlation Memory Index & Angle of Coherent Place Fields Rotation
% Angles: 0, 90, 180, and 270 degrees

% Notes: The reason that we rotate one session of rate map and calculate four distributions
% is that Kinsky 2018 paper examined a phenomenon called coherent rotation of place fields, which
% explains the phenomenon that an animal might be disoriented in a space but still remember it as the same 
% spatial context (so that the majority place fields rotate by a fixed angle); in our case, 
% since we are using square arena, we only need to consider angles of 90 degree multiples

% Calculate RM Correlation at different map rotation angles
Rate_Map_Corr=nan(size(Cell_Active_Both_Session_Smoothed_Ca_Event_Rate_Map1,3),4);
for n=1:size(Cell_Active_Both_Session_Smoothed_Ca_Event_Rate_Map1,3)
    Rate_Map_Corr(n,1)=XCorrCalc(Cell_Active_Both_Session_Smoothed_Ca_Event_Rate_Map1(:,:,n), Cell_Active_Both_Session_Smoothed_Ca_Event_Rate_Map2(:,:,n)); % no rotation
end

for n=1:size(Cell_Active_Both_Session_Smoothed_Ca_Event_Rate_Map1,3)
    Rate_Map_Corr(n,2)=XCorrCalc(Cell_Active_Both_Session_Smoothed_Ca_Event_Rate_Map1(:,:,n), rot90(Cell_Active_Both_Session_Smoothed_Ca_Event_Rate_Map2(:,:,n),1)); % 90 degrees rotation
end

for n=1:size(Cell_Active_Both_Session_Smoothed_Ca_Event_Rate_Map1,3)
    Rate_Map_Corr(n,3)=XCorrCalc(Cell_Active_Both_Session_Smoothed_Ca_Event_Rate_Map1(:,:,n), rot90(Cell_Active_Both_Session_Smoothed_Ca_Event_Rate_Map2(:,:,n),2)); % 180 degrees rotation
end

for n=1:size(Cell_Active_Both_Session_Smoothed_Ca_Event_Rate_Map1,3)
    Rate_Map_Corr(n,4)=XCorrCalc(Cell_Active_Both_Session_Smoothed_Ca_Event_Rate_Map1(:,:,n), rot90(Cell_Active_Both_Session_Smoothed_Ca_Event_Rate_Map2(:,:,n),3)); % 270 degrees rotation
end

% Calculate the mean rate map correlation value for each rotation angle
Mean_Rate_Map_Corr=nan(1,4);
Mean_Rate_Map_Corr(1,1)=mean(Rate_Map_Corr(:,1),'omitnan');
Mean_Rate_Map_Corr(1,2)=mean(Rate_Map_Corr(:,2),'omitnan');
Mean_Rate_Map_Corr(1,3)=mean(Rate_Map_Corr(:,3),'omitnan');
Mean_Rate_Map_Corr(1,4)=mean(Rate_Map_Corr(:,4),'omitnan');

% Calculate rate map correlation memory index
if max(Mean_Rate_Map_Corr,[],'all')>=ap.Ca_Rate_Map_Corr_Threshold 
   RateMap_Corr_Memory_Index=max(Mean_Rate_Map_Corr,[],'all'); 
else 
   RateMap_Corr_Memory_Index=Mean_Rate_Map_Corr(1,1); 
end

% Determine angle of coherent rotation
RateMap_Corr_Coherent_Rot=0; 
for n=1:4
    if Mean_Rate_Map_Corr(1,n)==max(Mean_Rate_Map_Corr,[],'all')
        RateMap_Corr_Coherent_Rot=(n-1)*90; 
    end
end
if RateMap_Corr_Memory_Index<ap.Ca_Rate_Map_Corr_Threshold
    RateMap_Corr_Coherent_Rot=0; 
end


% Binning the real RM Correlation Distribution
BinCounts_RateMapCorr = histcounts(Rate_Map_Corr(:,(RateMap_Corr_Coherent_Rot/90)+1),'BinEdges',(ap.Ca_Rate_Map_Corr_Lower_Bound:ap.Ca_Rate_Map_Corr_Interval:ap.Ca_Rate_Map_Corr_Upper_Bound));
Rate_Map_Correlation_Proportion=BinCounts_RateMapCorr/sum(BinCounts_RateMapCorr);

% Binning the Shuffled RM Correlation Distribution
BinCounts_RateMap_Shuffle_h = histcounts(Shuffled_Rate_Map_Corr_values,'BinEdges',(ap.Ca_Rate_Map_Corr_Lower_Bound:ap.Ca_Rate_Map_Corr_Interval:ap.Ca_Rate_Map_Corr_Upper_Bound));
BinCounts_RateMapCorr_Shuffle_Mean1=BinCounts_RateMap_Shuffle_h/(size(Shuffled_Rate_Map_Corr_values,3));
Shuffled_Rate_Map_Correlation_Proportion=BinCounts_RateMapCorr_Shuffle_Mean1/sum(BinCounts_RateMapCorr_Shuffle_Mean1);

% Optional Graphing
% % bar((ap.Ca_Rate_Map_Corr_Lower_Bound+0.5*ap.Ca_Rate_Map_Corr_Interval:ap.Ca_Rate_Map_Corr_Interval:ap.Ca_Rate_Map_Corr_Upper_Bound),Rate_Map_Correlation_Proportion,1,'b','FaceAlpha',0.6);
% % ylim([0 0.3]);
% % xlabel('Rate Map Correlation values', 'FontSize', 11);
% % ylabel('proportion of cells', 'FontSize', 11);
% % hold on;
% % plot((ap.Ca_Rate_Map_Corr_Lower_Bound+ap.Ca_Rate_Map_Corr_Interval:ap.Ca_Rate_Map_Corr_Interval:ap.Ca_Rate_Map_Corr_Upper_Bound),Shuffled_Rate_Map_Correlation_Proportion,'--k','lineWidth',1.5);
% % hold off;
% % legend('Proportion of cells','Shuffle');

end