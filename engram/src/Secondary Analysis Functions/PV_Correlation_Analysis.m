function [PV_Corr_Values,PV_Corr_Memory_Index,PV_Corr_Coherent_Rot,PV_Proportion,Shuffled_PV_Proportion]=PV_Correlation_Analysis(Num_Place_Cell, Shuffled_PV_Corr, BothSession_Smoothed_Ca_Rate_Map1, BothSession_Smoothed_Ca_Rate_Map2, PlaceCell_Summary1, PlaceCell_Summary2, ap)
%% Reject non-place cells if there are enough place cells

% calculate Ca map matrix size first to make code less messy:
size1=size(BothSession_Smoothed_Ca_Rate_Map1,1);
size2=size(BothSession_Smoothed_Ca_Rate_Map1,2);

% update to fix bugs
if Num_Place_Cell >= ap.Num_Place_Cell_Threshold
    for n=1:size(BothSession_Smoothed_Ca_Rate_Map1,3)
        if PlaceCell_Summary1(n,1) == 0 && PlaceCell_Summary2(n,1) == 0 
            BothSession_Smoothed_Ca_Rate_Map1(:,:,n)=0;
        end
    end
    for n=1:size(BothSession_Smoothed_Ca_Rate_Map2,3)
        if PlaceCell_Summary1(n,1) == 0 && PlaceCell_Summary2(n,1) == 0 
            BothSession_Smoothed_Ca_Rate_Map2(:,:,n)=0;
        end
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for n=1:size(BothSession_Smoothed_Ca_Rate_Map1,3)
    if BothSession_Smoothed_Ca_Rate_Map1(:,:,n)==zeros(size1,size2) %getting rid of some zero matrices (kept the cell names for other analyses to keep track of cell identities)
        BothSession_Smoothed_Ca_Rate_Map1(:,:,n)=nan;
    end
end
BothSession_Smoothed_Ca_Rate_Map1(isnan(BothSession_Smoothed_Ca_Rate_Map1))=[];
BothSession_Smoothed_Ca_Rate_Map1=reshape(BothSession_Smoothed_Ca_Rate_Map1,size1,size2,[]);



for n=1:size(BothSession_Smoothed_Ca_Rate_Map2,3)
    if BothSession_Smoothed_Ca_Rate_Map2(:,:,n)==zeros(size1,size2)
        BothSession_Smoothed_Ca_Rate_Map2(:,:,n)=nan;
    end
end
BothSession_Smoothed_Ca_Rate_Map2(isnan(BothSession_Smoothed_Ca_Rate_Map2))=[];
BothSession_Smoothed_Ca_Rate_Map2=reshape(BothSession_Smoothed_Ca_Rate_Map2,size1,size2,[]);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Keep sizes of the two matrices the same
if size(BothSession_Smoothed_Ca_Rate_Map1,3)>size(BothSession_Smoothed_Ca_Rate_Map2,3)
    BothSession_Smoothed_Ca_Rate_Map2(1:size(BothSession_Smoothed_Ca_Rate_Map2,1),1:size(BothSession_Smoothed_Ca_Rate_Map2,2),end+1:size(BothSession_Smoothed_Ca_Rate_Map1,3))=zeros(size(BothSession_Smoothed_Ca_Rate_Map2,1),size(BothSession_Smoothed_Ca_Rate_Map2,2),(size(BothSession_Smoothed_Ca_Rate_Map1,3)-size(BothSession_Smoothed_Ca_Rate_Map2,3)));
end

if size(BothSession_Smoothed_Ca_Rate_Map1,3)<size(BothSession_Smoothed_Ca_Rate_Map2,3)
    BothSession_Smoothed_Ca_Rate_Map1(1:size(BothSession_Smoothed_Ca_Rate_Map1,1),1:size(BothSession_Smoothed_Ca_Rate_Map1,2),end+1:size(BothSession_Smoothed_Ca_Rate_Map2,3))=zeros(size(BothSession_Smoothed_Ca_Rate_Map1,1),size(BothSession_Smoothed_Ca_Rate_Map1,2),(size(BothSession_Smoothed_Ca_Rate_Map2,3)-size(BothSession_Smoothed_Ca_Rate_Map1,3)));
end


%% Calculate PV_Correlation Values at 0, 90, 180, and 270 degrees rotation

PV_Corr_Values(1,:,1)=PV_Correlation(BothSession_Smoothed_Ca_Rate_Map1, BothSession_Smoothed_Ca_Rate_Map2);
PV_Corr_Values(1,:,2)=PV_Correlation(BothSession_Smoothed_Ca_Rate_Map1, rot90(BothSession_Smoothed_Ca_Rate_Map2,1));
PV_Corr_Values(1,:,3)=PV_Correlation(BothSession_Smoothed_Ca_Rate_Map1, rot90(BothSession_Smoothed_Ca_Rate_Map2,2));
PV_Corr_Values(1,:,4)=PV_Correlation(BothSession_Smoothed_Ca_Rate_Map1, rot90(BothSession_Smoothed_Ca_Rate_Map2,3));


%% Determine the angle of coherent rotation from PV_Correlation Analysis

Mean_PV_Corr=nan(1,4);
Mean_PV_Corr(1,1)=mean(PV_Corr_Values(1,:,1),'omitnan');
Mean_PV_Corr(1,2)=mean(PV_Corr_Values(1,:,2),'omitnan');
Mean_PV_Corr(1,3)=mean(PV_Corr_Values(1,:,3),'omitnan');
Mean_PV_Corr(1,4)=mean(PV_Corr_Values(1,:,4),'omitnan');


% calculate PV correlation memory index
if max(Mean_PV_Corr,[],'all')>=ap.PV_Corr_Threshold % if mean PV correlation is greater than 0.2, which is arbitrarily considered as the lowest threshold of having some memory
   PV_Corr_Memory_Index=max(Mean_PV_Corr,[],'all'); % use the maximum value in rate map correlation 
else 
   PV_Corr_Memory_Index=Mean_PV_Corr(1,1); % otherwise if all angles have no significant memory we just say the rotation is 0 degrees (no way to determine if there's no memory at all)
end


% determine angle of coherent rotation
PV_Corr_Coherent_Rot=0; % predefine variable to avoid errors
for n=1:4
    if Mean_PV_Corr(1,n)==max(Mean_PV_Corr,[],'all')
        PV_Corr_Coherent_Rot=(n-1)*90; % determine the coherent angle of rotation from rate map correlation (which angle generates the maximum value of mean rate map correlation)
    end
end

if PV_Corr_Memory_Index<ap.PV_Corr_Threshold
    PV_Corr_Coherent_Rot=0; % but if rate map corr memory index is less than 0.1 we will just keep it zero
end


%% Make a histogram of the distribution of PV correlation & shuffle (y-axis as proportion)

% replaced histogram function with hiscounts to avoid MATLAB crash
BinCounts_h = histcounts(Shuffled_PV_Corr,'BinEdges',(ap.PV_Corr_Lower_Bound:ap.PV_Corr_Interval:ap.PV_Corr_Upper_Bound));
BinCounts_h_Mean1=BinCounts_h/(size(Shuffled_PV_Corr,3));
Shuffled_PV_Proportion=BinCounts_h_Mean1/sum(BinCounts_h_Mean1);

BinCounts_PV = histcounts(PV_Corr_Values(1,:,(PV_Corr_Coherent_Rot/90)+1),'BinEdges',(ap.PV_Corr_Lower_Bound:ap.PV_Corr_Interval:ap.PV_Corr_Upper_Bound));
PV_Proportion=BinCounts_PV/sum(BinCounts_PV);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% graph
% % bar((Lower_Bound+0.5*Interval:Interval:Upper_Bound),PV_Proportion,1,'b','FaceAlpha',0.6);
% % ylim([0 0.3]);
% % xlabel('PV-correlation values', 'FontSize', 11);
% % ylabel('proportion of spatial bins', 'FontSize', 11);
% % hold on;
% % plot((Lower_Bound+Interval:Interval:Upper_Bound),Shuffled_PV_Proportion,'--k','lineWidth',1.5);
% % hold off;
% % legend('Proportion of spatial bins','Shuffle')

end