function [Smoothed_Ca_Rate_Map1_ForGraph, Smoothed_Ca_Rate_Map2_ForGraph]=Generate_Smoothed_Ca_Event_Maps_For_Graphing(RawFiringData1, RawFiringData2, RawXYData1, RawXYData2, ds, ap)
%% Step 1: Generate Occupancy_Map for Session1 and Session2 

A=min(RawXYData1(:,2),[],'omitnan');
B=max(RawXYData1(:,2),[],'omitnan');
C=min(RawXYData1(:,3),[],'omitnan');
D=max(RawXYData1(:,3),[],'omitnan');

Map1_LowX_Bound=fix(A)-ap.Arena_Add_Edge;
Map1_HighX_Bound=fix(B)+ap.Arena_Add_Edge;
Map1_LowY_Bound=fix(C)-ap.Arena_Add_Edge;
Map1_HighY_Bound=fix(D)+ap.Arena_Add_Edge;

Occupancy_Map1=hist3(RawXYData1(:,2:3),'ctrs',{Map1_LowX_Bound:(Map1_HighX_Bound-Map1_LowX_Bound)/ap.Ca_Map_Graph_Division:Map1_HighX_Bound  Map1_LowY_Bound:(Map1_HighY_Bound-Map1_LowY_Bound)/ap.Ca_Map_Graph_Division:Map1_HighY_Bound});

Occupancy_Map1=rot90(Occupancy_Map1,1); % Rotation correction (to match the behavioral tracking videos)

Corrected_Occupancy_Map1=Occupancy_Map1.*(1/ds.Noldus_Sampling_Frequency); % convert the values in matrix to be measured in seconds

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


E=min(RawXYData2(:,2),[],'omitnan');
F=max(RawXYData2(:,2),[],'omitnan');
G=min(RawXYData2(:,3),[],'omitnan');
H=max(RawXYData2(:,3),[],'omitnan');

Map2_LowX_Bound=fix(E)-ap.Arena_Add_Edge;
Map2_HighX_Bound=fix(F)+ap.Arena_Add_Edge;
Map2_LowY_Bound=fix(G)-ap.Arena_Add_Edge;
Map2_HighY_Bound=fix(H)+ap.Arena_Add_Edge;

Occupancy_Map2=hist3(RawXYData2(:,2:3),'ctrs',{Map2_LowX_Bound:(Map2_HighX_Bound-Map2_LowX_Bound)/ap.Ca_Map_Graph_Division:Map2_HighX_Bound  Map2_LowY_Bound:(Map2_HighY_Bound-Map2_LowY_Bound)/ap.Ca_Map_Graph_Division:Map2_HighY_Bound});

Occupancy_Map2=rot90(Occupancy_Map2,1);

Corrected_Occupancy_Map2=Occupancy_Map2.*(1/ds.Noldus_Sampling_Frequency);

%% Step 2: Generate Ca_Event_Map for Session1 and Session2 

n=(unique(RawFiringData1(:,2)))';
MaxFiring1=max(sum(RawFiringData1(:,2)==n));

RawFiring3DMatrix1=nan(MaxFiring1, size(RawFiringData1,2), size(unique(RawFiringData1(:,2)),1));
for n=(unique(RawFiringData1(:,2)))'
    y=sum(RawFiringData1(:,2)==n);
    RawFiring3DMatrix1(:,:,n)=reshape([RawFiringData1(cat(1,RawFiringData1(:,2))==n,:);NaN(MaxFiring1-y,size(RawFiringData1,2))],MaxFiring1,size(RawFiringData1,2)); % reorganize the raw firing data for analysis
end
RawFiring3DMatrix1(RawFiring3DMatrix1==0)=nan;


Ca_Event_Map1=nan(size(Occupancy_Map1,1),size(Occupancy_Map1,2),size(RawFiring3DMatrix1,3));
for k=1:size(RawFiring3DMatrix1,3)
    Ca_Event_Map1(:,:,k)=hist3(RawFiring3DMatrix1(:,4:5,k),'ctrs',{Map1_LowX_Bound:(Map1_HighX_Bound-Map1_LowX_Bound)/ap.Ca_Map_Graph_Division:Map1_HighX_Bound  Map1_LowY_Bound:(Map1_HighY_Bound-Map1_LowY_Bound)/ap.Ca_Map_Graph_Division:Map1_HighY_Bound}); % generate Ca_Event_Map, where each value per bin is the # of calcium events happened in that bin
end
Ca_Event_Map1=rot90(Ca_Event_Map1,1); %Rotational Correction



n=(unique(RawFiringData2(:,2)))';
MaxFiring2=max(sum(RawFiringData2(:,2)==n));

RawFiring3DMatrix2=nan(MaxFiring2, size(RawFiringData1,2), size(unique(RawFiringData2(:,2)),1));
for n=(unique(RawFiringData2(:,2)))'
    y=sum(RawFiringData2(:,2)==n);
    RawFiring3DMatrix2(:,:,n)=reshape([RawFiringData2(cat(1,RawFiringData2(:,2))==n,:);NaN(MaxFiring2-y,size(RawFiringData1,2))],MaxFiring2,size(RawFiringData1,2));
end


Ca_Event_Map2=nan(size(Occupancy_Map2,1),size(Occupancy_Map2,2),size(RawFiring3DMatrix2,3));
for k=1:size(RawFiring3DMatrix2,3)
    Ca_Event_Map2(:,:,k)=hist3(RawFiring3DMatrix2(:,4:5,k),'ctrs',{Map2_LowX_Bound:(Map2_HighX_Bound-Map2_LowX_Bound)/ap.Ca_Map_Graph_Division:Map2_HighX_Bound  Map2_LowY_Bound:(Map2_HighY_Bound-Map2_LowY_Bound)/ap.Ca_Map_Graph_Division:Map2_HighY_Bound}); 
end
Ca_Event_Map2=rot90(Ca_Event_Map2,1); %Rotational Correction



%% Step 3: Calculate the Unsmoothed_Ca_Rate_Map for Session1 and Session2 

q1=Ca_Event_Map1./Corrected_Occupancy_Map1;
Unsmoothed_Ca_Rate_Map1=nan(size(Ca_Event_Map1,1),size(Ca_Event_Map1,2),size(Ca_Event_Map1,3));
for n = (unique(RawFiringData1(:,2)))'
    Unsmoothed_Ca_Rate_Map1(:,:,n)=q1(:,:,n);
end
Unsmoothed_Ca_Rate_Map1(isnan(Unsmoothed_Ca_Rate_Map1))=0;


q2=Ca_Event_Map2./Corrected_Occupancy_Map2;
Unsmoothed_Ca_Rate_Map2=nan(size(Ca_Event_Map2,1),size(Ca_Event_Map2,2),size(Ca_Event_Map2,3));
for n = (unique(RawFiringData2(:,2)))'
    Unsmoothed_Ca_Rate_Map2(:,:,n)=q2(:,:,n);
end
Unsmoothed_Ca_Rate_Map2(isnan(Unsmoothed_Ca_Rate_Map2))=0;


% Keep size same
if size(Unsmoothed_Ca_Rate_Map1,3)<size(Unsmoothed_Ca_Rate_Map2,3)
    Unsmoothed_Ca_Rate_Map1(1:size(Unsmoothed_Ca_Rate_Map1,1),1:size(Unsmoothed_Ca_Rate_Map1,2),(end+1):size(Unsmoothed_Ca_Rate_Map2,3))=zeros(size(Unsmoothed_Ca_Rate_Map1,1),size(Unsmoothed_Ca_Rate_Map1,2),size(Unsmoothed_Ca_Rate_Map2,3)-size(Unsmoothed_Ca_Rate_Map1,3)); % keep size the same for shuffling convenience
end



%% Step 4: Calculate the Smoothed_Ca_Rate_Map_ForGraph for Session1 and Session2 (50-by-50 matrix)

Smoothed_Ca_Rate_Map1_ForGraph=zeros(size(Unsmoothed_Ca_Rate_Map1,1), size(Unsmoothed_Ca_Rate_Map1,2), size(Unsmoothed_Ca_Rate_Map1,3));
Smoothed_Ca_Rate_Map2_ForGraph=zeros(size(Unsmoothed_Ca_Rate_Map2,1), size(Unsmoothed_Ca_Rate_Map2,2), size(Unsmoothed_Ca_Rate_Map2,3));

for n=1:size(Unsmoothed_Ca_Rate_Map1,3)
       Smoothed_Ca_Rate_Map1_ForGraph(:,:,n)=filter2(ap.Gaussian_Filter_For_Graph,Unsmoothed_Ca_Rate_Map1(:,:,n));
end

for n=1:size(Unsmoothed_Ca_Rate_Map2,3)
      Smoothed_Ca_Rate_Map2_ForGraph(:,:,n)=filter2(ap.Gaussian_Filter_For_Graph,Unsmoothed_Ca_Rate_Map2(:,:,n));
end

Smoothed_Ca_Rate_Map1_ForGraph(isnan(Smoothed_Ca_Rate_Map1_ForGraph))=0;
Smoothed_Ca_Rate_Map1_ForGraph(isinf(Smoothed_Ca_Rate_Map1_ForGraph))=0;
Smoothed_Ca_Rate_Map2_ForGraph(isnan(Smoothed_Ca_Rate_Map2_ForGraph))=0;
Smoothed_Ca_Rate_Map2_ForGraph(isinf(Smoothed_Ca_Rate_Map2_ForGraph))=0;


%% Scripts for generating graphs of firing positions on behavioral traces & smoothed Ca-event-rate heat maps

% % n=;
% % subplot(2,2,1);
% % plot(RawXYData1(:,2),RawXYData1(:,3));
% % xlim([min(RawXYData1(:,2))-1,max(RawXYData1(:,2)+1)]);
% % ylim([min(RawXYData1(:,3))-1,max(RawXYData1(:,3)+1)]);
% % hold on
% % scatter(RawFiring3DMatrix1(:,2,n),RawFiring3DMatrix1(:,3,n),20*RawFiring3DMatrix1(:,5,n),'filled');
% % % % xlabel('X-position');
% % % % ylabel('Y-position');
% % % % k=num2str(n);
% % % % name=cat(2,'Session1 Cell',k,' calcium events on mouse trace');
% % % % title(name);
% % hold off
% % 
% % 
% % subplot(2,2,2);
% % plot(RawXYData2(:,2),RawXYData2(:,3));
% % hold on
% % scatter(RawFiring3DMatrix2(:,2,n),RawFiring3DMatrix2(:,3,n),20*RawFiring3DMatrix2(:,5,n),'filled');
% % % % xlabel('X-position');
% % % % ylabel('Y-position');
% % % % k=num2str(n);
% % % % name=cat(2,'Session2 Cell',k,' calcium events on mouse trace');
% % % % title(name);
% % hold off
% % 
% % 
% % subplot(2,2,3);
% % heatmap(Smoothed_Ca_Rate_Map1_ForGraph(:,:,n),'CellLabelColor','none'); 
% % colormap jet; 
% % Ax = gca;
% % Ax.XDisplayLabels = nan(size(Ax.XDisplayData));
% % Ax.YDisplayLabels = nan(size(Ax.YDisplayData));
% % % % k=num2str(n);
% % % % name=cat(2,'Session1 Cell',k,' smoothed calcium event rate map');
% % % % title(name);
% % 
% % 
% % subplot(2,2,4);
% % heatmap(Smoothed_Ca_Rate_Map2_ForGraph(:,:,n),'CellLabelColor','none'); 
% % colormap jet;
% % Ax = gca;
% % Ax.XDisplayLabels = nan(size(Ax.XDisplayData));
% % Ax.YDisplayLabels = nan(size(Ax.YDisplayData));
% % % % k=num2str(n);
% % % % name=cat(2,'Session2 Cell',k,' smoothed calcium event rate map');
% % % % title(name);


end