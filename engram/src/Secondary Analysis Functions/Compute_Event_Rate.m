function [Session1_Event_Rates] = Compute_Event_Rate(Inscopix, Start_T, ds)

Cell_Events = Inscopix;
Cell_Events(:,2) = Cell_Events(:,2)+1;
Cell_Events(:,1) = Cell_Events(:,1)-Start_T;

% Check if any timestamp is outside of what we expect
if any(Cell_Events(:,1)<0) || any(Cell_Events(:,1)>ds.Total_T_Recording)
    disp('corrected timestamps are outside of expected time range')
else
    disp('corrected timestamps are within expected time range, good to go')
end

% Compute Calcium Event Rate Distribution for Each Recording
Session1_Event_Rates=nan(max(unique(Cell_Events(:,2))),1);
for n=(unique(Cell_Events(:,2)))'
    Session1_Event_Rates(n,1)=sum(Cell_Events(:,2)==n)/ds.Total_T_Recording; 
end

end