function [Behavior_Tracking, Cell_Events] = XDay_RawData_Preparation(Inscopix, Noldus, Start_T, ds, ap)

%This function preprocesses the raw data (Insocpix, Noldus, Start_T) that 
%are extracted for cross day analysis

%%
% Prep
Behavior_Tracking=Noldus(1:(ds.Total_T_Recording/(1/ds.Noldus_Sampling_Frequency))+1,2:4); 
Cell_Events = Inscopix;
Cell_Events(:,2) = Cell_Events(:,2)+1;
Cell_Events(:,1) = Cell_Events(:,1)-Start_T;

% Check if any timestamp is outside of what we expect
if any(Cell_Events(:,1)<0) || any(Cell_Events(:,1)>ds.Total_T_Recording)
    disp('corrected timestamps are outside of expected time range')
else
    disp('corrected timestamps are within expected time range, good to go')
end

%%
% Linearly adjust each Ca2+ event according to its peak amplitude, measured in dF/noise (refer to IDPS documentation for further details)

% adjust the amplitude of calcium events by a factor
Cell_Events(:,3) = ap.Amp_factor_PV*Cell_Events(:,3);

for n=1:size(Cell_Events,1)
    if Cell_Events(n,3)<=1
        Cell_Events(n,3)=1; 
    elseif Cell_Events(n,3)>1 && rem(Cell_Events(n,3),1)<0.5
       Cell_Events(n,3)=Cell_Events(n,3)-rem(Cell_Events(n,3),1); 
    else
       Cell_Events(n,3)=Cell_Events(n,3)+(1-rem(Cell_Events(n,3),1)); 
    end
end

R1=reshape(Cell_Events(:,1),[1,1,size(Cell_Events,1)]);
R2=reshape(Cell_Events(:,2),[1,1,size(Cell_Events,1)]);
R3=reshape(Cell_Events(:,3),[1,1,size(Cell_Events,1)]);

CatR=cat(2,R1,R2,R3);

CorrectedCat=nan(max(Cell_Events(:,3),[],'all'), size(CatR,2), size(CatR,3));
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

Cell_Events=cat(2,C1,C2,C3); 

% fix bugs of going beyond total recording time:
for n=1:size(Cell_Events,1)
    if Cell_Events(n,1)>ds.Total_T_Recording
        Cell_Events(n,:)=nan;
    end
end
Cell_Events(isnan(Cell_Events))=[];
Cell_Events=reshape(Cell_Events,[],3);

%%
% Process Noldus Behavioral Tracking Data
% Generate a filter: smoothing the speed of the mouse by convolution with a triangular filter 
smooth_win = triang(round(ds.Noldus_Sampling_Frequency * ap.Behavior_Filter_Size));
smooth_win = smooth_win / sum(smooth_win);
% Calculate the displacement of the mouse between each two adjust
% timestamps for session1
Displacement=zeros(size(Behavior_Tracking,1),1);
for n=1:size(Behavior_Tracking,1)
    if n+1<size(Behavior_Tracking,1)
        Displacement(n,1)=sqrt((Behavior_Tracking(n+1,2)-Behavior_Tracking(n,2))^2+(Behavior_Tracking(n+1,3)-Behavior_Tracking(n,3))^2); % calculate displacement at each 0.04sec interval
    end
end
% Fix potential errors: sometimes displacement calculated is nan because of very occassinoal
% Noldus tracking errors
Displacement(isnan(Displacement))=0;
Displacement(isinf(Displacement))=0;
% calculate the speed during each 0.04 time interval
Speed=Displacement./(1/ds.Noldus_Sampling_Frequency); 
% Smooth the speed of animal during session1
Speed = conv(Speed, smooth_win, 'same');

% Finishing up this section by concatenating speed info of the animal with
% the original behavioral tracking matrices
Behavior_Tracking=cat(2, Behavior_Tracking, Speed); 

%%
% Compute the XY-positions of the Mouse Given a Calcium Event Timestamp & Exclude Events During Immobility
% Find position of mouse for each calcium event during session1
% Notes about the organization of Cell_Events1&2 matrices after this step:
% Column 1 = timestamps of calcium event
% Column 2 = cell identity
% Column 3 = event amplitude
% Column 4&5 = X&Y position of the mouse at the event
% Column 6 = speed of the mouse at the event
[Cell_Events]=Find_XY_Position_Single_Session(Cell_Events, Behavior_Tracking, ds.Noldus_Sampling_Frequency);
% Exclude events during immobility
for n=1:size(Cell_Events,1)
    if Cell_Events(n,6)>ap.High_Mobility_Threshold || Cell_Events(n,6)<ap.Low_Mobility_Threshold
        Cell_Events(n,1:6)=nan(1,6);
    end
end
% Clean up
Cell_Events(any(isnan(Cell_Events),2),:) = []; 

end