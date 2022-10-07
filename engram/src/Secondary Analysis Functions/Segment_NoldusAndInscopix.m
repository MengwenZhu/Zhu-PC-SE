function [Segmented_Noldus, Segmented_Inscopix] = Segment_NoldusAndInscopix(Noldus, Inscopix, Segment_Start_T, Segment_End_T)

% Compute the segment of Noldus data you want
Segmented_Noldus = nan(size(Noldus,1), size(Noldus,2));
for n = 1:size(Noldus, 1)
    if Noldus(n,1)>=Segment_Start_T && Noldus(n,1)<=Segment_End_T
        Segmented_Noldus(n,:) = Noldus(n,:);
    end
end
% Clean up
for n = 1:size(Segmented_Noldus,1)
    if isnan(Segmented_Noldus(n,1))
        Segmented_Noldus(n,:) = 999999;
    end
end
Segmented_Noldus(Segmented_Noldus==999999)=[];
Segmented_Noldus = reshape(Segmented_Noldus, [], size(Noldus,2));

% Compute the segment of Inscopix data you want
Segmented_Inscopix = nan(size(Inscopix,1), size(Inscopix,2));
for n = 1:size(Inscopix, 1)
    if Inscopix(n,1)>=Segment_Start_T && Inscopix(n,1)<=Segment_End_T
        Segmented_Inscopix(n,:) = Inscopix(n,:);
    end
end
% Clean up
for n = 1:size(Segmented_Inscopix,1)
    if isnan(Segmented_Inscopix(n,1))
        Segmented_Inscopix(n,:) = 999999;
    end
end
Segmented_Inscopix(Segmented_Inscopix==999999)=[];
Segmented_Inscopix = reshape(Segmented_Inscopix, [], size(Inscopix,2));

end