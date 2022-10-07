function [Shuffled_Rate_Map_Corr]=Shuffle_Smoothed_Rate_Map_Correlation(Num_Place_Cell, BothSession_Smoothed_Ca_Rate_Map1, BothSession_Smoothed_Ca_Rate_Map2, PlaceCell_Summary1, PlaceCell_Summary2, ap)
%% Shuffle the correlation between smoothed rate maps by shuffling the identities of the bins

% Reject cells that are non-place cell during both sessions if and only if
% the number of place cells observed in both sessions exceeds a certain
% number
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

% this shuffle requires the 'Matrix_Perm_Shuffling' function, which
% randomly permutate values in a matrix
Shuffled_Smoothed_Ca_Rate_Map1=zeros(1,size(BothSession_Smoothed_Ca_Rate_Map1,1)*size(BothSession_Smoothed_Ca_Rate_Map1,2),size(BothSession_Smoothed_Ca_Rate_Map1,3));
for n=1:size(BothSession_Smoothed_Ca_Rate_Map1,3)
    Shuffled_Smoothed_Ca_Rate_Map1(1,:,n)=Matrix_Perm_Shuffling(BothSession_Smoothed_Ca_Rate_Map1(:,:,n)); % (AM) shuffle each Ca rate matrix by randomly permutating values within that matrix
end

Shuffled_Smoothed_Ca_Rate_Map2=zeros(1,size(BothSession_Smoothed_Ca_Rate_Map2,1)*size(BothSession_Smoothed_Ca_Rate_Map2,2),size(BothSession_Smoothed_Ca_Rate_Map2,3));
for n=1:size(BothSession_Smoothed_Ca_Rate_Map2,3)
    Shuffled_Smoothed_Ca_Rate_Map2(1,:,n)=Matrix_Perm_Shuffling(BothSession_Smoothed_Ca_Rate_Map2(:,:,n)); % same thing for PM
end

Shuffled_Rate_Map_Corr=nan(size(BothSession_Smoothed_Ca_Rate_Map1,3),1);
for n=1:size(BothSession_Smoothed_Ca_Rate_Map1,3)
    Shuffled_Rate_Map_Corr(n,1)=XCorrCalc(Shuffled_Smoothed_Ca_Rate_Map1(1,:,n),Shuffled_Smoothed_Ca_Rate_Map2(1,:,n)); % calculate the rate map correlations between shuffled rate maps to get null distribution
end


end