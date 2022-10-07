function [Shuffled_PV_Corr, Shuffled_PV_Mean]=Shuffle_PV(Num_Place_Cell, BothSession_Smoothed_Ca_Rate_Map1, BothSession_Smoothed_Ca_Rate_Map2, PlaceCell_Summary1, PlaceCell_Summary2, ap)
%% Shuffle PV correlation between smoothed rate maps by shuffling the identities of pixels for each spatial location

% fix some hardware coding parameters
size1=size(BothSession_Smoothed_Ca_Rate_Map1,1);
size2=size(BothSession_Smoothed_Ca_Rate_Map1,2);

% Use place cells if there are enought place cells; if aren't enough, use
% all cells instead
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

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% in this case, need to get rid of intervening zero matrices because they
% will affect the PV correlation values, so just keep cells that are real 
for n=1:size(BothSession_Smoothed_Ca_Rate_Map1,3)
    if BothSession_Smoothed_Ca_Rate_Map1(:,:,n)==zeros(size(BothSession_Smoothed_Ca_Rate_Map1,1),size(BothSession_Smoothed_Ca_Rate_Map1,2))
        BothSession_Smoothed_Ca_Rate_Map1(:,:,n)=nan;
    end
end
BothSession_Smoothed_Ca_Rate_Map1(isnan(BothSession_Smoothed_Ca_Rate_Map1))=[];
BothSession_Smoothed_Ca_Rate_Map1=reshape(BothSession_Smoothed_Ca_Rate_Map1, size1, size2, []);



for n=1:size(BothSession_Smoothed_Ca_Rate_Map2,3)
    if BothSession_Smoothed_Ca_Rate_Map2(:,:,n)==zeros(size(BothSession_Smoothed_Ca_Rate_Map2,1),size(BothSession_Smoothed_Ca_Rate_Map2,2))
        BothSession_Smoothed_Ca_Rate_Map2(:,:,n)=nan;
    end
end
BothSession_Smoothed_Ca_Rate_Map2(isnan(BothSession_Smoothed_Ca_Rate_Map2))=[];
BothSession_Smoothed_Ca_Rate_Map2=reshape(BothSession_Smoothed_Ca_Rate_Map2, size1, size2, []);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


Shuffled_Smoothed_Ca_Rate_Map1=zeros(size(BothSession_Smoothed_Ca_Rate_Map1,1),size(BothSession_Smoothed_Ca_Rate_Map1,2),size(BothSession_Smoothed_Ca_Rate_Map1,3));
for n=1:size(BothSession_Smoothed_Ca_Rate_Map1,1)
    for m=1:size(BothSession_Smoothed_Ca_Rate_Map1,2)
        Shuffled_Smoothed_Ca_Rate_Map1(n,m,:)=reshape(Matrix_Perm_Shuffling(BothSession_Smoothed_Ca_Rate_Map1(n,m,:)),1,1,[]);
    end
end

Shuffled_Smoothed_Ca_Rate_Map2=zeros(size(BothSession_Smoothed_Ca_Rate_Map2,1),size(BothSession_Smoothed_Ca_Rate_Map2,2),size(BothSession_Smoothed_Ca_Rate_Map2,3));
for n=1:size(BothSession_Smoothed_Ca_Rate_Map2,1)
    for m=1:size(BothSession_Smoothed_Ca_Rate_Map2,2)
        Shuffled_Smoothed_Ca_Rate_Map2(n,m,:)=reshape(Matrix_Perm_Shuffling(BothSession_Smoothed_Ca_Rate_Map2(n,m,:)),1,1,[]);
    end
end



if size(Shuffled_Smoothed_Ca_Rate_Map1,3)>size(Shuffled_Smoothed_Ca_Rate_Map2,3)
    Shuffled_Smoothed_Ca_Rate_Map2(1:size(Shuffled_Smoothed_Ca_Rate_Map2,1),1:size(Shuffled_Smoothed_Ca_Rate_Map2,2),end+1:size(Shuffled_Smoothed_Ca_Rate_Map1,3))=zeros(size(Shuffled_Smoothed_Ca_Rate_Map2,1),size(Shuffled_Smoothed_Ca_Rate_Map2,2),(size(Shuffled_Smoothed_Ca_Rate_Map1,3)-size(Shuffled_Smoothed_Ca_Rate_Map2,3)));
end

if size(Shuffled_Smoothed_Ca_Rate_Map1,3)<size(Shuffled_Smoothed_Ca_Rate_Map2,3)
    Shuffled_Smoothed_Ca_Rate_Map1(1:size(Shuffled_Smoothed_Ca_Rate_Map1,1),1:size(Shuffled_Smoothed_Ca_Rate_Map1,2),end+1:size(Shuffled_Smoothed_Ca_Rate_Map2,3))=zeros(size(Shuffled_Smoothed_Ca_Rate_Map1,1),size(Shuffled_Smoothed_Ca_Rate_Map1,2),(size(Shuffled_Smoothed_Ca_Rate_Map2,3)-size(Shuffled_Smoothed_Ca_Rate_Map1,3)));
end



Shuffled_PV_Corr=PV_Correlation(Shuffled_Smoothed_Ca_Rate_Map1,Shuffled_Smoothed_Ca_Rate_Map2);
Shuffled_PV_Mean=mean(Shuffled_PV_Corr, 'all', 'omitnan');

end