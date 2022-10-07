function [PlaceCell_Summary1, PlaceCell_Summary2] = Find_Placeness_Summary_For_Cellset_Pair(META, Placeness_MetaData)

% This function generates a pair of matrices that reflect place cell
% information for a session pair to be analyzed in secondary analysis

% metadata is a row of secondary analysis metadata excel sheet

% Placeness_MetaData is a structure that locates within a mat file in
% metadata directory, which must be loaded first


if strcmp(META.exptParadigm, '4h_memory_test')
    for n = 1:size(Placeness_MetaData,2)
        if strcmp(Placeness_MetaData(n).animalName, META.animalName) && strcmp(Placeness_MetaData(n).exptDate1, META.exptDate1) && strcmp(Placeness_MetaData(n).session, 'AM')
            PlaceCell_Summary1 = Placeness_MetaData(n).Place_Cell_Summary;
            PlaceCell_Summary2 = Placeness_MetaData(n+1).Place_Cell_Summary;
        end
    end
elseif strcmp(META.exptParadigm, '24h_memory_test')
    for n = 1:size(Placeness_MetaData,2)
        if strcmp(Placeness_MetaData(n).animalName, META.animalName) && strcmp(Placeness_MetaData(n).exptDate1, META.exptDate1) && strcmp(Placeness_MetaData(n).exptDate2, META.exptDate2) && strcmp(Placeness_MetaData(n).session, 'day1')
            PlaceCell_Summary1 = Placeness_MetaData(n).Place_Cell_Summary;
            PlaceCell_Summary2 = Placeness_MetaData(n+1).Place_Cell_Summary;
        end
    end
end



end