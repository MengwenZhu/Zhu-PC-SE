function [Data_Summary1, Data_Summary2] = runPlacenessAnalysis(Inscopix1, Noldus1, Start_T1, Inscopix2, Noldus2, Start_T2, META, ds, ap)
% Just a function that combines a pair of MI shuffle tests into a single
% function

% perform place cell analysis for earlier one of the recording pair (AM/day1)
Data_Summary1 = Mutual_Info_Shuffle_Test_For_One_Cellset(Inscopix1, Noldus1, Start_T1, META, ds, ap);
% perform place cell analysis for later one of the recording pair (PM/day2)
Data_Summary2 = Mutual_Info_Shuffle_Test_For_One_Cellset(Inscopix2, Noldus2, Start_T2, META, ds, ap);

% organize the output strcutures into a structure that could be
% referred to in later secondary analysis
if strcmp(META.exptParadigm, '4h_memory_test')
    Data_Summary1.session = 'AM';
    Data_Summary2.session = 'PM';
elseif strcmp(META.exptParadigm, '24h_memory_test')
    Data_Summary1.session = 'day1';
    Data_Summary2.session = 'day2';
end

end