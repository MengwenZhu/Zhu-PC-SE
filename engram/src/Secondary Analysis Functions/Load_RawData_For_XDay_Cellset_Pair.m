function [Date, Inscopix, Noldus, Start_T, Dose, Context] = Load_RawData_For_XDay_Cellset_Pair(META, ds)
% This function loads inscopix event csv file and corresponding noldus file
% from PM/day2 sessions for XDay analysis
% Also reports the starting times of recordings

% About input: META should be one row of the overall metadata excel sheet for
% secondary analysis

InscopixFilePaths = getInscopixRawDataPathFromMetadata(META, ds.CNMFe_Inscopix_rawdataPath);
NoldusFilePaths = getNoldusRawDataPathFromMetadata(META, ds.CNMFe_Inscopix_rawdataPath);

if strcmp(META.exptParadigm, '4h_memory_test')
    Date = META.exptDate1; 
    Inscopix = importRawData(InscopixFilePaths.PM);
    Noldus = importXYData(NoldusFilePaths.PM);
    Start_T = META.PM_sessionStart_CaImg_posix;
    Dose = META.session2_dose;
    Context = META.session2_context;
elseif strcmp(META.exptParadigm, '24h_memory_test')
    Date = META.exptDate2; 
    Inscopix = importRawData(InscopixFilePaths.day2);
    Noldus = importXYData(NoldusFilePaths.day2);
    Start_T = META.Day2_sessionStart_CaImg_posix;
    Dose = META.session2_dose;
    Context = META.session2_context;
end


end