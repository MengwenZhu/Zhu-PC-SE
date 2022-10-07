function [Inscopix1, Inscopix2, Noldus1, Noldus2, Start_T1, Start_T2] = Load_RawData_For_Cellset_Pair(META, ds)
% This function loads inscopix event csv file and corresponding noldus file as two numeric matrices 
% Also reports the starting times of the pair of recordings

% About input: META should be one row of the overall metadata excel sheet for
% secondary analysis

% Generate directories to an AM/PM or Day1/Day2 pair of Inscopix and Noldus files
% Field names either AM/PM or day1/day2, depending on the exptParadigm

InscopixFilePaths = getInscopixRawDataPathFromMetadata(META, ds.CNMFe_Inscopix_rawdataPath);
NoldusFilePaths = getNoldusRawDataPathFromMetadata(META, ds.CNMFe_Inscopix_rawdataPath);

if strcmp(META.exptParadigm, '4h_memory_test')
    Inscopix1 = importRawData(InscopixFilePaths.AM);
    Inscopix2 = importRawData(InscopixFilePaths.PM);
    Noldus1 = importXYData(NoldusFilePaths.AM);
    Noldus2 = importXYData(NoldusFilePaths.PM);
    Start_T1 = META.AM_sessionStart_CaImg_posix;
    Start_T2 = META.PM_sessionStart_CaImg_posix;
elseif strcmp(META.exptParadigm, '24h_memory_test')
    Inscopix1 = importRawData(InscopixFilePaths.day1);
    Inscopix2 = importRawData(InscopixFilePaths.day2);
    Noldus1 = importXYData(NoldusFilePaths.day1);
    Noldus2 = importXYData(NoldusFilePaths.day2);
    Start_T1 = META.Day1_sessionStart_CaImg_posix;
    Start_T2 = META.Day2_sessionStart_CaImg_posix;
end


end
