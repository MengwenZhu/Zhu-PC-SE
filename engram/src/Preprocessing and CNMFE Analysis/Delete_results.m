% this scrip tries to delete previously-generated primary analysis files so
% that you could re-run the primary analysis. Don't run this if you want to
% keep the primary analysis results (non-reversible).

% get all analysis parameters & file paths
[ds, ap] = getdsap;

% extract metadata excel sheet that contains all recording information 
metadata = readtable(fullfile(ds.metadataPath, ds.metadataFileName_forPrimaryAnalysis));

% start looping through and deleting files
numRows = height(metadata);
for ii = 1:numRows
    % get paths and names to recordings
    [CNMFefilePaths, filePaths, fileNames] = getCaImagingRawRecordingPathAndNameFromMetadata(table2struct(metadata(ii, :)), ...
        ds.CNMFe_Inscopix_rawdataPath);
    % get a row of metadata table
    META = table2struct(metadata(ii,:));
    % if CNMFe Files folder exists, delete it
    if exist(CNMFefilePaths.path, 'dir')
        rmdir(CNMFefilePaths.path, 's')
    end
    % depending on the type of experimental paradigm the folders have
    % different names
    if strcmp(META.exptParadigm, '24h_memory_test')
        % get the list of file info within target folder
        files1 = dir(filePaths.day1);
        files2 = dir(filePaths.day2);
        % extract file names within target folder
        filenames1 = {files1.name};
        filenames2 = {files2.name};
        % remove unnecessary names
        filenames1=filenames1(~ismember(filenames1,{'.','..'}));
        filenames2=filenames2(~ismember(filenames2,{'.','..'}));
        % try to delete all files except the original triggered recording
        for n = 1:size(filenames1,2)
            % to remove files or folders you need different strategies
            if ~contains(filenames1(1,n), 'video_trig') && exist(fullfile(filePaths.day1, filenames1(1,n))) == 2
                delete(fullfile(filePaths.day1, filenames1(1,n)));
            elseif ~contains(filenames1(1,n), 'video_trig') && exist(fullfile(filePaths.day1, filenames1(1,n))) == 7
                rmdir(fullfile(filePaths.day1, filenames1(1,n)), 's');
            end
        end
        for n = 1:size(filenames2,2)
            if ~contains(filenames2(1,n), 'video_trig') && exist(fullfile(filePaths.day2, filenames2(1,n))) == 2
                delete(fullfile(filePaths.day2, filenames2(1,n)));
            elseif ~contains(filenames2(1,n), 'video_trig') && exist(fullfile(filePaths.day2, filenames2(1,n))) == 7
                rmdir(fullfile(filePaths.day2, filenames2(1,n)), 's');
            end
        end
    elseif strcmp(META.exptParadigm, '4h_memory_test')
        % get the list of file info within target folder
        files1 = dir(filePaths.AM);
        files2 = dir(filePaths.PM);
        % extract file names within target folder
        filenames1 = {files1.name};
        filenames2 = {files2.name};
        % remove unnecessary names
        filenames1=filenames1(~ismember(filenames1,{'.','..'}));
        filenames2=filenames2(~ismember(filenames2,{'.','..'}));
        % try to delete all files except the original triggered recording
        for n = 1:size(filenames1,2)
            % to remove files or folders you need different strategies
            if ~contains(filenames1(1,n), 'video_trig') && exist(fullfile(filePaths.AM, filenames1(1,n))) == 2
                delete(fullfile(filePaths.AM, filenames1(1,n)));
            elseif ~contains(filenames1(1,n), 'video_trig') && exist(fullfile(filePaths.AM, filenames1(1,n))) == 7
                rmdir(fullfile(filePaths.AM, filenames1(1,n)), 's');
            end
        end
        for n = 1:size(filenames2,2)
            if ~contains(filenames2(1,n), 'video_trig') && exist(fullfile(filePaths.PM, filenames2(1,n))) == 2
                delete(fullfile(filePaths.PM, filenames2(1,n)));
            elseif ~contains(filenames2(1,n), 'video_trig') && exist(fullfile(filePaths.PM, filenames2(1,n))) == 7
                rmdir(fullfile(filePaths.PM, filenames2(1,n)), 's');
            end
        end
    end
end