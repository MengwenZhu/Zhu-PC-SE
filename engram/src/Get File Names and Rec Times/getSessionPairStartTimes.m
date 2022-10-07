function getSessionPairStartTimes(ds, ap)

% getSessionPairStartTimes(ds, ap) loops over the rows (=AM/PM or day1/day2 session pairs) of the
% metadata file as specified in struct ds ('data set'), tries to locate the
% denoised_cellset files of the session pair, and determines the session start times.
% The times are converted to posix time (number of seconds (including
% fractional seconds) elapsed since 00:00:00 1-Jan-1970 UTC). All time info is
% appended to the metadata and written to a new metadata file (for secondary analysis) in the
% original file's directory.
% Input args ds ('data set') and ap ('analysis parameters') can be loaded
% via
%       [ds, ap] = getdsap;

% The current secondary analysis workflow will not require a computation of Tcorrection
% because inscopix input raw data will be all exported using posix time,
% so we just need to subtract start posix times from each timestamp to
% correct for that specific event set

% read metadata
metadata = readtable(fullfile(ds.metadataPath, ds.metadataFileName_forPrimaryAnalysis));
numRows = height(metadata);


% loop over entries (experiments)
for ii = 1:numRows
    % construct and verify full paths to CNMFe raw data files
    filePaths = getCNMFeRawDataFilePathFromMetadata(table2struct(metadata(ii, :)), ...
        ds.CNMFe_Inscopix_rawdataPath);
    if strcmp(metadata.exptParadigm{ii}, '4h_memory_test')
        % add new columns to metadata or overwrite existing ones
        %         metadata.sessionStart_CaImg_AM = NaT(numRows, 1);
        %         metadata.sessionStart_CaImg_PM = NaT(numRows, 1);
        tod = {'AM', 'PM'};
        % verify that struct filePaths has the requested fields
        if ~isequal(fieldnames(filePaths)', tod)
            error(['expecting fieldnames ', char(join(tod, ', '))])
        end
        % loop over tods
        for ix = 1:length(tod)
            if ~ismissing(filePaths.(tod{ix}))
                sessionTime = getCNMFeSessionTime(filePaths.(tod{ix}), ap);
                metadata{ii, strcat('sessionStart_CaImg_', tod(ix))} = sessionTime.start;
            end
        end
        % convert session times to posixtime (unit is s) and add to table
        metadata.AM_sessionStart_CaImg_posix = posixtime(metadata.sessionStart_CaImg_AM);
        metadata.PM_sessionStart_CaImg_posix = posixtime(metadata.sessionStart_CaImg_PM);
    elseif strcmp(metadata.exptParadigm{ii}, '24h_memory_test')
        % add new columns to metadata or overwrite existing ones
        %         metadata.sessionStart_CaImg_day1 = NaT(numRows, 1);
        %         metadata.sessionStart_CaImg_day2 = NaT(numRows, 1);
        tod = {'day1', 'day2'};
        % verify that struct filePaths has the requested fields
        if ~isequal(fieldnames(filePaths)', tod)
            error(['expecting fieldnames ', char(join(tod, ', '))])
        end
        % loop over tods
        for ix = 1:length(tod)
            if ~ismissing(filePaths.(tod{ix}))
                sessionTime = getCNMFeSessionTime(filePaths.(tod{ix}), ap);
                metadata{ii, strcat('sessionStart_CaImg_', tod(ix))} = sessionTime.start;
            end
        end
        % convert session times to posixtime (unit is s) and add to table
        metadata.Day1_sessionStart_CaImg_posix = posixtime(metadata.sessionStart_CaImg_day1);
        metadata.Day2_sessionStart_CaImg_posix = posixtime(metadata.sessionStart_CaImg_day2);
    end
end

% write amended table to file with ' with start times' appended to file name
[~, name, ext] = fileparts(ds.metadataFileName_forPrimaryAnalysis);
name = strrep(name, 'Primary', 'Secondary');
newFn = fullfile(ds.metadataPath, [name, ext]);
disp(['writing new metadata table to ', newFn]);
writetable(metadata, newFn);
