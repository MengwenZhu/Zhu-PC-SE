% This script aims to loop through all the inscopix and noldus raw data
% files to get information about the spatial firing characteristics of
% cells. Each animal should get its own structure that summarizes
% everything about their place cells

% get all basic parameters and directories
[ds, ap] = getdsap;

% read metadata for secondary analysis
metadata = readtable(fullfile(ds.metadataPath, ds.metadataFileName_forSecondaryAnalysis));
numRows = height(metadata);

% define a structure that contains all info about placeness of cells from
% all animals
Placeness_MetaData = [];

% time the analysis
t_start = tic;
for ii = 1:numRows
    META = table2struct(metadata(ii,:));
    try
        % load data
        [Inscopix1, Inscopix2, Noldus1, Noldus2, Start_T1, Start_T2] = Load_RawData_For_Cellset_Pair(META, ds);
        isDataRead = true;
    catch
        isDataRead = false;
        formatSpec = "Retrieving imaging and tracking data for %s, %s/%s failed";
        msg = compose(formatSpec, META.animalName, META.exptDate1, META.exptDate2);
        warndlg(msg)
    end
    if isDataRead     
        if ap.exec_type == "serial"
            % --- serial execution of analysis ---
            formatSpec = "Running placeness analysis on %s, %s/%s ...";
            disp(compose(formatSpec, META.animalName, META.exptDate1, META.exptDate2));
            % run analysis
            [Data_Summary1, Data_Summary2] = runPlacenessAnalysis(Inscopix1, Noldus1, Start_T1, Inscopix2, Noldus2, Start_T2, META, ds, ap);
            % append
            Placeness_MetaData = [Placeness_MetaData, Data_Summary1, Data_Summary2];
        elseif ap.exec_type == "parallel"
            % --- parallel execution of analysis ---
            futures(ii) = parfeval(@runPlacenessAnalysis, 2, Inscopix1, Noldus1, Start_T1, Inscopix2, Noldus2, Start_T2, META, ds, ap);
            formatSpec = "--- Submitted %s, %s/%s (job ID %i) to the analysis queue";
            disp(compose(formatSpec, META.animalName, ...
                META.exptDate1, META.exptDate2, futures(ii).ID));
        else
            error('bad value for ap.exec_type');
        end
    end
end

if ap.exec_type == "parallel"
    % Fetch results as they become available. Note that the array of
    % futures may contain entries with an 'unavailable' state, namely, the
    % jobs for which the data were faulty and which consequently have not
    % been parfeval'd; these must be omitted
    isJobNotFinished = true(1, length(futures));
    % preallocate index for canceled jobs
    isJobCanceled = false(1, length(futures));
    while any(isJobNotFinished)
        % identify jobs which have finished and not been read, so in
        % principle are ready to be fetched
        isJobToBeFetched = strcmp({futures.State}, 'finished') & ~[futures.Read];
        % of those, kick the ones having been canceled (which are also in
        % state 'finished' but have a nonempty Error property (which we
        % unfortunately have to query in a loop))
        for ix = 1:numel(isJobToBeFetched)
            if isJobToBeFetched(ix) && ~isempty(futures(ix).Error)
                isJobCanceled(ix) = true;
                isJobToBeFetched(ix) = false;
            end
        end
        % fetch results
        for jix = find(isJobToBeFetched)
            [~, Data_Summary1, Data_Summary2] = fetchNext(futures(jix));
            Placeness_MetaData = [Placeness_MetaData, Data_Summary1, Data_Summary2];
            formatSpec = "Successfully fetched data for job ID %i";
            disp(compose(formatSpec, futures(jix).ID));
        end
        % identify unfinished jobs - note that failed and unavailable jobs
        % are also defined as finished
        isJobNotFinished = strcmp({futures.State}, 'running') | strcmp({futures.State}, 'queued');
        % brief pause to prevent the while loop from consuming too many
        % resources
        pause(0.5)
    end
    % report jobs with errors (not unavailable ones, which have been
    % reported before)
    for f = futures
        if ~isempty(f.Error)
            args = f.InputArguments;
            meta = args{7};
            formatSpec = "Analysis for %s, %s/%s, failed";
            msg = compose(formatSpec, meta.animalName, meta.exptDate1, meta.exptDate2);
            warndlg(msg)
        end
    end
    % after all is done, cancel
    cancel(futures);
end

% fix bugs (have to sort placeness metadata)
[~,index] = sortrows({Placeness_MetaData.exptDate1}.'); Placeness_MetaData = Placeness_MetaData(index); clear index
[~,index] = sortrows({Placeness_MetaData.animalName}.'); Placeness_MetaData = Placeness_MetaData(index); clear index

% save the placeness structure within metadata directory and load it later
% for other analysis
save(fullfile(ds.metadataPath, 'Placeness_MetaData.mat'), '-v7.3', 'Placeness_MetaData')
% report
t_stop = toc(t_start);
disp("Place cell analysis took " + string(t_stop) + " s")
% clear all variables
clear