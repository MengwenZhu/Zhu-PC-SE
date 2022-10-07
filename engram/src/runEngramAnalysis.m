% Main script for running Engram Cell Analysis based on the
% exampleAnimalTable spreadsheet

% load ds ('data set') and ap ('analysis parameters')
[ds, ap] = getdsap;

% read metadata
Info_Table = readtable(fullfile(ds.metadataPath, ds.metadataFileName_forSecondaryAnalysis));

% Load structure that contains place cell information for all animals
load(fullfile(ds.metadataPath, 'Placeness_MetaData.mat'));

% loop through table to run analysis on each date/treatment (loop backwards
% for automatic preallocation)
Data_Summary = [];

% start the clock after the prelude above
t_start = tic;
% for ii = height(Info_Table):-1:1
for ii = 1:height(Info_Table)
    META = table2struct(Info_Table(ii,:));
    animalName = Info_Table.animalName(ii);
    genotype = Info_Table.genotype(ii);
    drug = Info_Table.drug{ii};
    session1_dose = Info_Table.session1_dose(ii);
    session2_dose = Info_Table.session2_dose(ii);
    session1_context = Info_Table.session1_context(ii);
    session2_context = Info_Table.session2_context(ii);
    exptParadigm = Info_Table.exptParadigm(ii);
    % Have to be careful with which date to report 
    if strcmp(META.exptParadigm, '4h_memory_test')
        dateToFind = Info_Table.exptDate1{ii};
    elseif strcmp(META.exptParadigm, '24h_memory_test')
        dateToFind = string(strcat(Info_Table.exptDate1{ii}, {' '}, 'and', {' '}, Info_Table.exptDate2{ii}));
    end

    % try to load data
    try
        [Inscopix1, Inscopix2, Noldus1, Noldus2, Start_T1, Start_T2] = Load_RawData_For_Cellset_Pair(META, ds); 
        [PlaceCell_Summary1, PlaceCell_Summary2] = Find_Placeness_Summary_For_Cellset_Pair(META, Placeness_MetaData);
        isDataRead = true;
    catch
        isDataRead = false;
        formatSpec = "Retrieving imaging and tracking data for %s, %s failed";
        msg = compose(formatSpec, animalName{:}, dateToFind);
        warndlg(msg)
    end
    if isDataRead
        if ap.exec_type == "serial"
            % --- serial execution of analysis ---
            formatSpec = "Running engram analysis on %s, %s ...";
            disp(compose(formatSpec, animalName{:}, dateToFind));
            % main analysis function
            data = runEngram_Cell_Analysis(Inscopix1, Inscopix2, Noldus1, Noldus2, Start_T1, Start_T2, PlaceCell_Summary1, PlaceCell_Summary2, dateToFind, animalName, drug, session1_dose, session2_dose, session1_context, session2_context, genotype, exptParadigm, ds, ap);
            % collect results
            if isempty(Data_Summary)
                Data_Summary = data;
            elseif ~isempty(Data_Summary)
                Data_Summary = [Data_Summary, data];
            end
        elseif ap.exec_type == "parallel"
            % --- parallel execution of analysis ---
            futures(ii) = parfeval(@runEngram_Cell_Analysis, 1, Inscopix1, Inscopix2, Noldus1, Noldus2, Start_T1, Start_T2, PlaceCell_Summary1, PlaceCell_Summary2, dateToFind, animalName, drug, session1_dose, session2_dose, session1_context, session2_context, genotype, exptParadigm, ds, ap);
            formatSpec = "--- Submitted %s, %s (job ID %i) to the analysis queue";
            disp(compose(formatSpec, animalName{:}, dateToFind, futures(ii).ID));
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
            [~, data] = fetchNext(futures(jix));
            if isempty(Data_Summary)
                Data_Summary = data;
            elseif ~isempty(Data_Summary)
                Data_Summary = [Data_Summary, data];
            end
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
            animalName = args{10};
            dateToFind = args{9};
            formatSpec = "Analysis for %s, %s, failed";
            msg = compose(formatSpec, animalName{:}, dateToFind);
            warndlg(msg)
        end
    end
    % after all is done, cancel
    cancel(futures);
end

% save workspace
save(fullfile(ds.metadataPath, 'Data Summary.mat'), '-v7.3', 'Data_Summary')
% report
t_stop = toc(t_start);
disp("Secondary analysis took " + string(t_stop) + " s")
% clear all variables
clear