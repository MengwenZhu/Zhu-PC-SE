function Comprehensive_LR(CNMFe_directory)
% This function runs comprehensive longitudinal registration for all
% denoised_cnmfe_cellsets from one animal, and then detect & export events
% into multiple excel csv files for secondary analysis

% Each recording with get its own LR file, own ED file, and own event
% export csv file, where each cell is given a globally registered name

% Input CNMFe_directory must be character, cannot be string

% If you really need to change analysis parameters here, you can modify the
% function later. Right now there is no need to modify anything here.

%% set directories
data_dir = CNMFe_directory; % have to define the CNMFe File folder for each animal here to run everything below 
output_dir = data_dir;
cd(output_dir)

%% remove previous files if you re-run the analysis 
% especially when you add new data for same animal, you must re-run LR to register all of them
if exist(fullfile(data_dir, 'LR export events'),'dir')
    rmdir(fullfile(data_dir, 'LR export events'), 's') % remove the original exported events directory
end

A = dir(data_dir);
for ii = 1:height(A)
    if A(ii).bytes > 0 && contains(A(ii).name,'LR')
        delete(A(ii).name);
    elseif A(ii).bytes > 0 && contains(A(ii).name,'ED')
        delete(A(ii).name);
    end
end

% Preparation
proc_cs_files = {};
for ii = 1:height(A)
    if A(ii).bytes > 0 && ~contains(A(ii).name,'LR') && ~contains(A(ii).name,'ED') % only select those that are denoised_cellset, in case you re-run the function
        proc_cs_files(ii) = cellstr(A(ii).name);
    end
end
proc_cs_files = proc_cs_files(~cellfun('isempty',proc_cs_files));

%append the file name to match the cnmfe cellset
LR_input_files = cellfun(@(x) fullfile(data_dir, x), proc_cs_files, 'UniformOutput', false);
LR_output_files = isx.make_output_file_paths(LR_input_files, output_dir, 'LR');

%% Run longitudinal registration function for denoised cnmfe_cellsets within a specific CNMFe folder

% run the following if there is no denoised_cellset_LR file
if ~exist(string(LR_output_files(1,1)), 'file')
    isx.longitudinal_registration( ...
         LR_input_files, LR_output_files, 'min_correlation', 0.5, 'accepted_cells_only', 1);
end

%% Run event detection
ED_input_files = LR_output_files;
ED_output_files = isx.make_output_file_paths(LR_input_files, output_dir, 'ED');

% run the following if there is no denoised_cellset_ED file
if ~exist(string(ED_output_files(1,1)), 'file')
    isx.event_detection(ED_input_files, ED_output_files, 'threshold', 4, 'tau', 0.2, 'event_time_ref', 'beginning', ...
        'ignore_negative_transients', true);
end

%% Export detected events into csv files for secondary analysis 
% Previously this is done by exporting all detected events into a
% comprehensive LR csv file, but it takes infinitely long for large dataset

% EE_input_files = ED_output_files;
% EE_output_file = fullfile(data_dir,'Comprehensive_LR_Event_Summary.csv');
% isx.export_event_set_to_csv(EE_input_files, EE_output_file, 'time_ref', 'unix', 'sparse_output', false, 'write_amplitude', true);

% Alternative strategy: export events for each recording, cells are named
% using their global name so that you know which cell is which across
% recordings

% create folder that contains all the exported csv files
if ~exist(string(fullfile(data_dir, 'LR export events')),'dir')
    mkdir('LR export events');
end

% make output directories
Event_Export_Directory = fullfile(CNMFe_directory, 'LR export events');
EE_input_files = ED_output_files;

csv_file_names = cell(1, size(proc_cs_files,2));
for kk = 1:size(proc_cs_files,2)
    csv_file_names(1,kk) = replace(proc_cs_files(1,kk), 'denoised_cellset.isxd', 'events.csv');
end

EE_output_files = cell(1, size(proc_cs_files,2));
for jj = 1:size(proc_cs_files,2)
    EE_output_files(1,jj) = fullfile(Event_Export_Directory, csv_file_names(1,jj));
end

% loop through all the ED files and export them as csv files
for qq = 1:size(EE_input_files, 2)
    if ~exist(string(fullfile(Event_Export_Directory, csv_file_names(1,qq))), 'file')
        isx.export_event_set_to_csv(EE_input_files(1,qq), char(EE_output_files(1,qq)), 'time_ref', 'unix', 'sparse_output', false, 'write_amplitude', true);
    end
end

end