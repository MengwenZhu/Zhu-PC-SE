% script that automatically run through primary analysis for all recordings
% listed on metadata excel sheet

% At first you will get a denoised_cnmfe_cellset for each recording,
% located within CNMFe Files in each animal's folder, and then these CNMFe
% files will be longitudinally registered, event detected, and the result
% will be exported as a single csv file, which will be used for secondary
% batch analysis.

% currently CNMFe is already running in parallel mode, so there is no
% obvious way to speed up this part

% P.S. if anything runs too slowly, it most likely is due to bad CNMFe
% analysis parameters within the excel metadata sheet. Try to manually
% investigate the CNMFe cellsets within IDPS and re-optimize your analysis
% parameters

%% Preparation

% get all analysis parameters & file paths
[ds, ap] = getdsap;

% extract metadata excel sheet that contains all recording information 
metadata = readtable(fullfile(ds.metadataPath, ds.metadataFileName_forPrimaryAnalysis));

%% Run through all recordings and get denoised_cnmfe_cellset

% start timing the analysis
t_start = tic;

% loop through all the recordings (can't parallelize this because the 
% processing of one recording by CNMFe is already parallelized)
numRows = height(metadata);
for ii = 1:numRows
    % get paths and names to recordings
    [CNMFefilePaths, filePaths, fileNames] = getCaImagingRawRecordingPathAndNameFromMetadata(table2struct(metadata(ii, :)), ...
        ds.CNMFe_Inscopix_rawdataPath);
        Meta_Data = table2struct(metadata(ii, :));
        Info_Table = struct();
        % right now there are two experimental paradigms: recording
        % separated by 4h (AM/PM pair) & those separated by 24h
        % (day1/day2 pair)
        if strcmp(metadata.exptParadigm{ii}, '4h_memory_test') 
            tod = {'AM', 'PM'};
            for ix = 1:length(tod)
                Info_Table.Recording_Directory = filePaths.(string(tod(ix)));
                Info_Table.Recording_Name = fileNames.(string(tod(ix)));
                Info_Table.CNMFe_Directory = CNMFefilePaths.path;
                % get primary analysis parameters that are specific to each
                % mouse (you find those values by manually testing out parameters within IDPS 
                % and then enter into metadata excel sheet)
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                Info_Table.spat_ds = Meta_Data.spatial_ds;
                Info_Table.temp_ds = Meta_Data.temporal_ds;
                Info_Table.cell_diameter = Meta_Data.cell_diameter;
                Info_Table.min_corr = Meta_Data.min_corr; 
                Info_Table.min_pnr = Meta_Data.min_pnr; 
                Info_Table.spike_snr_threshold = Meta_Data.spike_snr_threshold;
                Info_Table.noise_range_min = Meta_Data.noise_range_min;
                Info_Table.noise_range_max = Meta_Data.noise_range_max;
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                % assign each denoised cellset a unique name before you generate it
                if strcmp(string(tod(ix)), 'AM')
                Info_Table.Denoised_CNMFe_Name = strcat(convertCharsToStrings(Meta_Data.exptDate1), {' '}, string(tod(ix)), {' '}, 'denoised_cellset.isxd');
                elseif strcmp(string(tod(ix)), 'PM')
                Info_Table.Denoised_CNMFe_Name = strcat(convertCharsToStrings(Meta_Data.exptDate1), {' '}, string(tod(ix)), {' '}, 'denoised_cellset.isxd');  
                end
                % run analysis if and only if denoised cellset hasn't been generated
                if ~isfile(fullfile(Info_Table.CNMFe_Directory,Info_Table.Denoised_CNMFe_Name))
                   Primary_Analysis_V1(Info_Table);
                end
           end
        elseif strcmp(metadata.exptParadigm{ii}, '24h_memory_test')
                tod = {'day1', 'day2'};
                for ix = 1:length(tod)
                Info_Table.Recording_Directory = filePaths.(string(tod(ix)));
                Info_Table.Recording_Name = fileNames.(string(tod(ix)));
                Info_Table.CNMFe_Directory = CNMFefilePaths.path; 
                % get primary analysis parameters that are specific to each
                % mouse
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                Info_Table.spat_ds = Meta_Data.spatial_ds;
                Info_Table.temp_ds = Meta_Data.temporal_ds;
                Info_Table.cell_diameter = Meta_Data.cell_diameter;
                Info_Table.min_corr = Meta_Data.min_corr; 
                Info_Table.min_pnr = Meta_Data.min_pnr;                
                Info_Table.spike_snr_threshold = Meta_Data.spike_snr_threshold;
                Info_Table.noise_range_min = Meta_Data.noise_range_min;
                Info_Table.noise_range_max = Meta_Data.noise_range_max;
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                % name cellsets
                if strcmp(string(tod(ix)), 'day1')
                Info_Table.Denoised_CNMFe_Name = strcat(convertCharsToStrings(Meta_Data.exptDate1), {' '}, string(tod(ix)), {' '}, 'denoised_cellset.isxd');
                elseif strcmp(string(tod(ix)), 'day2')
                Info_Table.Denoised_CNMFe_Name = strcat(convertCharsToStrings(Meta_Data.exptDate2), {' '}, string(tod(ix)), {' '}, 'denoised_cellset.isxd');  
                end
                % run analysis if denoised cellset hasn't been generated
                if ~isfile(fullfile(Info_Table.CNMFe_Directory,Info_Table.Denoised_CNMFe_Name))
                   Primary_Analysis_V1(Info_Table);
                end
                end
        end
end

%% LR all denoised_cnmfe_cellset for each animal, detect events, and export csv files for secondary analysis

% first get directories to each animal's CNMFe folder
LR_Paths = struct();
for jj = 1:numRows
     [LR_Paths(jj).path, ~, ~] = getCaImagingRawRecordingPathAndNameFromMetadata(table2struct(metadata(jj,:)), ds.CNMFe_Inscopix_rawdataPath);
end

for kk = 1:size(LR_Paths,2)
    Update_Path(kk,1) = LR_Paths(kk).path;
end

% get unique values within the path structure
[~, idx] = unique([Update_Path.path].', 'rows', 'stable');
Update_Path = Update_Path(idx);

% input the directory to get automated longitudinal registration
for n = 1:size(Update_Path,1)
    Comprehensive_LR(convertStringsToChars(Update_Path(n).path));
end

% stop timing of the whole analysis and report time elapsed
t_stop = toc(t_start);
disp("Primary analysis took " + string(t_stop) + " s")

% clear all variables
clear